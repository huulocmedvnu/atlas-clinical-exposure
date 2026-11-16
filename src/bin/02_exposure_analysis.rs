//! Step 2. Exposure analysis of the Tahoe-100M summaries written by `tahoe_summary`: share of conditions and of DE calls above
//! clinical exposure, background of inert compounds, overlap between mechanisms, mechanism similarity, compound signatures.
//! Usage: exposure_analysis <summary dir> <data dir> <out dir> <results json>
//! Writes TSV tables to <out dir> and all test statistics to <results json>.
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use arrow::array::*;
use arrow::compute::cast;
use arrow::datatypes::DataType;
use arrow::record_batch::RecordBatch;
use atlas_exposure::*;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;
use rayon::prelude::*;
use serde_json::{json, Map, Value};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Mutex;

const DOSES: [f64; 3] = [0.05, 0.5, 5.0];
const CUTS: [f32; 5] = [0.25, 0.5, 0.75, 1.0, 1.5];
const CUTN: [&str; 5] = ["0.25", "0.5", "0.75", "1.0", "1.5"];
const BASIS: [&str; 3] = ["unbound", "total", "fbs"];
const SETS: [&str; 3] = ["liston", "label", "all"];
const STRATA: [(&str, u8); 4] = [("all", 9), ("straddling", 1), ("always_above", 2), ("all_reachable", 0)];

fn dose_idx(c: f32) -> u8 { if c < 0.2 { 0 } else if c < 2.0 { 1 } else { 2 } }

#[derive(Default)]
struct Intern { map: HashMap<String, u32>, names: Vec<String> }
impl Intern {
    fn get(&mut self, s: &str) -> u32 {
        if let Some(&i) = self.map.get(s) { return i; }
        let i = self.names.len() as u32; self.map.insert(s.to_string(), i); self.names.push(s.to_string()); i
    }
}
/// Per-thread cache in front of a shared interner.
struct Cache<'a> { g: &'a Mutex<Intern>, l: HashMap<String, u32>, last: (String, u32) }
impl<'a> Cache<'a> {
    fn new(g: &'a Mutex<Intern>) -> Self { Cache { g, l: HashMap::new(), last: (String::new(), u32::MAX) } }
    fn get(&mut self, s: &str) -> u32 {
        if self.last.1 != u32::MAX && self.last.0 == s { return self.last.1; }
        let i = match self.l.get(s) { Some(&i) => i, None => { let i = self.g.lock().unwrap().get(s); self.l.insert(s.to_string(), i); i } };
        self.last = (s.to_string(), i); i
    }
}

struct Drug { name: String, cmax: f64, cmax_u: f64, moa: String, generic: String, inert: bool, full: bool, moa_use: bool, src: String }
impl Drug {
    fn pk(&self) -> bool { self.cmax_u.is_finite() && self.cmax.is_finite() }
    fn thr_scale(&self, b: usize) -> (f64, f64) {
        let fu = (self.cmax_u / self.cmax).min(1.0);
        match b { 0 => (self.cmax_u, 1.0), 1 => (self.cmax, 1.0), _ => (self.cmax_u, 1.0 / (1.0 + 0.1 * (1.0 - fu) / fu)) }
    }
    fn reach(&self, b: usize, d: usize) -> bool { let (t, s) = self.thr_scale(b); DOSES[d] * s <= t * (1.0 + 1e-9) }
    /// 1 = straddling, 2 = always above, 0 = all doses reachable
    fn stratum(&self, b: usize) -> u8 { if !self.reach(b, 0) { 2 } else if !self.reach(b, 2) { 1 } else { 0 } }
    fn in_set(&self, s: &str) -> bool { s == "all" || self.src == s }
}

fn read(path: &str, cols: &[&str]) -> Result<Vec<RecordBatch>> {
    let bld = ParquetRecordBatchReaderBuilder::try_new(File::open(path)?)?;
    let sch = bld.parquet_schema();
    let idx: Vec<usize> = cols.iter().map(|c| sch.columns().iter().position(|x| x.name() == *c).ok_or(anyhow!("{path}: no column {c}"))).collect::<Result<_>>()?;
    let mask = ProjectionMask::leaves(sch, idx);
    Ok(bld.with_projection(mask).with_batch_size(1 << 16).build()?.collect::<std::result::Result<Vec<_>, _>>()?)
}
fn sc(b: &RecordBatch, n: &str) -> Result<StringArray> {
    let a = cast(b.column_by_name(n).ok_or(anyhow!("no {n}"))?, &DataType::Utf8)?;
    Ok(a.as_any().downcast_ref::<StringArray>().unwrap().clone())
}
fn fc(b: &RecordBatch, n: &str) -> Result<Float32Array> {
    let a = cast(b.column_by_name(n).ok_or(anyhow!("no {n}"))?, &DataType::Float32)?;
    Ok(a.as_any().downcast_ref::<Float32Array>().unwrap().clone())
}
fn ic(b: &RecordBatch, n: &str) -> Result<Int64Array> {
    let a = cast(b.column_by_name(n).ok_or(anyhow!("no {n}"))?, &DataType::Int64)?;
    Ok(a.as_any().downcast_ref::<Int64Array>().unwrap().clone())
}
fn list(dir: &str, pre: &str) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir).unwrap().filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with(pre) && n.ends_with(".parquet")).map(|n| format!("{dir}/{n}")).collect();
    v.sort(); v
}
fn inter(a: &[u32], b: &[u32]) -> Vec<u32> {
    let (mut i, mut j, mut o) = (0, 0, Vec::new());
    while i < a.len() && j < b.len() { if a[i] == b[j] { o.push(a[i]); i += 1; j += 1 } else if a[i] < b[j] { i += 1 } else { j += 1 } } o
}
fn union(a: &[u32], b: &[u32]) -> Vec<u32> { let mut o = [a, b].concat(); o.sort_unstable(); o.dedup(); o }
fn diff(a: &[u32], b: &[u32]) -> Vec<u32> {
    let (mut j, mut o) = (0, Vec::new());
    for &x in a { while j < b.len() && b[j] < x { j += 1 } if j >= b.len() || b[j] != x { o.push(x) } } o
}
fn jaccard(a: &[u32], b: &[u32]) -> f64 { let i = inter(a, b).len() as f64; i / (a.len() as f64 + b.len() as f64 - i) }
fn out(dir: &str, name: &str, header: &str) -> Result<BufWriter<File>> {
    let mut w = BufWriter::new(File::create(format!("{dir}/{name}"))?); writeln!(w, "{header}")?; Ok(w)
}
fn summ(v: &[f64]) -> Value { json!({"n": v.len(), "median": jnum(med(v)), "q1": jnum(quantile(v, 0.25)), "q3": jnum(quantile(v, 0.75))}) }

type Key = (u16, u8, u16); // drug, dose, line

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    let (p2, inp, od, rj) = (a[1].as_str(), a[2].as_str(), a[3].as_str(), a[4].as_str());
    std::fs::create_dir_all(od)?;
    let t0 = std::time::Instant::now();
    let mut res = Map::new();
    // ---------- inputs
    let drugs: Vec<Drug> = load_tahoe_compounds(inp)?.into_iter().map(|c| Drug { name: c.name, cmax: c.cmax, cmax_u: c.cmax_u, moa: c.moa, generic: c.generic, inert: c.inert, full: c.full, moa_use: c.moa_use, src: c.src }).collect();
    let didx: HashMap<&str, u16> = drugs.iter().enumerate().map(|(i, d)| (d.name.as_str(), i as u16)).collect();
    let (lines, plates, genes) = (Mutex::new(Intern::default()), Mutex::new(Intern::default()), Mutex::new(Intern::default()));

    // ---------- per-contrast counts
    struct CRow { drug: u16, dose: u8, line: u16, plate: u16, n: [f64; 5] }
    let cfiles = list(p2, "c");
    let crow: Vec<CRow> = cfiles.par_iter().map(|f| -> Result<Vec<CRow>> {
        let mut v = Vec::new();
        for b in read(f, &["drug", "concentration", "Cell_Name_Vevo", "plate", "n_de_0.25", "n_de_0.5", "n_de_0.75", "n_de_1.0", "n_de_1.5"])? {
            let (dr, co, li, pl) = (sc(&b, "drug")?, fc(&b, "concentration")?, sc(&b, "Cell_Name_Vevo")?, sc(&b, "plate")?);
            let ns: Vec<Int64Array> = CUTN.iter().map(|c| ic(&b, &format!("n_de_{c}"))).collect::<Result<_>>()?;
            for r in 0..b.num_rows() {
                let Some(&d) = didx.get(dr.value(r)) else { continue };
                let mut n = [0.0; 5]; for k in 0..5 { n[k] = ns[k].value(r) as f64 }
                v.push(CRow { drug: d, dose: dose_idx(co.value(r)), line: lines.lock().unwrap().get(li.value(r)) as u16,
                              plate: plates.lock().unwrap().get(pl.value(r)) as u16, n });
            }
        }
        Ok(v)
    }).collect::<Result<Vec<_>>>()?.into_iter().flatten().collect();
    let nl = lines.lock().unwrap().names.len();
    eprintln!("contrast rows {} lines {} [{:.0}s]", crow.len(), nl, t0.elapsed().as_secs_f64());
    // inert floor per plate x line and cut
    let mut fl: HashMap<(u16, u16), [Vec<f64>; 5]> = HashMap::new();
    for r in crow.iter().filter(|r| drugs[r.drug as usize].inert) { let e = fl.entry((r.plate, r.line)).or_default(); for k in 0..5 { e[k].push(r.n[k]) } }
    let floor: HashMap<(u16, u16), [f64; 5]> = fl.into_iter().map(|(k, v)| { let mut o = [0.0; 5]; for c in 0..5 { o[c] = quantile(&v[c], 0.95) } (k, o) }).collect();
    #[derive(Default, Clone)] struct Acc { n: [f64; 5], ex: [f64; 5], c: f64 }
    let (mut kf, mut ka): (HashMap<Key, Acc>, HashMap<Key, Acc>) = (HashMap::new(), HashMap::new());
    for r in &crow {
        let key = (r.drug, r.dose, r.line);
        let e = ka.entry(key).or_default(); e.c += 1.0; for k in 0..5 { e.n[k] += r.n[k] }
        if let Some(f) = floor.get(&(r.plate, r.line)) { let e = kf.entry(key).or_default(); e.c += 1.0; for k in 0..5 { e.n[k] += r.n[k]; e.ex[k] += (r.n[k] - f[k]).max(0.0) } }
    }
    for m in [&mut kf, &mut ka] { for e in m.values_mut() { for k in 0..5 { e.n[k] /= e.c; e.ex[k] /= e.c } } }
    let ln = lines.lock().unwrap().names.clone();
    res.insert("data".into(), json!({"summary_files": cfiles.len(), "contrast_rows_by_plate": crow.len(), "contrasts": ka.len(), "cell_lines": nl,
        "compounds": ka.keys().map(|k| k.0).collect::<HashSet<_>>().len(), "plates": plates.lock().unwrap().names.len(), "plate_line_floors": floor.len(),
        "inert_median_n_de_by_dose": (0..3u8).map(|d| jnum(med(&crow.iter().filter(|r| drugs[r.drug as usize].inert && r.dose == d).map(|r| r.n[3]).collect::<Vec<_>>()))).collect::<Vec<_>>(),
        "all_median_n_de_by_dose": (0..3u8).map(|d| jnum(med(&ka.iter().filter(|(k, _)| k.1 == d).map(|(_, e)| e.n[3]).collect::<Vec<_>>()))).collect::<Vec<_>>()}));

    // per-drug summary and the selection check
    let mut w = out(od, "drug_summary.tsv", "drug\tmoa\tpk_source\tmedian_n_de_5uM\tmedian_log10_de5\tmedian_slope")?;
    let mut ds: Vec<(usize, f64, f64, f64)> = vec![];
    for (i, d) in drugs.iter().enumerate() {
        let (mut n5, mut y5, mut sl) = (vec![], vec![], vec![]);
        for l in 0..nl as u16 {
            let y: Vec<Option<f64>> = (0..3).map(|k| ka.get(&(i as u16, k, l)).map(|e| (e.n[3] + 1.0).log10())).collect();
            if let Some(v) = y[2] { y5.push(v); n5.push(ka[&(i as u16, 2, l)].n[3]) }
            if y.iter().all(|v| v.is_some()) { sl.push((y[2].unwrap() - y[0].unwrap()) / 2.0) } // least-squares slope over three log10-spaced doses
        }
        if y5.is_empty() { continue }
        let r = (i, med(&n5), med(&y5), med(&sl)); writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}", d.name, d.moa, d.src, r.1, r.2, r.3)?; ds.push(r);
    }
    drop(w);
    {
        let is_m = |r: &&(usize, f64, f64, f64)| drugs[r.0].src == "liston" && drugs[r.0].pk();
        let m: Vec<&(usize, f64, f64, f64)> = ds.iter().filter(is_m).collect(); let o: Vec<&(usize, f64, f64, f64)> = ds.iter().filter(|r| !is_m(r)).collect();
        let cls: HashSet<&str> = m.iter().map(|r| drugs[r.0].moa.as_str()).filter(|s| *s != "unclear").collect();
        let oc: Vec<&(usize, f64, f64, f64)> = o.iter().cloned().filter(|r| cls.contains(drugs[r.0].moa.as_str())).collect();
        let col = |v: &Vec<&(usize, f64, f64, f64)>, k: usize| -> Vec<f64> { v.iter().map(|r| match k { 1 => r.1, 2 => r.2, _ => r.3 }).filter(|x| x.is_finite()).collect() };
        let test = |x: &Vec<&(usize, f64, f64, f64)>, k: usize| { let (a, b) = (col(&m, k), col(x, k)); json!({"ks_p": jnum(ks_2samp(&a, &b).1), "mann_whitney_p": jnum(mann_whitney(&a, &b).1)}) };
        res.insert("selection".into(), json!({"matched": summ(&col(&m, 1)), "other_all": summ(&col(&o, 1)), "other_same_classes": summ(&col(&oc, 1)), "n_classes": cls.len(),
            "de5_all": test(&o, 2), "de5_same_classes": test(&oc, 2), "slope_all": test(&o, 3), "slope_same_classes": test(&oc, 3)}));
    }
    let mut w = out(od, "contrasts.tsv", "drug\tpk_source\tdose\tline\tn_plates\tn_de_1.0\tn_plates_with_floor\texcess_1.0")?;
    for (k, e) in &ka { let d = &drugs[k.0 as usize]; if !d.pk() && !d.inert { continue }
        let (c, x) = kf.get(k).map(|f| (f.c, f.ex[3])).unwrap_or((0.0, f64::NAN));
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", d.name, d.src, DOSES[k.1 as usize], ln[k.2 as usize], e.c, e.n[3], c, x)?; }
    drop(w);
    let mut w = out(od, "strata.tsv", "drug\tgeneric\tmoa\tpk_source\tcmax_uM\tcmax_unbound_uM\tbasis\tscale\tstratum\treach_0.05\treach_0.5\treach_5")?;
    for d in drugs.iter().filter(|d| d.pk()) { for b in 0..3 {
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", d.name, d.generic, d.moa, d.src, d.cmax, d.cmax_u, BASIS[b], d.thr_scale(b).1, d.stratum(b), d.reach(b, 0) as u8, d.reach(b, 1) as u8, d.reach(b, 2) as u8)?;
    } }
    drop(w);

    // ---------- DE gene sets (level = highest |log2FC| cutoff passed), union over plates
    let sets: HashMap<Key, Vec<u32>> = list(p2, "g").par_iter().map(|f| -> Result<HashMap<Key, Vec<u32>>> {
        let (mut gc, mut lc) = (Cache::new(&genes), Cache::new(&lines)); let mut m: HashMap<Key, Vec<u32>> = HashMap::new();
        for b in read(f, &["drug", "concentration", "Cell_Name_Vevo", "gene_name", "log2FoldChange"])? {
            let (dr, co, li, ge, lf) = (sc(&b, "drug")?, fc(&b, "concentration")?, sc(&b, "Cell_Name_Vevo")?, sc(&b, "gene_name")?, fc(&b, "log2FoldChange")?);
            for r in 0..b.num_rows() {
                let Some(&d) = didx.get(dr.value(r)) else { continue };
                let al = lf.value(r).abs(); let lev = CUTS.iter().filter(|c| al >= **c).count() as u32; if lev == 0 { continue }
                m.entry((d, dose_idx(co.value(r)), lc.get(li.value(r)) as u16)).or_default().push(gc.get(ge.value(r)) << 3 | (lev - 1));
            }
        }
        Ok(m)
    }).try_reduce(HashMap::new, |mut a, b| { for (k, v) in b { a.entry(k).or_default().extend(v) } Ok(a) })?;
    let sets: HashMap<Key, Vec<u32>> = sets.into_par_iter().map(|(k, mut v)| {
        v.sort_unstable(); let mut o: Vec<u32> = Vec::with_capacity(v.len());
        for x in v { if let Some(l) = o.last_mut() { if *l >> 3 == x >> 3 { *l = x; continue } } o.push(x) }
        (k, o)
    }).collect();
    eprintln!("gene sets {} [{:.0}s]", sets.len(), t0.elapsed().as_secs_f64());
    let empty: Vec<u32> = vec![];
    let set_at = |k: &Key, cut: usize| -> Vec<u32> { sets.get(k).unwrap_or(&empty).iter().filter(|x| (**x & 7) as usize >= cut).map(|x| x >> 3).collect() };

    // ---------- sweep: cut x basis x PK set x stratum
    let has = |m: &HashMap<Key, Acc>, i: u16| (0..nl as u16).any(|l| m.contains_key(&(i, 2, l)));
    let pk: Vec<u16> = (0..drugs.len() as u16).filter(|&i| drugs[i as usize].pk() && has(&kf, i)).collect();
    let pkall: Vec<u16> = (0..drugs.len() as u16).filter(|&i| drugs[i as usize].pk() && has(&ka, i)).collect();
    let mut w = out(od, "sweep.tsv", "cut\tbasis\tpk_set\tstratum\tn_compounds\tshare_lowest_dose_above\tshare_5uM_above_10x\tde_share_above\tde_share_above_lo\tde_share_above_hi\texcess_share_above\texcess_share_above_lo\texcess_share_above_hi\tn_floor_compounds\trecovery_n_compounds\trecovery_median_ratio\trecovery_share_ratio_le_0.5\trecovery_median_gene_share")?;
    let mut wq = out(od, "recovery_per_compound.tsv", "cut\tbasis\tdrug\tpk_source\ttop_reachable_uM\tn_lines\tratio_excess\tgene_recovery")?;
    let mut sweep: Vec<Value> = vec![];
    let mut rng = Rng(20261116);
    for cut in 0..5 { for b in 0..3 {
        let mut per: Vec<(u16, f64, f64)> = vec![];
        for &i in &pk {
            let d = &drugs[i as usize]; if d.stratum(b) != 1 { continue }
            let top = if d.reach(b, 1) { 1u8 } else { 0u8 };
            let (mut rs, mut sh) = (vec![], vec![]);
            for l in 0..nl as u16 {
                let (Some(x), Some(y)) = (kf.get(&(i, top, l)), kf.get(&(i, 2, l))) else { continue };
                if y.ex[cut] <= 0.0 { continue }
                rs.push(x.ex[cut] / y.ex[cut]);
                let s5 = set_at(&(i, 2, l), cut); if !s5.is_empty() { sh.push(inter(&s5, &set_at(&(i, top, l), cut)).len() as f64 / s5.len() as f64) }
            }
            if rs.is_empty() { continue }
            let (n, r, s) = (rs.len(), med(&rs), med(&sh));
            writeln!(wq, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", CUTN[cut], BASIS[b], d.name, d.src, DOSES[top as usize], n, r, s)?;
            per.push((i, r, s));
        }
        // per-drug totals for the share of DE calls from unreachable doses
        let tot = |m: &HashMap<Key, Acc>, ex: bool| -> HashMap<u16, (f64, f64)> { let mut t: HashMap<u16, (f64, f64)> = HashMap::new();
            for (k, e) in m { let d = &drugs[k.0 as usize]; if !d.pk() { continue } let v = if ex { e.ex[cut] } else { e.n[cut] };
                let x = t.entry(k.0).or_default(); x.1 += v; if !d.reach(b, k.1 as usize) { x.0 += v } } t };
        let (ta, tf) = (tot(&ka, false), tot(&kf, true));
        let share = |t: &HashMap<u16, (f64, f64)>, sub: &[u16], rng: &mut Rng| -> (f64, f64, f64) {
            let v: Vec<(f64, f64)> = sub.iter().filter_map(|i| t.get(i).cloned()).collect(); if v.is_empty() { return (f64::NAN, f64::NAN, f64::NAN) }
            let f = |x: &[(f64, f64)]| x.iter().map(|p| p.0).sum::<f64>() / x.iter().map(|p| p.1).sum::<f64>();
            let bs: Vec<f64> = (0..2000).map(|_| { let s: Vec<(f64, f64)> = (0..v.len()).map(|_| v[(rng.next() % v.len() as u64) as usize]).collect(); f(&s) }).collect();
            (f(&v), quantile(&bs, 0.025), quantile(&bs, 0.975)) };
        for set in SETS { for (sn, st) in STRATA {
            let f = |i: &u16| drugs[*i as usize].in_set(set) && (st == 9 || drugs[*i as usize].stratum(b) == st);
            let sa: Vec<u16> = pkall.iter().cloned().filter(f).collect(); let sf: Vec<u16> = pk.iter().cloned().filter(f).collect();
            if sa.is_empty() { continue }
            let q1a = sa.iter().filter(|&&i| !drugs[i as usize].reach(b, 0)).count() as f64 / sa.len() as f64;
            let q1b = sa.iter().filter(|&&i| { let (t, s) = drugs[i as usize].thr_scale(b); 5.0 * s > 10.0 * t }).count() as f64 / sa.len() as f64;
            let (qr, qe) = (share(&ta, &sa, &mut rng), share(&tf, &sf, &mut rng));
            let (mut cu, mut ct, mut nu, mut nr) = (0.0, 0.0, 0.0, 0.0);
            for (k, e) in &ka { if !sa.contains(&k.0) { continue } ct += 1.0; if drugs[k.0 as usize].reach(b, k.1 as usize) { nr += e.n[cut] } else { cu += 1.0; nu += e.n[cut] } }
            let (cond_share, mean_un, mean_re) = (cu / ct, nu / cu, nr / (ct - cu));
            let r: Vec<f64> = per.iter().filter(|p| sf.contains(&p.0)).map(|p| p.1).collect();
            let s: Vec<f64> = per.iter().filter(|p| sf.contains(&p.0)).map(|p| p.2).collect();
            let le = r.iter().filter(|x| **x <= 0.5).count() as f64 / r.len().max(1) as f64;
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", CUTN[cut], BASIS[b], set, sn, sa.len(), q1a, q1b, qr.0, qr.1, qr.2, qe.0, qe.1, qe.2, sf.len(), r.len(), med(&r), le, med(&s))?;
            sweep.push(json!({"cut": CUTN[cut], "basis": BASIS[b], "pk_set": set, "stratum": sn, "n_compounds": sa.len(), "share_lowest_dose_above": q1a, "share_5uM_above_10x": q1b,
                "share_contrasts_unreachable": jnum(cond_share), "mean_n_de_unreachable": jnum(mean_un), "mean_n_de_reachable": jnum(mean_re),
                "de_share_above": jnum(qr.0), "de_share_above_ci": [jnum(qr.1), jnum(qr.2)], "excess_share_above": jnum(qe.0), "excess_share_above_ci": [jnum(qe.1), jnum(qe.2)], "n_floor_compounds": sf.len(),
                "recovery_n_compounds": r.len(), "recovery_median_ratio": jnum(med(&r)), "recovery_share_ratio_le_0.5": le, "recovery_median_gene_share": jnum(med(&s))}));
        } }
    } }
    drop(w); drop(wq);
    res.insert("sweep".into(), Value::Array(sweep));
    {   // response by exposure relative to unbound Cmax (Liston set, |log2FC| >= 1, contrasts on plates with inert compounds)
        const BN: [&str; 5] = ["<=0.1x", "0.1-1x", "1-10x", "10-100x", ">100x"];
        let binof = |d: &Drug, k: usize| { let r = (DOSES[k] / d.cmax_u).log10(); if r <= -1.0 { 0 } else if r <= 0.0 { 1 } else if r <= 1.0 { 2 } else if r <= 2.0 { 3 } else { 4 } };
        let mut bins: Vec<Vec<(f64, f64)>> = vec![vec![]; 5]; let mut rs: [Vec<(f64, f64)>; 2] = [vec![], vec![]];
        for (k, e) in &kf { let d = &drugs[k.0 as usize]; if !(d.pk() && d.src == "liston") { continue }
            bins[binof(d, k.1 as usize)].push((e.n[3], e.ex[3])); rs[d.reach(0, k.1 as usize) as usize].push((e.n[3], e.ex[3])) }
        let inert_above = { let v: Vec<f64> = kf.iter().filter(|(k, _)| drugs[k.0 as usize].inert).map(|(_, e)| (e.ex[3] > 0.0) as u8 as f64).collect(); v.iter().sum::<f64>() / v.len() as f64 };
        let st = |v: &Vec<(f64, f64)>| json!({"n": v.len(), "share_above_floor": v.iter().filter(|x| x.1 > 0.0).count() as f64 / v.len().max(1) as f64, "median_n_de": jnum(med(&v.iter().map(|x| x.0).collect::<Vec<_>>())),
            "mean_n_de": v.iter().map(|x| x.0).sum::<f64>() / v.len().max(1) as f64, "mean_excess": v.iter().map(|x| x.1).sum::<f64>() / v.len().max(1) as f64});
        // within straddling compounds: share of lines above the floor at the highest reachable dose and at 5 uM
        let (mut lo, mut hi) = (vec![], vec![]);
        for &i in &pk { let d = &drugs[i as usize]; if d.src != "liston" || d.stratum(0) != 1 { continue } let top = if d.reach(0, 1) { 1u8 } else { 0u8 };
            let f = |k: u8| { let v: Vec<f64> = (0..nl as u16).filter_map(|l| kf.get(&(i, k, l)).map(|e| (e.ex[3] > 0.0) as u8 as f64)).collect(); v.iter().sum::<f64>() / v.len().max(1) as f64 };
            lo.push(f(top)); hi.push(f(2)) }
        let wt = wilcoxon(&hi, &lo);
        res.insert("exposure_response".into(), json!({"by_exposure": BN.iter().zip(&bins).map(|(n, v)| (n.to_string(), st(v))).collect::<Map<_, _>>(), "reachable": st(&rs[1]), "unreachable": st(&rs[0]),
            "inert_share_above_floor": inert_above, "inert_median_n_de": jnum(med(&kf.iter().filter(|(k, _)| drugs[k.0 as usize].inert).map(|(_, e)| e.n[3]).collect::<Vec<_>>())),
            "inert_median_n_de_by_dose": (0..3u8).map(|dz| jnum(med(&kf.iter().filter(|(k, _)| drugs[k.0 as usize].inert && k.1 == dz).map(|(_, e)| e.n[3]).collect::<Vec<_>>()))).collect::<Vec<_>>(), "straddling": {"n_compounds": lo.len(), "median_share_above_floor_highest_reachable": jnum(med(&lo)), "median_share_above_floor_5uM": jnum(med(&hi)), "wilcoxon_p_two_sided": jnum(wt.0)}}));
    }

    // ---------- cross-mechanism overlap of high-dose-only genes (|log2FC| >= 1), unbound and FBS basis, by stratum
    let stress: Vec<(String, HashSet<u32>)> = {
        let mut m: Vec<(String, HashSet<u32>)> = vec![]; let mut g = genes.lock().unwrap();
        for l in std::fs::read_to_string(format!("{inp}/gene_sets/stress.tsv"))?.lines() { let f: Vec<&str> = l.split('\t').collect(); if f.len() < 2 { continue }
            let id = g.get(f[1]); match m.iter_mut().find(|x| x.0 == f[0]) { Some(x) => { x.1.insert(id); } None => m.push((f[0].to_string(), HashSet::from([id]))) } }
        let all: HashSet<u32> = m.iter().flat_map(|x| x.1.iter().cloned()).collect(); m.push(("ANY_STRESS".into(), all)); m
    };
    let mut seen = HashSet::new();
    let a3d: Vec<u16> = pkall.iter().cloned().filter(|&i| { let d = &drugs[i as usize]; d.src == "liston" && d.moa != "unclear" && seen.insert(d.generic.clone()) }).collect();
    let mut wp = out(od, "overlap_pairs.tsv", "basis\tstratum\tline\tdrug_a\tdrug_b\tmoa_a\tmoa_b\tj_high_only\tj_low\tn_high_a\tn_high_b\tn_low_a\tn_low_b")?;
    let mut wu = out(od, "overlap_sets.tsv", "basis\tstratum\tdrug\tline\tn_low\tn_high_only\tn_shared\tn_5uM")?;
    let mut a3 = Map::new();
    for b in [0usize, 2] { for st in [1u8, 2] {
        struct U { d: u16, l: u16, lo: Vec<u32>, hi: Vec<u32>, sh: Vec<u32> }
        let mut us: Vec<U> = vec![];
        for &i in &a3d {
            let d = &drugs[i as usize]; if d.stratum(b) != st { continue }
            let low: Vec<u8> = if st == 1 { (0..2u8).filter(|&k| d.reach(b, k as usize)).collect() } else { vec![0, 1] };
            let top = if st == 1 { *low.last().unwrap() } else { 0 };
            for l in 0..nl as u16 {
                if !ka.contains_key(&(i, 2, l)) { continue }
                let s5 = set_at(&(i, 2, l), 3); let mut un: Vec<u32> = vec![]; for &k in &low { un = union(&un, &set_at(&(i, k, l), 3)) }
                let u = U { d: i, l, lo: set_at(&(i, top, l), 3), hi: diff(&s5, &un), sh: inter(&s5, &un) };
                writeln!(wu, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", BASIS[b], st, d.name, ln[l as usize], u.lo.len(), u.hi.len(), u.sh.len(), s5.len())?;
                us.push(u);
            }
        }
        let (mut jh, mut jl) = (vec![], vec![]); let mut dset = HashSet::new(); let mut lset = HashSet::new();
        for (x, u) in us.iter().enumerate() { for v in &us[x + 1..] {
            let (da, db) = (&drugs[u.d as usize], &drugs[v.d as usize]);
            if u.l != v.l || da.moa == db.moa || u.hi.len().min(v.hi.len()).min(u.lo.len()).min(v.lo.len()) < 10 { continue }
            let (h, l) = (jaccard(&u.hi, &v.hi), jaccard(&u.lo, &v.lo)); jh.push(h); jl.push(l); dset.insert(u.d); dset.insert(v.d); lset.insert(u.l);
            writeln!(wp, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", BASIS[b], st, ln[u.l as usize], da.name, db.name, da.moa, db.moa, h, l, u.hi.len(), v.hi.len(), u.lo.len(), v.lo.len())?;
        } }
        let wt = wilcoxon(&jh, &jl); let mut sj = Map::new();
        for (gs, set) in &stress {
            let cnt = |pick: usize| -> (u64, u64) { let (mut a_, mut t_) = (0u64, 0u64);
                for u in &us { let s = match pick { 0 => &u.hi, 1 => &u.sh, _ => &u.lo }; t_ += s.len() as u64; a_ += s.iter().filter(|g| set.contains(g)).count() as u64 } (a_, t_) };
            let (h, s, l) = (cnt(0), cnt(1), cnt(2)); let fz = fisher(h.0, h.1 - h.0, s.0, s.1 - s.0);
            sj.insert(gs.clone(), json!({"high_only": [h.0, h.1], "shared": [s.0, s.1], "low": [l.0, l.1], "share_high_only": h.0 as f64 / h.1.max(1) as f64,
                "share_shared": s.0 as f64 / s.1.max(1) as f64, "share_low": l.0 as f64 / l.1.max(1) as f64, "odds_ratio_high_only_vs_shared": jnum(fz.0), "fisher_p": jnum(fz.1)}));
        }
        a3.insert(format!("{}_stratum{}", BASIS[b], st), json!({"n_pairs": jh.len(), "n_drugs": dset.len(), "n_lines": lset.len(), "median_j_high_only": jnum(med(&jh)), "median_j_low": jnum(med(&jl)),
            "ratio_of_medians": jnum(med(&jh) / med(&jl)), "wilcoxon_p_high_gt_low": jnum(wt.1), "wilcoxon_p_low_gt_high": jnum(wt.2), "stress": Value::Object(sj)}));
    } }
    drop(wp); drop(wu);
    res.insert("overlap".into(), Value::Object(a3));
    eprintln!("exposure tables done [{:.0}s]", t0.elapsed().as_secs_f64());

    // ---------- full vectors: signatures for the patient comparison and mechanism discriminability
    let fd: Vec<u16> = (0..drugs.len() as u16).filter(|&i| drugs[i as usize].full).collect();
    let fpos: HashMap<u16, usize> = fd.iter().enumerate().map(|(k, &i)| (i, k)).collect();
    struct FR { f: u8, dose: u8, line: u16, gene: u32, lfc: f32, padj: f32 }
    let frows: Vec<Vec<FR>> = list(p2, "f").par_iter().map(|f| -> Result<Vec<FR>> {
        let (mut gc, mut lc) = (Cache::new(&genes), Cache::new(&lines)); let mut v = Vec::new();
        for b in read(f, &["drug", "concentration", "Cell_Name_Vevo", "gene_name", "log2FoldChange", "padj"])? {
            let (dr, co, li, ge, lf, pa) = (sc(&b, "drug")?, fc(&b, "concentration")?, sc(&b, "Cell_Name_Vevo")?, sc(&b, "gene_name")?, fc(&b, "log2FoldChange")?, fc(&b, "padj")?);
            for r in 0..b.num_rows() {
                let Some(k) = didx.get(dr.value(r)).and_then(|d| fpos.get(d)) else { continue };
                let (l, p) = (lf.value(r), if pa.is_null(r) { f32::NAN } else { pa.value(r) }); if !l.is_finite() { continue }
                v.push(FR { f: *k as u8, dose: dose_idx(co.value(r)), line: lc.get(li.value(r)) as u16, gene: gc.get(ge.value(r)), lfc: l, padj: p });
            }
        }
        Ok(v)
    }).collect::<Result<_>>()?;
    let (ng, nlf) = (genes.lock().unwrap().names.len(), lines.lock().unwrap().names.len());
    let nc = fd.len() * 3 * nlf;
    let ix = |f: usize, d: usize, l: usize| (f * 3 + d) * nlf + l;
    let (mut sl, mut sp, mut cn, mut cp) = (vec![0f32; nc * ng], vec![0f32; nc * ng], vec![0u8; nc * ng], vec![0u8; nc * ng]);
    for v in &frows { for r in v { let j = ix(r.f as usize, r.dose as usize, r.line as usize) * ng + r.gene as usize;
        sl[j] += r.lfc; cn[j] += 1; if r.padj.is_finite() { sp[j] += r.padj; cp[j] += 1 } } }
    drop(frows);
    sl.par_iter_mut().zip(cn.par_iter()).for_each(|(s, &c)| if c > 0 { *s /= c as f32 } else { *s = f32::NAN });
    sp.par_iter_mut().zip(cp.par_iter()).for_each(|(s, &c)| if c > 0 { *s /= c as f32 } else { *s = f32::NAN });
    eprintln!("full vectors loaded: {} drugs, {} genes [{:.0}s]", fd.len(), ng, t0.elapsed().as_secs_f64());
    let gn = genes.lock().unwrap().names.clone(); let ln = lines.lock().unwrap().names.clone();
    let rows: Vec<String> = (0..fd.len() * 3).into_par_iter().map(|fdz| {
        let (f, d) = (fdz / 3, fdz % 3); let mut s = String::new();
        for g in 0..ng {
            let (mut up, mut dn) = (0, 0); let mut v: Vec<f64> = vec![];
            for l in 0..nlf { let j = ix(f, d, l) * ng + g; if cn[j] == 0 { continue } v.push(sl[j] as f64);
                if sp[j] < 0.05 { if sl[j] >= 1.0 { up += 1 } else if sl[j] <= -1.0 { dn += 1 } } }
            if v.is_empty() { continue } let n = v.len();
            s.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\n", drugs[fd[f] as usize].name, DOSES[d], gn[g], n, up, dn, median(&mut v)));
        }
        s
    }).collect();
    let mut w = out(od, "tahoe_signatures.tsv", "drug\tdose\tgene\tn_lines\tn_up\tn_dn\tmedian_lfc")?; for s in rows { w.write_all(s.as_bytes())? } drop(w);

    // mechanism discriminability
    let mu: Vec<usize> = (0..fd.len()).filter(|&k| drugs[fd[k] as usize].moa_use).collect();
    let present: Vec<bool> = (0..nc).into_par_iter().map(|c| (0..ng).filter(|&g| cn[c * ng + g] > 0).count() >= 1000).collect();
    let cs: Vec<usize> = mu.iter().flat_map(|&f| (0..3).flat_map(move |d| (0..nlf).map(move |l| (f, d, l)))).map(|(f, d, l)| ix(f, d, l)).filter(|&c| present[c]).collect();
    let mut gv: Vec<(f64, usize)> = (0..ng).into_par_iter().filter_map(|g| {
        let v: Vec<f64> = cs.iter().filter(|&&c| cn[c * ng + g] > 0).map(|&c| sl[c * ng + g] as f64).collect();
        if (v.len() as f64) < 0.8 * cs.len() as f64 { return None }
        let m = v.iter().sum::<f64>() / v.len() as f64; Some((v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64, g))
    }).collect();
    gv.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap()); let top: Vec<usize> = gv.iter().take(1000).map(|x| x.1).collect();
    let mut w = out(od, "moa.tsv", "line\tdose\tn_drugs\tn_genes\tn_within\tn_between\tauroc\tmean_within\tmean_between")?;
    let mut wp = out(od, "moa_pairs.tsv", "line\tdose\tdrug_a\tdrug_b\tsame_moa\tspearman\tbin_a\tbin_b")?;
    const BIN: [&str; 4] = ["<=1x", "1-10x", "10-100x", ">100x"];
    let bin = |d: &Drug, dose: usize| -> usize { let r = (DOSES[dose] / d.cmax_u).log10(); if r <= 0.0 { 0 } else if r <= 1.0 { 1 } else if r <= 2.0 { 2 } else { 3 } };
    // per dose and per exposure bin: per-line within and between similarities
    let mut by_dose: Vec<Vec<(f64, f64, f64)>> = vec![vec![(f64::NAN, f64::NAN, f64::NAN); nlf]; 3];
    let mut by_bin: Vec<Vec<(Vec<f64>, Vec<f64>)>> = vec![vec![(vec![], vec![]); nlf]; 5]; // bins 0-3 and 4 = above 10x
    let mut ngenes = vec![];
    for l in 0..nlf { for d in 0..3 {
        let dsx: Vec<usize> = mu.iter().cloned().filter(|&f| present[ix(f, d, l)]).collect(); if dsx.len() < 4 { continue }
        let gs: Vec<usize> = top.iter().cloned().filter(|&g| dsx.iter().all(|&f| cn[ix(f, d, l) * ng + g] > 0)).collect(); if gs.len() < 100 { continue }
        let rk: Vec<Vec<f64>> = dsx.iter().map(|&f| ranks(&gs.iter().map(|&g| sl[ix(f, d, l) * ng + g] as f64).collect::<Vec<f64>>()).0).collect();
        let (mut wi, mut be) = (vec![], vec![]);
        for x in 0..dsx.len() { for y in x + 1..dsx.len() {
            let (da, db) = (&drugs[fd[dsx[x]] as usize], &drugs[fd[dsx[y]] as usize]); let r = pearson(&rk[x], &rk[y]); let same = da.moa == db.moa;
            let (ba, bb) = (bin(da, d), bin(db, d));
            writeln!(wp, "{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{}", ln[l], DOSES[d], da.name, db.name, same as u8, r, BIN[ba], BIN[bb])?;
            if same { wi.push(r) } else { be.push(r) }
            if ba == bb { let e = &mut by_bin[ba][l]; if same { e.0.push(r) } else { e.1.push(r) } }
            if ba >= 2 && bb >= 2 { let e = &mut by_bin[4][l]; if same { e.0.push(r) } else { e.1.push(r) } }
        } }
        let m = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len().max(1) as f64;
        by_dose[d][l] = (auroc(&wi, &be), m(&wi), m(&be)); ngenes.push(gs.len() as f64);
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", ln[l], DOSES[d], dsx.len(), gs.len(), wi.len(), be.len(), by_dose[d][l].0, by_dose[d][l].1, by_dose[d][l].2)?;
    } }
    drop(w); drop(wp);
    {
        let ok: Vec<usize> = (0..nlf).filter(|&l| (0..3).all(|d| by_dose[d][l].0.is_finite())).collect();
        let col = |d: usize, k: usize| -> Vec<f64> { ok.iter().map(|&l| match k { 0 => by_dose[d][l].0, 1 => by_dose[d][l].1, _ => by_dose[d][l].2 }).collect() };
        let half = vec![0.5; ok.len()];
        let mut m = Map::new();
        m.insert("n_lines".into(), json!(ok.len())); m.insert("n_drugs".into(), json!(mu.len())); m.insert("median_genes".into(), jnum(med(&ngenes)));
        let mut cl: HashMap<&str, usize> = HashMap::new(); for &k in &mu { *cl.entry(drugs[fd[k] as usize].moa.as_str()).or_default() += 1 }
        m.insert("classes".into(), json!(cl));
        for d in 0..3 { m.insert(format!("dose_{}", DOSES[d]), json!({"median_auroc": jnum(med(&col(d, 0))), "auroc_q1": jnum(quantile(&col(d, 0), 0.25)), "auroc_q3": jnum(quantile(&col(d, 0), 0.75)),
            "auroc_vs_half_p": jnum(wilcoxon(&col(d, 0), &half).0), "median_within": jnum(med(&col(d, 1))), "median_between": jnum(med(&col(d, 2)))})); }
        let (wa, wb, ww) = (wilcoxon(&col(2, 0), &col(0, 0)), wilcoxon(&col(2, 2), &col(0, 2)), wilcoxon(&col(2, 1), &col(0, 1)));
        m.insert("auroc_5_vs_0.05".into(), json!({"p_two_sided": jnum(wa.0), "p_5_lower": jnum(wa.2), "p_5_higher": jnum(wa.1)}));
        m.insert("between_5_vs_0.05".into(), json!({"p_two_sided": jnum(wb.0), "p_5_higher": jnum(wb.1)}));
        m.insert("within_5_vs_0.05".into(), json!({"p_two_sided": jnum(ww.0), "p_5_higher": jnum(ww.1)}));
        let mut bj = Map::new(); let mut per_line: Vec<Vec<f64>> = vec![];
        for (bi, name) in ["<=1x", "1-10x", "10-100x", ">100x", ">10x"].iter().enumerate() {
            let (wi, be): (Vec<f64>, Vec<f64>) = (by_bin[bi].iter().flat_map(|e| e.0.clone()).collect(), by_bin[bi].iter().flat_map(|e| e.1.clone()).collect());
            let pl: Vec<f64> = (0..nlf).map(|l| { let e = &by_bin[bi][l]; if e.0.len() >= 2 && e.1.len() >= 2 { auroc(&e.0, &e.1) } else { f64::NAN } }).collect();
            let f: Vec<f64> = pl.iter().cloned().filter(|x| x.is_finite()).collect();
            bj.insert(name.to_string(), json!({"n_within": wi.len(), "n_between": be.len(), "pooled_auroc": jnum(auroc(&wi, &be)), "mean_within": jnum(wi.iter().sum::<f64>() / wi.len().max(1) as f64),
                "mean_between": jnum(be.iter().sum::<f64>() / be.len().max(1) as f64), "n_lines": f.len(), "median_line_auroc": jnum(med(&f)), "line_auroc_q1": jnum(quantile(&f, 0.25)),
                "line_auroc_q3": jnum(quantile(&f, 0.75)), "line_auroc_vs_half_p": jnum(wilcoxon(&f, &vec![0.5; f.len()]).0)}));
            per_line.push(pl);
        }
        let both: Vec<usize> = (0..nlf).filter(|&l| per_line[0][l].is_finite() && per_line[4][l].is_finite()).collect();
        let pw = wilcoxon(&both.iter().map(|&l| per_line[4][l]).collect::<Vec<_>>(), &both.iter().map(|&l| per_line[0][l]).collect::<Vec<_>>());
        bj.insert(">10x_vs_<=1x".into(), json!({"n_lines": both.len(), "p_two_sided": jnum(pw.0), "p_above_10x_higher": jnum(pw.1)}));
        m.insert("by_exposure".into(), Value::Object(bj));
        res.insert("mechanism".into(), Value::Object(m));
    }
    std::fs::write(rj, serde_json::to_string_pretty(&Value::Object(res))?)?;
    eprintln!("done [{:.0}s]", t0.elapsed().as_secs_f64());
    Ok(())
}
