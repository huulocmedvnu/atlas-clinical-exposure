//! Step 8. Agreement of Tahoe-100M signatures with the change in malignant cells of single-cell patient cohorts.
//! Usage:
//!   sc_cohorts signature <step 1 output dir> <compound> <out tsv>      signature table of one compound from its full profiles
//!   sc_cohorts feline <GSE158724 10x counts .txt.gz> <out tsv>         counts summed over cancer cells per sample
//!   sc_cohorts tenx <GSE292469 counts .h5> <cell metadata .tsv.gz> <out tsv>   counts summed over tumour cells per sample
//!   sc_cohorts compare <data dir> <ribociclib signature tsv> <feline sums> <osimertinib sums or -> <results json>
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use arrow::array::*;
use arrow::compute::cast;
use arrow::datatypes::DataType;
use atlas_exposure::*;
use flate2::read::MultiGzDecoder;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rayon::prelude::*;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

const DOSES: [f64; 3] = [0.05, 0.5, 5.0];
const MIN_CELLS: f64 = 30.0;
const MIN_CPM: f64 = 5.0;

fn gz_lines(path: &str) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path)?;
    Ok(if path.ends_with(".gz") { Box::new(BufReader::with_capacity(1 << 22, MultiGzDecoder::new(f))) } else { Box::new(BufReader::with_capacity(1 << 22, f)) })
}
fn wilson(k: f64, n: f64) -> (f64, f64) { if n == 0.0 { return (f64::NAN, f64::NAN) } let (z, p) = (1.959964, k / n); let d = 1.0 + z * z / n;
    let (c, h) = ((p + z * z / (2.0 * n)) / d, z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / d); (c - h, c + h) }
fn mean(v: &[f64]) -> f64 { v.iter().sum::<f64>() / v.len() as f64 }
fn var(v: &[f64]) -> f64 { let m = mean(v); v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0) }

/// Signature table of one compound: per dose and gene, cell lines measured, cell lines with up- and down-regulation, median log2FC.
fn signature(dir: &str, drug: &str, out: &str) -> Result<()> {
    let mut files: Vec<String> = std::fs::read_dir(dir)?.filter_map(|e| e.ok()).map(|e| e.file_name().to_string_lossy().to_string()).filter(|n| n.starts_with('f') && n.ends_with(".parquet")).map(|n| format!("{dir}/{n}")).collect();
    files.sort();
    let parts: Vec<Vec<(u8, String, String, f32, f32)>> = files.par_iter().map(|f| -> Result<Vec<(u8, String, String, f32, f32)>> {
        let mut v = vec![];
        for b in ParquetRecordBatchReaderBuilder::try_new(std::fs::File::open(f)?)?.with_batch_size(1 << 16).build()? { let b = b?;
            let st = |n: &str| -> Result<StringArray> { Ok(cast(b.column_by_name(n).ok_or(anyhow!("no {n}"))?, &DataType::Utf8)?.as_any().downcast_ref::<StringArray>().unwrap().clone()) };
            let fl = |n: &str| -> Result<Float32Array> { Ok(cast(b.column_by_name(n).ok_or(anyhow!("no {n}"))?, &DataType::Float32)?.as_any().downcast_ref::<Float32Array>().unwrap().clone()) };
            let (dr, co, li, ge, lf, pa) = (st("drug")?, fl("concentration")?, st("Cell_Name_Vevo")?, st("gene_name")?, fl("log2FoldChange")?, fl("padj")?);
            for r in 0..b.num_rows() { if dr.value(r) != drug || lf.is_null(r) || !lf.value(r).is_finite() { continue }
                let c = co.value(r) as f64; let k = if c < 0.2 { 0 } else if c < 2.0 { 1 } else { 2 };
                v.push((k, li.value(r).to_string(), ge.value(r).to_string(), lf.value(r), if pa.is_null(r) { f32::NAN } else { pa.value(r) })); } }
        Ok(v) }).collect::<Result<_>>()?;
    // mean over plates per dose, cell line and gene
    let mut acc: HashMap<(u8, String, String), (f64, f64, f64, f64)> = HashMap::new();
    for p in parts { for (k, l, g, lfc, padj) in p { let e = acc.entry((k, g, l)).or_insert((0.0, 0.0, 0.0, 0.0)); e.0 += lfc as f64; e.1 += 1.0; if padj.is_finite() { e.2 += padj as f64; e.3 += 1.0 } } }
    let mut per: BTreeMap<(u8, String), (u32, u32, Vec<f64>)> = BTreeMap::new();
    for ((k, g, _), (sl, n, sp, np)) in acc { let (l, p) = (sl / n, if np > 0.0 { sp / np } else { f64::NAN }); let e = per.entry((k, g)).or_insert((0, 0, vec![]));
        e.2.push(l); if p < 0.05 { if l >= 1.0 { e.0 += 1 } else if l <= -1.0 { e.1 += 1 } } }
    let mut w = BufWriter::new(std::fs::File::create(out)?); writeln!(w, "drug\tdose\tgene\tn_lines\tn_up\tn_dn\tmedian_lfc")?;
    for ((k, g), (up, dn, mut v)) in per { writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}", drug, DOSES[k as usize], g, v.len(), up, dn, median(&mut v))?; }
    Ok(())
}

fn write_sums(out: &str, samples: &[String], ncell: &[f64], genes: &[String], sums: &[Vec<f64>]) -> Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(out)?);
    writeln!(w, "gene\t{}", samples.join("\t"))?; writeln!(w, "n_cells\t{}", ncell.iter().map(|x| format!("{x}")).collect::<Vec<_>>().join("\t"))?;
    for (g, s) in genes.iter().zip(sums) { writeln!(w, "{}\t{}", g, s.iter().map(|x| format!("{x}")).collect::<Vec<_>>().join("\t"))?; }
    Ok(())
}
/// GSE158724: tab-separated genes x cells, cell names `<patient>_<time point>_<barcode>`, cancer cells only.
fn feline(path: &str, out: &str) -> Result<()> {
    let mut rd = gz_lines(path)?; let mut line = String::new(); rd.read_line(&mut line)?;
    let mut samples: Vec<String> = vec![]; let mut pos: HashMap<String, usize> = HashMap::new(); let mut col = vec![];
    for c in line.trim_end().split('\t').skip(1) { let mut it = c.splitn(3, '_'); let s = format!("{}_{}", it.next().unwrap_or(""), it.next().unwrap_or(""));
        let i = *pos.entry(s.clone()).or_insert_with(|| { samples.push(s.clone()); samples.len() - 1 }); col.push(i); }
    let mut ncell = vec![0.0; samples.len()]; for &i in &col { ncell[i] += 1.0 }
    let (mut genes, mut sums) = (vec![], vec![]);
    loop { line.clear(); if rd.read_line(&mut line)? == 0 { break }
        let mut it = line.trim_end().split('\t'); let g = it.next().unwrap_or("").to_string(); let mut s = vec![0.0; samples.len()];
        for (j, v) in it.enumerate() { if v != "0" { s[col[j]] += v.parse::<f64>().unwrap_or(0.0) } }
        genes.push(g); sums.push(s); }
    write_sums(out, &samples, &ncell, &genes, &sums)
}
/// GSE292469: 10x HDF5 matrix (cells in columns) and the per-nucleus table of the authors. Tumour cells of treatment-naive (TN) and residual disease (MRD) samples.
fn tenx(h5: &str, meta: &str, out: &str) -> Result<()> {
    let mut txt = String::new(); MultiGzDecoder::new(std::fs::File::open(meta)?).read_to_string(&mut txt)?;
    let mut it = txt.lines(); let h: Vec<&str> = it.next().ok_or(anyhow!("empty table"))?.split('\t').collect();
    let c = |n: &str| h.iter().position(|x| *x == n).ok_or(anyhow!("column {n} missing")); let (ci, cs, ct, cf) = (c("cell_id")?, c("sample_id")?, c("time_point")?, c("is_tumor_cell")?);
    let mut cell: HashMap<String, String> = HashMap::new();
    for l in it { let f: Vec<&str> = l.split('\t').collect(); if f[cf].eq_ignore_ascii_case("true") && (f[ct] == "TN" || f[ct] == "MRD") { cell.insert(f[ci].to_string(), f[cs].to_string()); } }
    let f = hdf5::File::open(h5)?; let grp = if f.link_exists("matrix") { f.group("matrix")? } else { f.group("/")? };
    let strs = |n: &str| -> Result<Vec<String>> { let d = grp.dataset(n)?;
        if let Ok(v) = d.read_raw::<hdf5::types::VarLenUnicode>() { return Ok(v.iter().map(|s| s.to_string()).collect()) }
        if let Ok(v) = d.read_raw::<hdf5::types::VarLenAscii>() { return Ok(v.iter().map(|s| s.to_string()).collect()) }
        macro_rules! fixed { ($($n:literal),*) => { $( if let Ok(v) = d.read_raw::<hdf5::types::FixedAscii<$n>>() { return Ok(v.iter().map(|s| s.to_string()).collect()) } if let Ok(v) = d.read_raw::<hdf5::types::FixedUnicode<$n>>() { return Ok(v.iter().map(|s| s.to_string()).collect()) } )* } }
        fixed!(16, 24, 32, 48, 64, 96, 128, 256);
        Err(anyhow!("{n}: unsupported string type")) };
    let barcodes = strs("barcodes")?; let genes = if grp.link_exists("features") { strs("features/name")? } else { strs("gene_names")? };
    let indptr: Vec<u64> = grp.dataset("indptr")?.read_raw()?; let (di, dd) = (grp.dataset("indices")?, grp.dataset("data")?);
    let mut samples: Vec<String> = vec![]; let mut pos: HashMap<String, usize> = HashMap::new();
    let col: Vec<Option<usize>> = barcodes.iter().map(|b| cell.get(b).map(|s| *pos.entry(s.clone()).or_insert_with(|| { samples.push(s.clone()); samples.len() - 1 }))).collect();
    eprintln!("{} nuclei in the matrix, {} tumour nuclei of TN and MRD samples matched, {} samples", barcodes.len(), col.iter().flatten().count(), samples.len());
    let mut ncell = vec![0.0; samples.len()]; let mut sums = vec![vec![0.0; samples.len()]; genes.len()];
    let step = 20000usize; let mut a = 0;
    while a < barcodes.len() { let b = (a + step).min(barcodes.len()); let (lo, hi) = (indptr[a] as usize, indptr[b] as usize);
        if hi > lo { let idx: Vec<i64> = di.read_slice_1d::<i64, _>(lo..hi)?.to_vec(); let dat: Vec<f64> = dd.read_slice_1d::<f64, _>(lo..hi)?.to_vec();
            for j in a..b { if let Some(s) = col[j] { ncell[s] += 1.0; for k in indptr[j] as usize - lo..indptr[j + 1] as usize - lo { sums[idx[k] as usize][s] += dat[k] } } } }
        a = b; }
    // sum rows with the same gene name
    let mut by: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for (g, s) in genes.iter().zip(sums) { let e = by.entry(g.clone()).or_insert_with(|| vec![0.0; samples.len()]); for (x, y) in e.iter_mut().zip(s) { *x += y } }
    write_sums(out, &samples, &ncell, &by.keys().cloned().collect::<Vec<_>>(), &by.values().cloned().collect::<Vec<_>>())
}

struct Pb { samples: Vec<String>, genes: Vec<String>, x: Vec<Vec<f64>> } // log2(CPM + 1), genes x samples, filtered
fn load_sums(path: &str) -> Result<Pb> {
    let rd = gz_lines(path)?; let mut it = rd.lines(); let h = it.next().ok_or(anyhow!("empty"))??; let all: Vec<String> = h.split('\t').skip(1).map(String::from).collect();
    let nc: Vec<f64> = it.next().ok_or(anyhow!("no n_cells"))??.split('\t').skip(1).map(|x| x.parse().unwrap_or(0.0)).collect();
    let keep: Vec<usize> = (0..all.len()).filter(|&i| nc[i] >= MIN_CELLS).collect();
    let (mut genes, mut raw) = (vec![], vec![]);
    for l in it { let l = l?; let mut f = l.split('\t'); let g = f.next().unwrap_or("").to_string(); let v: Vec<f64> = f.map(|x| x.parse().unwrap_or(0.0)).collect(); genes.push(g); raw.push(keep.iter().map(|&i| v[i]).collect::<Vec<f64>>()); }
    let tot: Vec<f64> = (0..keep.len()).map(|j| raw.iter().map(|r| r[j]).sum()).collect();
    let (mut g2, mut x) = (vec![], vec![]);
    for (g, r) in genes.into_iter().zip(raw) { let cpm: Vec<f64> = r.iter().zip(&tot).map(|(v, t)| v / t * 1e6).collect(); if mean(&cpm) >= MIN_CPM { g2.push(g); x.push(cpm.iter().map(|v| (v + 1.0).log2()).collect()) } }
    Ok(Pb { samples: keep.iter().map(|&i| all[i].clone()).collect(), genes: g2, x })
}
#[derive(Clone)] struct Chg { lfc: f64, t: f64, p: f64 }
fn one_sample(d: &[f64]) -> Option<Chg> { if d.len() < 3 { return None } let n = d.len() as f64; let (m, v) = (mean(d), var(d)); let t = m / (v / n).sqrt(); if !t.is_finite() { return None } Some(Chg { lfc: m, t, p: t_two_sided(t, n - 1.0) }) }
fn welch(a: &[f64], b: &[f64]) -> Option<Chg> { if a.len() < 3 || b.len() < 3 { return None } let (na, nb) = (a.len() as f64, b.len() as f64); let (va, vb) = (var(a) / na, var(b) / nb);
    let t = (mean(a) - mean(b)) / (va + vb).sqrt(); if !t.is_finite() { return None } let df = (va + vb).powi(2) / (va * va / (na - 1.0) + vb * vb / (nb - 1.0)); Some(Chg { lfc: mean(a) - mean(b), t, p: t_two_sided(t, df) }) }

type Sig = HashMap<u8, HashMap<String, (u32, u32, f64, u32)>>; // dose index -> gene -> up, down, median log2FC, cell lines
fn load_sig(path: &str, drug: &str, sig: &mut Sig) -> Result<()> {
    for l in gz_lines(path)?.lines().skip(1) { let l = l?; let f: Vec<&str> = l.split('\t').collect(); if f[0] != drug { continue }
        let d: f64 = f[1].parse()?; let k = if d < 0.2 { 0 } else if d < 2.0 { 1 } else { 2 };
        sig.entry(k).or_default().insert(f[2].to_string(), (f[4].parse()?, f[5].parse()?, f[6].parse()?, f[3].parse()?)); }
    Ok(())
}
/// Agreement of the atlas signatures of one compound with one patient comparison.
fn agreement(sig: &Sig, g: &BTreeMap<String, Chg>, hall: &BTreeMap<String, HashSet<String>>, with_mid: bool) -> Value {
    let de = |k: u8| -> HashSet<String> { sig.get(&k).map(|m| m.iter().filter(|(_, v)| (v.0 >= 3 && v.0 > v.1) || (v.1 >= 3 && v.1 > v.0)).map(|(x, _)| x.clone()).collect()).unwrap_or_default() };
    let (s0, s1, s2) = (de(0), de(1), de(2)); let hi: HashSet<String> = s2.iter().filter(|x| !s0.contains(*x) && !s1.contains(*x)).cloned().collect();
    let mut r = Map::new(); r.insert("n_genes".into(), json!(g.len())); r.insert("genes_changed_p_below_0.05".into(), json!(g.values().filter(|v| v.p < 0.05).count()));
    let mut cnt: Vec<[u64; 2]> = vec![]; let mut sets: Vec<(&str, &HashSet<String>, u8)> = vec![("low", &s0, 0)]; if with_mid { sets.push(("mid", &s1, 1)) } sets.push(("high_only", &hi, 2));
    for (nm, set, k) in &sets { let m = &sig[k]; let gs: Vec<&String> = set.iter().filter(|x| g.contains_key(*x)).collect();
        let ok = |x: &String| (m[x].2 > 0.0) == (g[x].lfc > 0.0) && m[x].2 != 0.0 && g[x].lfc != 0.0;
        let (a, n) = (gs.iter().filter(|x| ok(x)).count() as u64, gs.len() as u64); let ci = wilson(a as f64, n as f64); cnt.push([a, n - a]);
        let gc: Vec<&&String> = gs.iter().filter(|x| g[**x].p < 0.05).collect(); let ac = gc.iter().filter(|x| ok(x)).count() as f64; let cc = wilson(ac, gc.len() as f64);
        r.insert(nm.to_string(), json!({"n_signature": set.len(), "n_measured": n, "agree": a, "rate": jnum(a as f64 / n as f64), "ci": [jnum(ci.0), jnum(ci.1)], "present": ci.0 > 0.5,
            "patient_changed": {"n": gc.len(), "agree": ac, "rate": jnum(ac / gc.len() as f64), "ci": [jnum(cc.0), jnum(cc.1)]}})); }
    let (a, b) = (cnt[0], cnt[cnt.len() - 1]); r.insert("low_vs_high_only_fisher_p".into(), jnum(fisher(a[0], a[1], b[0], b[1]).1));
    for (k, nm) in [(0u8, "spearman_0.05uM"), (1, "spearman_0.5uM"), (2, "spearman_5uM")] { if let Some(m) = sig.get(&k) { let gs: Vec<&String> = m.keys().filter(|x| g.contains_key(*x)).collect();
        r.insert(nm.into(), jnum(spearman(&gs.iter().map(|x| m[*x].2).collect::<Vec<_>>(), &gs.iter().map(|x| g[*x].lfc).collect::<Vec<_>>()).0)); } }
    // gene set scores: mean over set genes minus mean over all genes (in vitro median log2FC of genes measured in >= 10 cell lines, patient t statistic)
    let mut hs = Map::new();
    for set in ["E2F_TARGETS", "G2M_CHECKPOINT"] { if let Some(genes) = hall.get(set) { let mut e = Map::new();
        let all: Vec<f64> = g.values().map(|v| v.t).collect(); let ins: Vec<f64> = g.iter().filter(|(x, _)| genes.contains(*x)).map(|(_, v)| v.t).collect(); e.insert("patient".into(), jnum(mean(&ins) - mean(&all))); e.insert("n_genes".into(), json!(ins.len()));
        for (k, nm) in [(0u8, "in_vitro_0.05uM"), (1, "in_vitro_0.5uM"), (2, "in_vitro_5uM")] { if let Some(m) = sig.get(&k) { let p: Vec<(&String, f64)> = m.iter().filter(|(x, v)| v.3 >= 10 && g.contains_key(*x)).map(|(x, v)| (x, v.2)).collect();
            let ins: Vec<f64> = p.iter().filter(|(x, _)| genes.contains(*x)).map(|x| x.1).collect(); e.insert(nm.into(), jnum(mean(&ins) - mean(&p.iter().map(|x| x.1).collect::<Vec<_>>()))); } }
        hs.insert(set.into(), Value::Object(e)); } }
    r.insert("gene_sets".into(), Value::Object(hs)); Value::Object(r)
}

fn compare(inp: &str, ribo: &str, fel: &str, osi: &str, rj: &str) -> Result<()> {
    let mut hall: BTreeMap<String, HashSet<String>> = BTreeMap::new();
    for l in std::fs::read_to_string(format!("{inp}/gene_sets/hallmark.tsv"))?.lines() { if let Some((s, g)) = l.split_once('\t') { hall.entry(s.trim_start_matches("HALLMARK_").to_string()).or_default().insert(g.to_string()); } }
    let mut res = Map::new();
    { // FELINE
        let pb = load_sums(fel)?; let mut arm: HashMap<String, String> = HashMap::new();
        for l in std::fs::read_to_string(format!("{inp}/processed/single_cell/feline_samples.tsv"))?.lines().skip(1) { let f: Vec<&str> = l.split('\t').collect(); arm.insert(f[1].to_string(), f[2].to_string()); }
        let pos: HashMap<&str, usize> = pb.samples.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
        let mut pts: Vec<(String, usize, usize)> = arm.keys().filter_map(|p| Some((p.clone(), *pos.get(format!("{p}_S").as_str())?, *pos.get(format!("{p}_M").as_str())?))).collect(); pts.sort();
        let (bc, a): (Vec<&(String, usize, usize)>, Vec<&(String, usize, usize)>) = pts.iter().partition(|p| arm[&p.0] != "A");
        let diff = |row: &Vec<f64>, v: &[&(String, usize, usize)]| -> Vec<f64> { v.iter().map(|p| row[p.2] - row[p.1]).collect() };
        let (mut prim, mut sec, mut let_) = (BTreeMap::new(), BTreeMap::new(), BTreeMap::new()); let mut dbc: Vec<(String, Vec<f64>)> = vec![];
        for (g, row) in pb.genes.iter().zip(&pb.x) { let (x, y) = (diff(row, &bc), diff(row, &a));
            if let Some(c) = welch(&x, &y) { prim.insert(g.clone(), c); } if let Some(c) = one_sample(&x) { sec.insert(g.clone(), c); } if let Some(c) = one_sample(&y) { let_.insert(g.clone(), c); }
            dbc.push((g.clone(), x)); }
        let mut sig: Sig = HashMap::new(); load_sig(ribo, "Ribociclib", &mut sig)?;
        let mut e = Map::new(); e.insert("series".into(), json!("GSE158724")); e.insert("pairs_ribociclib_arms".into(), json!(bc.len())); e.insert("pairs_letrozole_arm".into(), json!(a.len()));
        e.insert("primary_ribociclib_minus_letrozole".into(), agreement(&sig, &prim, &hall, true)); e.insert("secondary_within_ribociclib_arms".into(), agreement(&sig, &sec, &hall, true)); e.insert("letrozole_arm".into(), agreement(&sig, &let_, &hall, true));
        // agreement between random halves of the patients of the ribociclib arms
        let n = bc.len(); let mut rng = Rng(12); let mut rates = vec![];
        for _ in 0..200 { let h1 = rng.sample(n, n / 2); let in1: HashSet<usize> = h1.iter().cloned().collect(); let h2: Vec<usize> = (0..n).filter(|i| !in1.contains(i)).collect();
            let mut top: Vec<(f64, f64, usize)> = dbc.iter().enumerate().filter_map(|(i, (_, d))| one_sample(&h1.iter().map(|&j| d[j]).collect::<Vec<_>>()).map(|c| (c.p, c.lfc, i))).filter(|x| x.0.is_finite()).collect();
            top.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap()); top.truncate(500);
            let ag = top.iter().filter(|(_, l, i)| (*l > 0.0) == (mean(&h2.iter().map(|&j| dbc[*i].1[j]).collect::<Vec<_>>()) > 0.0)).count() as f64; rates.push(ag / top.len() as f64); }
        e.insert("between_halves".into(), json!({"mean_rate": mean(&rates), "q025": quantile(&rates, 0.025), "q975": quantile(&rates, 0.975)}));
        let p = &e["primary_ribociclib_minus_letrozole"]; let any = ["low", "mid", "high_only"].iter().any(|k| p[*k]["present"].as_bool().unwrap_or(false));
        e.insert("agreement_in_malignant_cells_present".into(), json!(any)); res.insert("ribociclib".into(), Value::Object(e));
    }
    if osi != "-" { // osimertinib, between patients
        let pb = load_sums(osi)?; let (mrd, tn): (Vec<usize>, Vec<usize>) = (0..pb.samples.len()).partition(|&i| pb.samples[i].starts_with("MRD"));
        let mut g = BTreeMap::new();
        for (gn, row) in pb.genes.iter().zip(&pb.x) { if let Some(c) = welch(&mrd.iter().map(|&i| row[i]).collect::<Vec<_>>(), &tn.iter().map(|&i| row[i]).collect::<Vec<_>>()) { g.insert(gn.clone(), c); } }
        let mut sig: Sig = HashMap::new(); load_sig(&format!("{inp}/processed/tahoe_signatures.tsv.gz"), "Osimertinib (mesylate)", &mut sig)?;
        res.insert("osimertinib".into(), json!({"series": "GSE292469", "samples_residual_disease": mrd.len(), "samples_treatment_naive": tn.len(), "residual_disease_minus_treatment_naive": agreement(&sig, &g, &hall, false)}));
    }
    std::fs::write(rj, serde_json::to_string_pretty(&Value::Object(res))?)?; Ok(())
}

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("signature") => signature(&a[2], &a[3], &a[4]),
        Some("feline") => feline(&a[2], &a[3]),
        Some("tenx") => tenx(&a[2], &a[3], &a[4]),
        Some("compare") => compare(&a[2], &a[3], &a[4], &a[5], &a[6]),
        _ => Err(anyhow!("usage: sc_cohorts signature|feline|tenx|compare ...")),
    }
}
