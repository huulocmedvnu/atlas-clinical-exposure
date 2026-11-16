//! Step 4. Agreement of Tahoe-100M signatures with paired on-treatment versus pre-treatment change in patient tumours.
//! Usage: patient_biopsies <GEO dir> <Tahoe signature table (.tsv or .tsv.gz)> <data dir> <out dir> <results json>
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use atlas_exposure::*;
use flate2::read::MultiGzDecoder;
use serde_json::{json, Map, Value};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufWriter, Read, Write};

fn gz(path: &str) -> Result<String> { let mut s = vec![]; MultiGzDecoder::new(std::fs::File::open(path)?).read_to_end(&mut s)?; Ok(String::from_utf8_lossy(&s).into_owned()) }
fn unq(s: &str) -> String { s.trim().trim_matches('"').to_string() }
fn num(s: &str) -> f64 { unq(s).parse().unwrap_or(f64::NAN) }
struct Mat { ids: Vec<String>, cols: Vec<String>, x: Vec<Vec<f64>> }

/// GEO series matrix: sample header fields and the expression table.
fn series(path: &str) -> Result<(HashMap<String, Vec<Vec<String>>>, Mat)> {
    let txt = gz(path)?; let mut hdr: HashMap<String, Vec<Vec<String>>> = HashMap::new(); let mut m = Mat { ids: vec![], cols: vec![], x: vec![] }; let mut inside = false;
    for l in txt.lines() {
        if l.starts_with("!series_matrix_table_begin") { inside = true; continue }
        if l.starts_with("!series_matrix_table_end") { break }
        let f: Vec<&str> = l.split('\t').collect();
        if inside { if m.cols.is_empty() { m.cols = f[1..].iter().map(|s| unq(s)).collect() } else { m.ids.push(unq(f[0])); m.x.push(f[1..].iter().map(|s| num(s)).collect()) } }
        else if l.starts_with("!Sample_") { hdr.entry(f[0].to_string()).or_default().push(f[1..].iter().map(|s| unq(s)).collect()) }
    }
    Ok((hdr, m))
}
/// GEO platform annotation: probe identifier to a single gene symbol.
fn annot(path: &str) -> Result<HashMap<String, String>> {
    let txt = gz(path)?; let mut it = txt.lines().skip_while(|l| !l.starts_with("ID\t"));
    let h: Vec<&str> = it.next().ok_or(anyhow!("{path}: no table"))?.split('\t').collect(); let k = h.iter().position(|c| *c == "Gene symbol").ok_or(anyhow!("no symbol column"))?;
    Ok(it.filter(|l| !l.starts_with('!')).filter_map(|l| { let f: Vec<&str> = l.split('\t').collect(); if f.len() <= k || f[k].is_empty() || f[k].contains("///") { None } else { Some((f[0].to_string(), f[k].to_string())) } }).collect())
}
fn chars(hdr: &HashMap<String, Vec<Vec<String>>>, key: &str) -> Result<Vec<String>> {
    for v in &hdr["!Sample_characteristics_ch1"] { if v[0].to_lowercase().starts_with(key) { return Ok(v.iter().map(|x| x.splitn(2, ':').nth(1).unwrap_or("").trim().to_string()).collect()) } }
    Err(anyhow!("no characteristic {key}"))
}
fn nanmean(v: &[f64]) -> f64 { let f: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect(); if f.is_empty() { f64::NAN } else { f.iter().sum::<f64>() / f.len() as f64 } }

struct Gene { lfc: f64, t: f64, p: f64 }
/// Mean within-patient difference (post minus pre) per gene with a one-sample t test.
fn paired(m: &Mat, pid: &[String], tp: &[String], pre: &str, post: &str, map: Option<&HashMap<String, String>>, log: bool) -> (BTreeMap<String, Gene>, usize) {
    let tr = |v: f64| if log { (v.max(0.0) + 1.0).log2() } else { v };
    let mut rows: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (r, id) in m.ids.iter().enumerate() { match map { Some(mp) => if let Some(g) = mp.get(id) { rows.entry(g.clone()).or_default().push(r) }, None => rows.entry(id.clone()).or_default().push(r) } }
    let col = |p: &str, t: &str| -> Vec<usize> { (0..m.cols.len()).filter(|&c| pid[c] == p && tp[c] == t).collect() };
    let mut pts: Vec<&String> = pid.iter().collect::<HashSet<_>>().into_iter().filter(|p| !col(p, pre).is_empty() && !col(p, post).is_empty()).collect(); pts.sort();
    let pc: Vec<(Vec<usize>, Vec<usize>)> = pts.iter().map(|p| (col(p, pre), col(p, post))).collect(); let need = 3.max(pts.len() / 2);
    let mut out = BTreeMap::new();
    for (g, rs) in rows {
        let gm = |cols: &Vec<usize>| nanmean(&cols.iter().map(|&c| nanmean(&rs.iter().map(|&r| tr(m.x[r][c])).collect::<Vec<_>>())).collect::<Vec<_>>());
        let d: Vec<f64> = pc.iter().map(|(a, b)| gm(b) - gm(a)).filter(|v| v.is_finite()).collect(); if d.len() < need { continue }
        let n = d.len() as f64; let mean = d.iter().sum::<f64>() / n; let sd = (d.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        let t = mean / (sd / n.sqrt()); out.insert(g, Gene { lfc: mean, t, p: t_two_sided(t, n - 1.0) });
    }
    (out, pts.len())
}
fn wilson(k: f64, n: f64) -> (f64, f64) { if n == 0.0 { return (f64::NAN, f64::NAN) } let (z, p) = (1.959964, k / n); let d = 1.0 + z * z / n;
    let (c, h) = ((p + z * z / (2.0 * n)) / d, z * (p * (1.0 - p) / n + z * z / (4.0 * n * n)).sqrt() / d); (c - h, c + h) }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect(); let (bd, sg, inp, an, rj) = (a[1].as_str(), a[2].as_str(), a[3].as_str(), a[4].as_str(), a[5].as_str()); std::fs::create_dir_all(an)?;
    let il = annot(&format!("{bd}/GPL10558.annot.gz"))?;
    let mut coh: Vec<(&str, &str, &str, BTreeMap<String, Gene>, usize)> = vec![]; // cohort, Tahoe drug, description, genes, pairs
    { let (h, m) = series(&format!("{bd}/GSE119262_series_matrix.txt.gz"))?;
      let ti: Vec<(String, String)> = h["!Sample_title"][0].iter().map(|s| { let t = s.split(": ").nth(1).unwrap_or(s); let (p, r) = t.split_once('.').unwrap(); (p.to_string(), r.split('_').next().unwrap().to_string()) }).collect();
      let (g, n) = paired(&m, &ti.iter().map(|x| x.0.clone()).collect::<Vec<_>>(), &ti.iter().map(|x| x.1.clone()).collect::<Vec<_>>(), "Pre", "Post", Some(&annot(&format!("{bd}/GPL6104.annot.gz"))?), true);
      coh.push(("everolimus", "Everolimus", "GSE119262, breast, 11 to 14 days", g, n)); }
    { let (h, m) = series(&format!("{bd}/GSE192464_series_matrix.txt.gz"))?;
      let (g, n) = paired(&m, &chars(&h, "patrient_id")?, &chars(&h, "specimen")?, "1_baseline", "2_after", Some(&il), true);
      coh.push(("afatinib", "Afatinib", "GSE192464, head and neck, 14 days", g, n)); }
    { let (h, m) = series(&format!("{bd}/GSE93204_series_matrix.txt.gz"))?; let ag = annot(&format!("{bd}/GPL6480.annot.gz"))?;
      let ti: Vec<(String, String)> = h["!Sample_title"][0].iter().map(|s| { let f: Vec<&str> = s.rsplitn(3, '_').collect(); (f[1].to_string(), f[0].to_string()) }).collect();
      let (pid, tp): (Vec<String>, Vec<String>) = ti.into_iter().unzip();
      let (g, n) = paired(&m, &pid, &tp, "C1D1", "C1D15", Some(&ag), false); coh.push(("palbociclib", "palbociclib", "GSE93204, breast, 15 days added to anastrozole", g, n));
      let (g, n) = paired(&m, &pid, &tp, "BL", "C1D1", Some(&ag), false); coh.push(("anastrozole_b", "Anastrozole", "GSE93204, breast, 4 weeks", g, n)); }
    { let txt = gz(&format!("{bd}/GSE153470_norm.txt.gz"))?; let mut it = txt.lines(); let h: Vec<&str> = it.next().unwrap().split('\t').collect();
      let keep: Vec<usize> = (1..h.len()).filter(|&c| !h[c].starts_with("Detection")).collect(); let mut m = Mat { ids: vec![], cols: keep.iter().map(|&c| unq(h[c])).collect(), x: vec![] };
      for l in it { let f: Vec<&str> = l.split('\t').collect(); if f.len() < h.len() { continue } m.ids.push(unq(f[0])); m.x.push(keep.iter().map(|&c| num(f[c])).collect()) }
      let pid: Vec<String> = m.cols.iter().map(|c| c[..c.len() - 1].to_string()).collect(); let tp: Vec<String> = m.cols.iter().map(|c| c[c.len() - 1..].to_string()).collect();
      let (g, n) = paired(&m, &pid, &tp, "A", "B", Some(&il), false); coh.push(("anastrozole", "Anastrozole", "GSE153470, breast, 2 weeks", g, n)); }
    { let (h, m) = series(&format!("{bd}/GSE58837_series_matrix.txt.gz"))?;
      let (g, n) = paired(&m, &chars(&h, "patient")?, &chars(&h, "timepoint")?, "t1", "t2", Some(&annot(&format!("{bd}/GPL6244.annot.gz"))?), false);
      coh.push(("sunitinib", "Sunitinib", "GSE58837, breast, before docetaxel", g, n)); }
    { let txt = gz(&format!("{bd}/GSE341672_counts.tsv.gz"))?; let mut it = txt.lines(); let h: Vec<&str> = it.next().unwrap().split('\t').collect();
      let mut sum: BTreeMap<String, Vec<f64>> = BTreeMap::new();
      for l in it { let f: Vec<&str> = l.split('\t').collect(); if f.len() < h.len() { continue } let e = sum.entry(unq(f[0])).or_insert_with(|| vec![0.0; h.len() - 1]); for c in 1..h.len() { let v = num(f[c]); if v.is_finite() { e[c - 1] += v } } }
      let tot: Vec<f64> = (0..h.len() - 1).map(|c| sum.values().map(|v| v[c]).sum()).collect();
      let m = Mat { ids: sum.keys().cloned().collect(), cols: h[1..].iter().map(|s| unq(s)).collect(), x: sum.values().map(|v| (0..v.len()).map(|c| v[c] / tot[c] * 1e6).collect()).collect() };
      let pid: Vec<String> = m.cols.iter().map(|c| c.split('_').next().unwrap().to_string()).collect(); let tp: Vec<String> = m.cols.iter().map(|c| c.splitn(2, '_').nth(1).unwrap_or("").to_string()).collect();
      let (g, n) = paired(&m, &pid, &tp, "B", "OT1", None, true); coh.push(("olaparib", "olaparib", "GSE341672, breast, 14 days, sequencing batch differs by time point", g, n)); }
    // Tahoe signatures
    let mut sig: HashMap<(String, u8), HashMap<String, (u32, u32, f64, u32)>> = HashMap::new();
    let want: HashSet<&str> = coh.iter().map(|c| c.1).collect();
    let sig_txt = if sg.ends_with(".gz") { gz(sg)? } else { std::fs::read_to_string(sg)? };
    for l in sig_txt.lines().skip(1) { let f: Vec<&str> = l.split('\t').collect(); if !want.contains(f[0]) { continue }
        let d: f64 = f[1].parse()?; let k = if d < 0.2 { 0 } else if d < 2.0 { 1 } else { 2 };
        sig.entry((f[0].to_string(), k)).or_default().insert(f[2].to_string(), (f[4].parse()?, f[5].parse()?, f[6].parse()?, f[3].parse()?)); }
    let stress: HashSet<String> = std::fs::read_to_string(format!("{inp}/gene_sets/stress.tsv"))?.lines().filter_map(|l| l.split('\t').nth(1).map(String::from)).collect();
    let cmax_u: HashMap<String, f64> = load_tahoe_compounds(inp)?.into_iter().map(|c| (c.name, c.cmax_u)).collect();
    let de = |d: &str, k: u8| -> HashSet<String> { sig.get(&(d.to_string(), k)).map(|m| m.iter().filter(|(_, v)| (v.0 >= 3 && v.0 > v.1) || (v.1 >= 3 && v.1 > v.0)).map(|(g, _)| g.clone()).collect()).unwrap_or_default() };
    let mut w = BufWriter::new(std::fs::File::create(format!("{an}/patient_lfc.tsv"))?); writeln!(w, "cohort\tgene\tlfc\tp")?;
    let mut hall: BTreeMap<String, HashSet<String>> = BTreeMap::new();
    for l in std::fs::read_to_string(format!("{inp}/gene_sets/hallmark.tsv"))?.lines() { if let Some((s, g)) = l.split_once('\t') { hall.entry(s.trim_start_matches("HALLMARK_").to_string()).or_default().insert(g.to_string()); } }
    let target = |d: &str| -> &[&str] { match d.to_lowercase().as_str() { "everolimus" => &["MTORC1_SIGNALING"], "palbociclib" => &["E2F_TARGETS", "G2M_CHECKPOINT"], "anastrozole" => &["ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE"], "sunitinib" => &["ANGIOGENESIS", "HYPOXIA"], _ => &[] } };
    let mut pw = BufWriter::new(std::fs::File::create(format!("{an}/pathway_scores.tsv"))?); writeln!(pw, "cohort\tset\tn_genes\tin_vitro_0.05uM\tin_vitro_5uM\tpatient")?;
    let mut topc: BTreeMap<(usize, &str, usize), [u64; 2]> = BTreeMap::new(); // (N, pool, signature) -> agree, total
    let mut direct = [[0u64; 2]; 2]; let mut pathc: BTreeMap<(&str, usize), [u64; 2]> = BTreeMap::new(); let mut rho_pos = [0usize; 2]; let mut hi_better = 0usize; let mut tgt = [[0u64; 2]; 2];
    let mut cj = Map::new(); let mut pool = [[0u64; 2]; 2]; let (mut ev, mut wins, mut cond_b) = (vec![], vec![], vec![]);
    for (name, drug, desc, g, npair) in &coh {
        for (gene, v) in g { writeln!(w, "{}\t{}\t{:.5}\t{:.4e}", name, gene, v.lfc, v.p)?; }
        let (s05, s5, s_5) = (de(drug, 0), de(drug, 2), de(drug, 1)); let supra: HashSet<String> = s5.iter().filter(|x| !s05.contains(*x) && !s_5.contains(*x)).cloned().collect();
        let cu = cmax_u[*drug]; let mut r = Map::new();
        r.insert("drug".into(), json!(drug)); r.insert("cohort".into(), json!(desc)); r.insert("n_pairs".into(), json!(npair)); r.insert("n_genes".into(), json!(g.len()));
        r.insert("unbound_cmax_uM".into(), jnum(cu)); r.insert("lowest_dose_over_unbound_cmax".into(), jnum(0.05 / cu)); r.insert("lowest_dose_reachable".into(), json!(0.05 <= cu));
        let mut cnt = [[0u64; 2]; 2];
        for (si, (nm, set, k)) in [("low", &s05, 0u8), ("high_only", &supra, 2u8)].iter().enumerate() {
            let m = &sig[&(drug.to_string(), *k)]; let gs: Vec<&String> = set.iter().filter(|x| g.contains_key(*x)).collect();
            let agree = gs.iter().filter(|x| (m[**x].2 > 0.0) == (g[**x].lfc > 0.0) && m[**x].2 != 0.0 && g[**x].lfc != 0.0).count() as u64; let n = gs.len() as u64; let ci = wilson(agree as f64, n as f64);
            cnt[si] = [agree, n - agree];
            r.insert(nm.to_string(), json!({"n_signature": set.len(), "n_measured": n, "agree": agree, "rate": jnum(agree as f64 / n as f64), "ci": [jnum(ci.0), jnum(ci.1)], "binomial_p_above_half": jnum(binom_greater(agree, n)),
                "stress_share": jnum(set.iter().filter(|x| stress.contains(*x)).count() as f64 / set.len().max(1) as f64)}));
        }
        let fz = fisher(cnt[0][0], cnt[0][1], cnt[1][0], cnt[1][1]); r.insert("low_vs_high_only".into(), json!({"odds_ratio": jnum(fz.0), "fisher_p": jnum(fz.1)}));
        for (k, nm) in [(0u8, "spearman_0.05uM"), (2u8, "spearman_5uM")] { let m = &sig[&(drug.to_string(), k)]; let gs: Vec<&String> = m.keys().filter(|x| g.contains_key(*x)).collect();
            let s = spearman(&gs.iter().map(|x| m[*x].2).collect::<Vec<_>>(), &gs.iter().map(|x| g[*x].lfc).collect::<Vec<_>>()); r.insert(nm.into(), json!({"rho": jnum(s.0), "p": jnum(s.1), "n": gs.len()})); }
        let pde: Vec<&String> = g.iter().filter(|(_, v)| v.p < 0.01).map(|(k, _)| k).collect(); let pshare = pde.iter().filter(|x| stress.contains(**x)).count() as f64 / pde.len().max(1) as f64;
        r.insert("patient_changed_genes".into(), json!(pde.len())); r.insert("patient_changed_stress_share".into(), jnum(pshare));
        let primary = !matches!(*name, "olaparib" | "anastrozole_b"); let is_direct = primary && *name != "anastrozole";
        r.insert("primary".into(), json!(primary)); r.insert("direct_target".into(), json!(is_direct));
        if is_direct { for i in 0..2 { for j in 0..2 { direct[i][j] += cnt[i][j] } } }
        // the N signature genes with the largest in vitro change
        let mut tj = Map::new();
        for nt in [50usize, 100] { let mut e = Map::new();
            for (si, (nm, set, k)) in [("low", &s05, 0u8), ("high_only", &supra, 2u8)].iter().enumerate() {
                let m = &sig[&(drug.to_string(), *k)]; let mut gs: Vec<&String> = set.iter().filter(|x| g.contains_key(*x)).collect();
                gs.sort_by(|a, b| m[*b].2.abs().partial_cmp(&m[*a].2.abs()).unwrap().then(a.cmp(b))); gs.truncate(nt);
                let agree = gs.iter().filter(|x| (m[**x].2 > 0.0) == (g[**x].lfc > 0.0)).count() as u64; let n = gs.len() as u64; let ci = wilson(agree as f64, n as f64);
                let sp = spearman(&gs.iter().map(|x| m[*x].2).collect::<Vec<_>>(), &gs.iter().map(|x| g[*x].lfc).collect::<Vec<_>>());
                e.insert(nm.to_string(), json!({"n": n, "agree": agree, "rate": jnum(agree as f64 / n as f64), "ci": [jnum(ci.0), jnum(ci.1)], "rho": jnum(sp.0), "rho_p": jnum(sp.1)}));
                for (pn, ok) in [("primary", primary), ("direct", is_direct)] { if ok { let c = topc.entry((nt, pn, si)).or_default(); c[0] += agree; c[1] += n } } }
            tj.insert(nt.to_string(), Value::Object(e)); }
        r.insert("top".into(), Value::Object(tj));
        // hallmark set scores: mean over set genes minus mean over all genes, in vitro (median log2FC, genes in >= 10 lines) and in patients (paired t)
        { let prof = |k: u8| -> HashMap<&String, f64> { sig[&(drug.to_string(), k)].iter().filter(|(x, v)| v.3 >= 10 && g.contains_key(*x)).map(|(x, v)| (x, v.2)).collect() };
          let (p0, p2) = (prof(0), prof(2)); let common: Vec<&String> = p0.keys().filter(|x| p2.contains_key(**x)).cloned().collect(); let nall = common.len() as f64;
          let (m0, m2, mt) = (common.iter().map(|x| p0[x]).sum::<f64>() / nall, common.iter().map(|x| p2[x]).sum::<f64>() / nall, common.iter().map(|x| g[*x].t).sum::<f64>() / nall);
          let mut sc: Vec<(&String, f64, f64, f64)> = vec![];
          for (sn, genes) in &hall { let gs: Vec<&&String> = common.iter().filter(|x| genes.contains(**x)).collect(); if gs.len() < 10 { continue } let n = gs.len() as f64;
              let v = (gs.iter().map(|x| p0[**x]).sum::<f64>() / n - m0, gs.iter().map(|x| p2[**x]).sum::<f64>() / n - m2, gs.iter().map(|x| g[**x].t).sum::<f64>() / n - mt);
              writeln!(pw, "{}\t{}\t{}\t{:.5}\t{:.5}\t{:.5}", name, sn, gs.len(), v.0, v.1, v.2)?; sc.push((sn, v.0, v.1, v.2)); }
          let pt: Vec<f64> = sc.iter().map(|x| x.3).collect(); let mut e = Map::new(); let mut rr = [0.0; 2];
          for (si, nm) in ["0.05uM", "5uM"].iter().enumerate() { let iv: Vec<f64> = sc.iter().map(|x| if si == 0 { x.1 } else { x.2 }).collect(); let sp = spearman(&iv, &pt);
              let agree = iv.iter().zip(&pt).filter(|(a, b)| (**a > 0.0) == (**b > 0.0)).count() as u64; rr[si] = sp.0;
              e.insert(nm.to_string(), json!({"n_sets": sc.len(), "rho": jnum(sp.0), "rho_p": jnum(sp.1), "agree": agree, "rate": agree as f64 / sc.len() as f64}));
              if primary { let c = pathc.entry(("primary", si)).or_default(); c[0] += agree; c[1] += sc.len() as u64; if sp.0 > 0.0 { rho_pos[si] += 1 } }
              if is_direct { let c = pathc.entry(("direct", si)).or_default(); c[0] += agree; c[1] += sc.len() as u64 } }
          if primary && rr[1] > rr[0] { hi_better += 1 }
          e.insert("rho_0.05uM_vs_5uM_in_vitro".into(), jnum(spearman(&sc.iter().map(|x| x.1).collect::<Vec<_>>(), &sc.iter().map(|x| x.2).collect::<Vec<_>>()).0));
          let mut te = Map::new();
          for t in target(drug) { if let Some(x) = sc.iter().find(|x| x.0 == t) { te.insert(t.to_string(), json!({"in_vitro_0.05uM": x.1, "in_vitro_5uM": x.2, "patient": x.3}));
              if primary { for (si, v) in [x.1, x.2].iter().enumerate() { tgt[si][0] += ((*v > 0.0) == (x.3 > 0.0)) as u64; tgt[si][1] += 1 } } } }
          e.insert("target_sets".into(), Value::Object(te)); e.insert("n_genes".into(), json!(common.len()));
          r.insert("pathway".into(), Value::Object(e)); }
        if primary { let (ra, rb) = (cnt[0][0] as f64 / (cnt[0][0] + cnt[0][1]) as f64, cnt[1][0] as f64 / (cnt[1][0] + cnt[1][1]) as f64);
            ev.push(*name); if ra > rb { wins.push(*name) } for i in 0..2 { for j in 0..2 { pool[i][j] += cnt[i][j] } }
            let ss = supra.iter().filter(|x| stress.contains(*x)).count() as f64 / supra.len().max(1) as f64; if ss > 0.0 && pshare < 0.5 * ss { cond_b.push(*name) } }
        cj.insert(name.to_string(), Value::Object(r));
    }
    let fz = fisher(pool[0][0], pool[0][1], pool[1][0], pool[1][1]);
    let (cl, ch) = (wilson(pool[0][0] as f64, (pool[0][0] + pool[0][1]) as f64), wilson(pool[1][0] as f64, (pool[1][0] + pool[1][1]) as f64));
    let cond_a = wins.len() * 2 > ev.len() && fz.2 < 0.05; let a2 = cond_b.len() * 2 > ev.len();
    let pooled = json!({"cohorts": ev, "low_higher_in": wins, "low": [pool[0][0], pool[0][0] + pool[0][1]], "high_only": [pool[1][0], pool[1][0] + pool[1][1]],
        "low_rate": pool[0][0] as f64 / (pool[0][0] + pool[0][1]) as f64, "low_ci": [cl.0, cl.1], "high_only_rate": pool[1][0] as f64 / (pool[1][0] + pool[1][1]) as f64, "high_only_ci": [ch.0, ch.1],
        "fisher_p_two_sided": jnum(fz.1), "fisher_p_low_higher": jnum(fz.2), "cohorts_with_lower_stress_share": cond_b, "prespecified_condition_met": cond_a || a2});
    drop(pw);
    let two = |a: [u64; 2], b: [u64; 2]| { let fz = fisher(a[0], a[1] - a[0], b[0], b[1] - b[0]); let (ca, cb) = (wilson(a[0] as f64, a[1] as f64), wilson(b[0] as f64, b[1] as f64));
        json!({"low": a, "low_rate": a[0] as f64 / a[1] as f64, "low_ci": [ca.0, ca.1], "high": b, "high_rate": b[0] as f64 / b[1] as f64, "high_ci": [cb.0, cb.1], "fisher_p": jnum(fz.1),
               "binomial_p_low_above_half": jnum(binom_greater(a[0], a[1])), "binomial_p_high_above_half": jnum(binom_greater(b[0], b[1]))}) };
    let mut agg = Map::new();
    agg.insert("direct_target_genes".into(), two([direct[0][0], direct[0][0] + direct[0][1]], [direct[1][0], direct[1][0] + direct[1][1]]));
    for nt in [50usize, 100] { for pn in ["primary", "direct"] { agg.insert(format!("top{nt}_{pn}"), two(topc[&(nt, pn, 0)], topc[&(nt, pn, 1)])); } }
    for pn in ["primary", "direct"] { agg.insert(format!("pathway_{pn}"), two(pathc[&(pn, 0)], pathc[&(pn, 1)])); }
    let np = ev.len(); let pp = &agg["pathway_primary"];
    let rescue = |si: usize, key: &str| pp[key].as_f64().unwrap_or(1.0) < 0.05 && rho_pos[si] * 2 > np;
    let (res_low, res_high) = (rescue(0, "binomial_p_low_above_half"), rescue(1, "binomial_p_high_above_half"));
    let outperf = hi_better * 2 > np && pp["fisher_p"].as_f64().unwrap_or(1.0) < 0.05 && pp["high_rate"].as_f64() > pp["low_rate"].as_f64();
    agg.insert("pathway_rule".into(), json!({"cohorts_rho_positive_0.05uM": rho_pos[0], "cohorts_rho_positive_5uM": rho_pos[1], "cohorts_5uM_rho_higher": hi_better, "n_primary": np,
        "rescue_0.05uM": res_low, "rescue_5uM": res_high, "5uM_outperforms": outperf, "target_sets_same_direction_0.05uM": tgt[0], "target_sets_same_direction_5uM": tgt[1]}));
    std::fs::write(rj, serde_json::to_string_pretty(&json!({"cohorts": cj, "pooled": pooled, "aggregated": agg}))?)?;
    Ok(())
}
