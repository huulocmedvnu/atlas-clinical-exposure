//! Shared statistics, input readers and small helpers for the analysis programs.
//!
//! Author: Thai Huu Loc Tran
use statrs::distribution::{ContinuousCDF, Normal, StudentsT};
use statrs::function::gamma::ln_gamma;

pub fn median(v: &mut Vec<f64>) -> f64 {
    v.retain(|x| x.is_finite());
    if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}
pub fn med(v: &[f64]) -> f64 { median(&mut v.to_vec()) }
/// Linear-interpolation quantile (numpy default).
pub fn quantile(v: &[f64], q: f64) -> f64 {
    let mut v: Vec<f64> = v.iter().cloned().filter(|x| x.is_finite()).collect(); if v.is_empty() { return f64::NAN; }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap()); let h = q * (v.len() - 1) as f64; let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
    v[lo] + (h - lo as f64) * (v[hi] - v[lo])
}
pub fn norm_sf(z: f64) -> f64 { Normal::new(0.0, 1.0).unwrap().sf(z) }
pub fn t_two_sided(t: f64, df: f64) -> f64 { if !t.is_finite() || df <= 0.0 { return f64::NAN; } 2.0 * StudentsT::new(0.0, 1.0, df).unwrap().sf(t.abs()) }
/// Average ranks (1-based) and the tie term sum(t^3 - t).
pub fn ranks(v: &[f64]) -> (Vec<f64>, f64) {
    let mut ix: Vec<usize> = (0..v.len()).collect(); ix.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
    let mut r = vec![0.0; v.len()]; let (mut i, mut tie) = (0, 0.0);
    while i < ix.len() { let mut j = i; while j + 1 < ix.len() && v[ix[j + 1]] == v[ix[i]] { j += 1 }
        let (m, t) = ((i + j) as f64 / 2.0 + 1.0, (j - i + 1) as f64); tie += t * t * t - t; for k in i..=j { r[ix[k]] = m } i = j + 1 }
    (r, tie)
}
pub fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64; let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for i in 0..a.len() { let (x, y) = (a[i] - ma, b[i] - mb); sab += x * y; saa += x * x; sbb += y * y }
    sab / (saa * sbb).sqrt()
}
/// Spearman correlation with the t-approximation P value.
pub fn spearman(a: &[f64], b: &[f64]) -> (f64, f64) {
    let r = pearson(&ranks(a).0, &ranks(b).0); let n = a.len() as f64;
    (r, t_two_sided(r * ((n - 2.0) / (1.0 - r * r)).sqrt(), n - 2.0))
}
/// Wilcoxon signed-rank test on paired differences x - y (zeros dropped, normal approximation with tie correction).
/// Returns (P two-sided, P for x > y, P for x < y).
pub fn wilcoxon(x: &[f64], y: &[f64]) -> (f64, f64, f64) {
    let d: Vec<f64> = x.iter().zip(y).map(|(a, b)| a - b).filter(|v| v.is_finite() && *v != 0.0).collect(); let n = d.len() as f64;
    if d.len() < 6 { return (f64::NAN, f64::NAN, f64::NAN); }
    let (r, tie) = ranks(&d.iter().map(|v| v.abs()).collect::<Vec<_>>());
    let wp: f64 = d.iter().zip(&r).filter(|(v, _)| **v > 0.0).map(|(_, r)| *r).sum();
    let z = (wp - n * (n + 1.0) / 4.0) / (n * (n + 1.0) * (2.0 * n + 1.0) / 24.0 - tie / 48.0).sqrt();
    (2.0 * norm_sf(z.abs()), norm_sf(z), norm_sf(-z))
}
/// Mann-Whitney U (two-sided, normal approximation with tie and continuity correction). Returns (AUROC of x over y, P).
pub fn mann_whitney(x: &[f64], y: &[f64]) -> (f64, f64) {
    let (n1, n2) = (x.len() as f64, y.len() as f64); if x.is_empty() || y.is_empty() { return (f64::NAN, f64::NAN); }
    let (r, tie) = ranks(&[x, y].concat()); let u = r[..x.len()].iter().sum::<f64>() - n1 * (n1 + 1.0) / 2.0; let n = n1 + n2;
    let sd = (n1 * n2 / 12.0 * ((n + 1.0) - tie / (n * (n - 1.0)))).sqrt();
    let z = ((u - n1 * n2 / 2.0).abs() - 0.5).max(0.0) / sd;
    (u / (n1 * n2), 2.0 * norm_sf(z))
}
/// Two-sample Kolmogorov-Smirnov test (asymptotic distribution). Returns (D, P).
pub fn ks_2samp(x: &[f64], y: &[f64]) -> (f64, f64) {
    let (mut a, mut b) = (x.to_vec(), y.to_vec()); a.sort_by(|p, q| p.partial_cmp(q).unwrap()); b.sort_by(|p, q| p.partial_cmp(q).unwrap());
    let (mut i, mut j, mut d) = (0, 0, 0.0f64);
    while i < a.len() && j < b.len() { let v = a[i].min(b[j]); while i < a.len() && a[i] <= v { i += 1 } while j < b.len() && b[j] <= v { j += 1 }
        d = d.max((i as f64 / a.len() as f64 - j as f64 / b.len() as f64).abs()) }
    let en = (a.len() * b.len()) as f64 / (a.len() + b.len()) as f64; let l = (en.sqrt() + 0.12 + 0.11 / en.sqrt()) * d;
    let p: f64 = (1..101).map(|k| { let k = k as f64; 2.0 * (-1f64).powf(k - 1.0) * (-2.0 * k * k * l * l).exp() }).sum();
    (d, p.clamp(0.0, 1.0))
}
fn ln_choose(n: f64, k: f64) -> f64 { ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0) }
/// Fisher exact test for [[a, b], [c, d]]. Returns (sample odds ratio, P two-sided, P for odds ratio > 1).
pub fn fisher(a: u64, b: u64, c: u64, d: u64) -> (f64, f64, f64) {
    let (r1, r2, c1) = ((a + b) as f64, (c + d) as f64, (a + c) as f64); let n = r1 + r2;
    let lp = |k: f64| ln_choose(r1, k) + ln_choose(r2, c1 - k) - ln_choose(n, c1);
    let (lo, hi) = ((c1 - r2).max(0.0) as u64, r1.min(c1) as u64); let obs = lp(a as f64);
    let (mut two, mut gt) = (0.0, 0.0);
    for k in lo..=hi { let l = lp(k as f64); if l <= obs + 1e-7 { two += l.exp() } if k >= a { gt += l.exp() } }
    ((a as f64 * d as f64) / (b as f64 * c as f64), two.min(1.0), gt.min(1.0))
}
/// One-sided binomial test, P(X >= k) under p = 0.5.
pub fn binom_greater(k: u64, n: u64) -> f64 {
    (k..=n).map(|i| (ln_choose(n as f64, i as f64) + n as f64 * 0.5f64.ln()).exp()).sum::<f64>().min(1.0)
}
/// Benjamini-Hochberg adjusted P values; NaN entries stay NaN and are not counted.
pub fn bh(p: &[f64]) -> Vec<f64> {
    let mut ix: Vec<usize> = (0..p.len()).filter(|&i| p[i].is_finite()).collect(); ix.sort_by(|&a, &b| p[a].partial_cmp(&p[b]).unwrap());
    let m = ix.len() as f64; let mut out = vec![f64::NAN; p.len()]; let mut prev = 1.0f64;
    for (rank, &i) in ix.iter().enumerate().rev() { prev = prev.min(p[i] * m / (rank + 1) as f64); out[i] = prev }
    out
}
pub fn auroc(w: &[f64], b: &[f64]) -> f64 { mann_whitney(w, b).0 }
/// Small deterministic generator (SplitMix64) for sampling without replacement.
pub struct Rng(pub u64);
impl Rng {
    pub fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z = self.0; z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }
    pub fn sample(&mut self, n: usize, k: usize) -> Vec<usize> { let mut v: Vec<usize> = (0..n).collect(); for i in 0..k.min(n) { let j = i + (self.next() % (n - i) as u64) as usize; v.swap(i, j) } v.truncate(k); v }
}
pub fn jnum(x: f64) -> serde_json::Value { if x.is_finite() { serde_json::json!(x) } else { serde_json::Value::Null } }

// ------------------------------------------------------------------ curated inputs (data/pk)
/// Reads a CSV file with optional double-quoted fields. Returns the header and the rows. Fields are not trimmed.
pub fn read_csv(path: &str) -> anyhow::Result<(Vec<String>, Vec<Vec<String>>)> {
    let txt = std::fs::read_to_string(path).map_err(|e| anyhow::anyhow!("{path}: {e}"))?; let mut rows: Vec<Vec<String>> = vec![]; let (mut row, mut f, mut q) = (vec![], String::new(), false);
    let mut it = txt.chars().peekable();
    while let Some(c) = it.next() {
        if q { if c == '"' { if it.peek() == Some(&'"') { f.push('"'); it.next(); } else { q = false } } else { f.push(c) } }
        else { match c { '"' => q = true, ',' => row.push(std::mem::take(&mut f)), '\r' => {}, '\n' => { row.push(std::mem::take(&mut f)); rows.push(std::mem::take(&mut row)) }, _ => f.push(c) } }
    }
    if !f.is_empty() || !row.is_empty() { row.push(f); rows.push(row) }
    let h = rows.remove(0); Ok((h, rows))
}
/// One Tahoe-100M compound entry. `src` is "liston" (Liston and Davis 2017), "label" (product label) or "none".
pub struct Compound { pub name: String, pub generic: String, pub moa: String, pub inert: bool, pub full: bool, pub moa_use: bool, pub cmax: f64, pub cmax_u: f64, pub src: String }
/// Tahoe-100M compounds in the order of `compound_moa_annotations.csv`, joined with the two pharmacokinetic tables.
pub fn load_tahoe_compounds(data: &str) -> anyhow::Result<Vec<Compound>> {
    use std::collections::HashMap;
    let col = |h: &[String], n: &str| h.iter().position(|c| c == n).ok_or(anyhow::anyhow!("column {n} missing"));
    let num = |s: &str| s.parse::<f64>().unwrap_or(f64::NAN);
    let mut pk: HashMap<String, (f64, f64, &str)> = HashMap::new();
    let (h, r) = read_csv(&format!("{data}/pk/liston_davis_2017_pk.csv"))?; let (ca, cc, ct, cu) = (col(&h, "atlas")?, col(&h, "compound")?, col(&h, "cmax_total_uM")?, col(&h, "cmax_unbound_uM")?);
    for x in r.iter().filter(|x| x[ca] == "Tahoe-100M") { pk.insert(x[cc].clone(), (num(&x[ct]), num(&x[cu]), "liston")); }
    let (h, r) = read_csv(&format!("{data}/pk/fda_ema_label_extension.csv"))?; let (cc, ct, cu) = (col(&h, "compound")?, col(&h, "cmax_total_uM")?, col(&h, "cmax_unbound_uM")?);
    for x in &r { pk.insert(x[cc].clone(), (num(&x[ct]), num(&x[cu]), "label")); }
    let (h, r) = read_csv(&format!("{data}/pk/compound_moa_annotations.csv"))?;
    let (cc, cg, cm, ci, cf, cs) = (col(&h, "compound")?, col(&h, "generic_name")?, col(&h, "mechanism_class")?, col(&h, "inert_reference")?, col(&h, "full_profile")?, col(&h, "mechanism_analysis")?);
    Ok(r.iter().map(|x| { let (cmax, cmax_u, src) = pk.get(&x[cc]).cloned().unwrap_or((f64::NAN, f64::NAN, "none"));
        Compound { name: x[cc].clone(), generic: x[cg].clone(), moa: x[cm].clone(), inert: x[ci] == "1", full: x[cf] == "1", moa_use: x[cs] == "1", cmax, cmax_u, src: src.to_string() } }).collect())
}
/// sci-Plex 3 compounds with pharmacokinetic data: name as in the atlas, generic name, total and unbound Cmax (uM).
pub fn load_sciplex_pk(data: &str) -> anyhow::Result<Vec<(String, String, f64, f64)>> {
    let (h, r) = read_csv(&format!("{data}/pk/liston_davis_2017_pk.csv"))?; let p = |n: &str| h.iter().position(|c| c == n).unwrap();
    Ok(r.iter().filter(|x| x[p("atlas")] == "sci-Plex 3").map(|x| (x[p("compound")].clone(), x[p("generic_name")].clone(), x[p("cmax_total_uM")].parse().unwrap(), x[p("cmax_unbound_uM")].parse().unwrap())).collect())
}

/// Exposure tier (1 to 4) and weight of a nominal dose, from total and unbound Cmax (all in uM).
pub fn exposure_tier(dose: f64, cmax: f64, cmax_u: f64) -> (u8, f64) {
    let tol = 1.0 + 1e-9; let fu = (cmax_u / cmax).min(1.0); let fum = 1.0 / (1.0 + 0.1 * (1.0 - fu) / fu); let r = dose / cmax_u;
    let t = if r <= tol { 1 } else if dose * fum / cmax_u <= tol { 2 } else if r <= 10.0 * tol { 3 } else { 4 };
    (t, (1.0 / r).powi(2).min(1.0))
}
