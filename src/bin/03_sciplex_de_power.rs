//! Step 3. sci-Plex 3: per-well sums, negative binomial differential expression, exposure metrics and power at clinical exposure.
//!   sciplex_de_power pseudobulk <h5ad> <data dir> <work dir>   reads the h5ad once and sums counts per well
//!   sciplex_de_power de <work dir> <data dir> <out dir> <results json>   tests, exposure metrics, power
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use atlas_exposure::*;
use hdf5::types::VarLenUnicode;
use rayon::prelude::*;
use serde_json::{json, Map, Value};
use statrs::function::gamma::ln_gamma;
use std::collections::{HashMap, HashSet};
use std::io::{BufWriter, Read, Write};

const DOSES: [f64; 4] = [0.01, 0.1, 1.0, 10.0];
const BASIS: [&str; 3] = ["unbound", "total", "fbs"];

struct Pk { cmax: f64, cmax_u: f64 }
impl Pk {
    fn thr_scale(&self, b: usize) -> (f64, f64) {
        let fu = (self.cmax_u / self.cmax).min(1.0);
        match b { 0 => (self.cmax_u, 1.0), 1 => (self.cmax, 1.0), _ => (self.cmax_u, 1.0 / (1.0 + 0.1 * (1.0 - fu) / fu)) }
    }
    fn reach(&self, b: usize, d: usize) -> bool { let (t, s) = self.thr_scale(b); DOSES[d] * s <= t * (1.0 + 1e-9) }
}
fn load_pk(inp: &str) -> Result<HashMap<String, Pk>> {
    Ok(load_sciplex_pk(inp)?.into_iter().map(|(d, _, cmax, cmax_u)| (d, Pk { cmax, cmax_u })).collect())
}
fn dose_idx(nm: f64) -> usize { DOSES.iter().position(|d| (d * 1000.0 - nm).abs() < 1e-6 * nm.max(1.0)).unwrap_or(9) }

// ------------------------------------------------------------------ pseudobulk
fn cat(f: &hdf5::File, name: &str) -> Result<(Vec<i32>, Vec<String>)> {
    let codes: Vec<i32> = f.dataset(&format!("{name}/codes"))?.read_raw()?;
    let cats: Vec<VarLenUnicode> = f.dataset(&format!("{name}/categories"))?.read_raw()?;
    Ok((codes, cats.into_iter().map(|s| s.as_str().to_string()).collect()))
}
fn pseudobulk(h5: &str, inp: &str, wd: &str) -> Result<()> {
    std::fs::create_dir_all(wd)?;
    let pk = load_pk(inp)?; let f = hdf5::File::open(h5)?; let t0 = std::time::Instant::now();
    let (pc, pn) = cat(&f, "obs/perturbation")?; let (lc, lnm) = cat(&f, "obs/cell_line")?; let (rc, rn) = cat(&f, "obs/replicate")?;
    let (plc, pln) = cat(&f, "obs/plate")?; let (wc, wn) = cat(&f, "obs/well")?;
    let dose: Vec<f64> = f.dataset("obs/dose_value")?.read_raw()?; let time: Vec<f64> = f.dataset("obs/time")?.read_raw()?;
    let (gc, gn) = cat(&f, "var/ensembl_id")?; let n = pc.len();
    let keep: Vec<bool> = (0..n).map(|i| time[i] == 24.0 && pc[i] >= 0 && lc[i] >= 0 && (pn[pc[i] as usize] == "control" || pk.contains_key(&pn[pc[i] as usize]))).collect();
    let key = |i: usize| format!("{}\t{}\t{}\t{}\t{}\t{}", lnm[lc[i] as usize], pn[pc[i] as usize], dose[i], rn[rc[i] as usize], pln[plc[i] as usize], wn[wc[i] as usize]);
    let mut wells: Vec<String> = (0..n).filter(|&i| keep[i]).map(key).collect::<HashSet<_>>().into_iter().collect(); wells.sort();
    let widx: HashMap<&str, usize> = wells.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let gi: Vec<i32> = (0..n).map(|i| if keep[i] { widx[key(i).as_str()] as i32 } else { -1 }).collect();
    let mut cells: Vec<Vec<usize>> = vec![vec![]; wells.len()]; for i in 0..n { if gi[i] >= 0 { cells[gi[i] as usize].push(i) } }
    // downsampling plan: a 10 uM well with more cells than the 10 nM well of the same drug, line and replicate keeps a random subset of that size
    let part = |s: &str| -> Vec<String> { s.split('\t').map(String::from).collect() };
    let mut by: HashMap<(String, String, String), HashMap<usize, usize>> = HashMap::new();
    for (w, s) in wells.iter().enumerate() { let p = part(s); if p[1] == "control" { continue } by.entry((p[0].clone(), p[1].clone(), p[3].clone())).or_default().insert(dose_idx(p[2].parse().unwrap()), w); }
    let mut drop_ds = vec![false; n]; let mut ds_wells: Vec<(usize, usize, usize)> = vec![]; let mut rng = Rng(7);
    let mut ks: Vec<_> = by.iter().collect(); ks.sort_by(|a, b| a.0.cmp(b.0));
    for (_, m) in ks { let (Some(&lo), Some(&hi)) = (m.get(&0), m.get(&3)) else { continue };
        if cells[hi].len() > cells[lo].len() { let pick: HashSet<usize> = rng.sample(cells[hi].len(), cells[lo].len()).into_iter().collect();
            for (j, &c) in cells[hi].iter().enumerate() { if !pick.contains(&j) { drop_ds[c] = true } } ds_wells.push((hi, cells[hi].len(), cells[lo].len())) } }
    eprintln!("cells {} kept {} wells {} downsampled wells {} [{:.0}s]", n, keep.iter().filter(|x| **x).count(), wells.len(), ds_wells.len(), t0.elapsed().as_secs_f64());
    let ng = gn.len(); let (mut acc, mut accd) = (vec![0f32; wells.len() * ng], vec![0f32; wells.len() * ng]);
    let ip: Vec<i64> = f.dataset("X/indptr")?.read_raw()?; let (dd, di) = (f.dataset("X/data")?, f.dataset("X/indices")?);
    let step = 20000;
    for s in (0..n).step_by(step) {
        let e = (s + step).min(n); let (a, b) = (ip[s] as usize, ip[e] as usize); if a == b { continue }
        let x: Vec<f32> = dd.read_slice_1d::<f32, _>(a..b)?.to_vec(); let j: Vec<i32> = di.read_slice_1d::<i32, _>(a..b)?.to_vec();
        for c in s..e { if gi[c] < 0 { continue } let w = gi[c] as usize * ng;
            for k in (ip[c] as usize - a)..(ip[c + 1] as usize - a) { let g = gc[j[k] as usize]; if g < 0 { continue }
                acc[w + g as usize] += x[k]; if !drop_ds[c] { accd[w + g as usize] += x[k] } } }
        if (s / step) % 10 == 0 { eprintln!("  cells {} [{:.0}s]", e, t0.elapsed().as_secs_f64()) }
    }
    let gk: Vec<usize> = (0..ng).filter(|&g| (0..wells.len()).map(|w| acc[w * ng + g] as f64).sum::<f64>() >= 50.0).collect();
    let mut w = BufWriter::new(std::fs::File::create(format!("{wd}/meta.tsv"))?);
    writeln!(w, "cell_line\tdrug\tdose_nM\treplicate\tplate\twell\tn_cells\tn_cells_downsampled")?;
    for (i, s) in wells.iter().enumerate() { let d = ds_wells.iter().find(|x| x.0 == i).map(|x| x.2).unwrap_or(cells[i].len()); writeln!(w, "{}\t{}\t{}", s, cells[i].len(), d)?; }
    std::fs::write(format!("{wd}/genes.txt"), gk.iter().map(|&g| gn[g].clone()).collect::<Vec<_>>().join("\n"))?;
    for (name, m) in [("counts.bin", &acc), ("counts_ds.bin", &accd)] {
        let mut w = BufWriter::new(std::fs::File::create(format!("{wd}/{name}"))?);
        for r in 0..wells.len() { for &g in &gk { w.write_all(&m[r * ng + g].to_le_bytes())? } } }
    eprintln!("wells {} genes {} [{:.0}s]", wells.len(), gk.len(), t0.elapsed().as_secs_f64());
    Ok(())
}

// ------------------------------------------------------------------ negative binomial model, one-way layout
fn trigamma(mut x: f64) -> f64 { let mut r = 0.0; while x < 6.0 { r += 1.0 / (x * x); x += 1.0 }
    let f = 1.0 / (x * x); r + 1.0 / x + f / 2.0 + f / x * (1.0 / 6.0 - f * (1.0 / 30.0 - f / 42.0)) }
/// Cox-Reid adjusted profile log-likelihood of the dispersion for one gene with fixed fitted means.
fn cr_ll(la: f64, y: &[f64], mu: &[f64], grp: &[usize], ngrp: usize, wbuf: &mut Vec<f64>) -> f64 {
    let a = la.exp(); let ia = 1.0 / a; let mut ll = 0.0; wbuf.clear(); wbuf.resize(ngrp, 0.0);
    for j in 0..y.len() { let m = mu[j]; if m <= 0.0 { continue }
        ll += ln_gamma(y[j] + ia) - ln_gamma(ia) - ln_gamma(y[j] + 1.0) + y[j] * (a * m / (1.0 + a * m)).ln() - ia * (1.0 + a * m).ln();
        wbuf[grp[j]] += m / (1.0 + a * m) }
    ll - 0.5 * wbuf.iter().filter(|w| **w > 0.0).map(|w| w.ln()).sum::<f64>()
}
fn maximise<F: FnMut(f64) -> f64>(mut f: F, lo: f64, hi: f64) -> f64 {
    let n = 24; let (mut bx, mut bv) = (lo, f64::NEG_INFINITY);
    for i in 0..=n { let x = lo + (hi - lo) * i as f64 / n as f64; let v = f(x); if v > bv { bv = v; bx = x } }
    let h = (hi - lo) / n as f64; let (mut a, mut b) = ((bx - h).max(lo), (bx + h).min(hi)); let gr = 0.618_033_988_75;
    let (mut c, mut d) = (b - gr * (b - a), a + gr * (b - a)); let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..40 { if fc > fd { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = f(c) } else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = f(d) } }
    0.5 * (a + b)
}
struct Fit { genes: Vec<usize>, base: Vec<f64>, a0: f64, a1: f64, lfc: Vec<Vec<f64>>, se: Vec<Vec<f64>>, padj: Vec<Vec<f64>> } // per non-reference group, per gene
/// y: counts [gene][sample]; grp: group of each sample, 0 = reference.
fn fit(y: &[Vec<f64>], grp: &[usize], ngrp: usize) -> Fit {
    let m = grp.len();
    let genes: Vec<usize> = (0..y.len()).filter(|&g| y[g].iter().filter(|v| **v >= 10.0).count() >= 2).collect();
    // median-of-ratios size factors; total-count scaling if too few genes are detected in every sample
    let lg: Vec<Option<f64>> = genes.iter().map(|&g| if y[g].iter().all(|v| *v > 0.0) { Some(y[g].iter().map(|v| v.ln()).sum::<f64>() / m as f64) } else { None }).collect();
    let nall = lg.iter().filter(|x| x.is_some()).count();
    let sf: Vec<f64> = if nall >= 50 { (0..m).map(|j| med(&genes.iter().zip(&lg).filter_map(|(&g, l)| l.map(|l| y[g][j].ln() - l)).collect::<Vec<_>>()).exp()).collect() }
        else { let t: Vec<f64> = (0..m).map(|j| genes.iter().map(|&g| y[g][j]).sum::<f64>()).collect(); let gm = (t.iter().map(|v| v.ln()).sum::<f64>() / m as f64).exp(); t.iter().map(|v| v / gm).collect() };
    let gsf: Vec<f64> = (0..ngrp).map(|k| (0..m).filter(|&j| grp[j] == k).map(|j| sf[j]).sum()).collect();
    let base: Vec<f64> = genes.iter().map(|&g| (0..m).map(|j| y[g][j] / sf[j]).sum::<f64>() / m as f64).collect();
    let pois = |g: usize| -> (Vec<f64>, Vec<f64>) { let mut s = vec![0.0; ngrp]; for j in 0..m { s[grp[j]] += y[g][j] }
        let q: Vec<f64> = (0..ngrp).map(|k| s[k] / gsf[k]).collect(); ((0..m).map(|j| sf[j] * q[grp[j]]).collect(), s) };
    let (lo, hi) = ((1e-8f64).ln(), (10.0f64).ln());
    let gw: Vec<f64> = genes.par_iter().map(|&g| { let (mu, _) = pois(g); let mut wb = vec![]; maximise(|la| cr_ll(la, &y[g], &mu, grp, ngrp, &mut wb), lo, hi) }).collect();
    // mean-dispersion trend alpha = a1 / mean + a0, fitted by iterative weighted least squares
    let use_: Vec<usize> = (0..genes.len()).filter(|&i| gw[i] > (1e-6f64).ln()).collect();
    let (mut a0, mut a1) = (0.1, 1.0);
    for _ in 0..10 {
        let (mut s00, mut s01, mut s11, mut t0, mut t1) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for &i in &use_ { let p = a0 + a1 / base[i]; let r = gw[i].exp() / p; if !(r > 1e-4 && r < 15.0) { continue }
            let (w, x) = (1.0 / (p * p), 1.0 / base[i]); s00 += w; s01 += w * x; s11 += w * x * x; t0 += w * gw[i].exp(); t1 += w * x * gw[i].exp() }
        let det = s00 * s11 - s01 * s01; if det.abs() < 1e-300 { break }
        let (n0, n1) = ((t0 * s11 - t1 * s01) / det, (t1 * s00 - t0 * s01) / det); if !(n0 > 0.0 && n1 > 0.0) { break }
        let done = ((n0 / a0).ln().abs() + (n1 / a1).ln().abs()) < 1e-6; a0 = n0; a1 = n1; if done { break }
    }
    let tr: Vec<f64> = base.iter().map(|b| (a0 + a1 / b).ln().clamp(lo, hi)).collect();
    let resid: Vec<f64> = use_.iter().map(|&i| gw[i] - tr[i]).collect(); let mr = med(&resid);
    let s_lr = 1.4826 * med(&resid.iter().map(|r| (r - mr).abs()).collect::<Vec<_>>());
    let prior = (s_lr * s_lr - trigamma((m as f64 - ngrp as f64) / 2.0)).max(0.25);
    let out: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = genes.par_iter().enumerate().map(|(i, &g)| {
        let (mu, s) = pois(g); let mut wb = vec![];
        let la = if gw[i] > tr[i] + 2.0 * s_lr { gw[i] } else { maximise(|la| cr_ll(la, &y[g], &mu, grp, ngrp, &mut wb) - (la - tr[i]).powi(2) / (2.0 * prior), lo, hi) };
        let a = la.exp();
        // group means: maximum likelihood under the final dispersion; a group with no counts gets half a count
        let (mut q, mut w) = (vec![0.0; ngrp], vec![0.0; ngrp]);
        for k in 0..ngrp { let mut qk = if s[k] > 0.0 { s[k] / gsf[k] } else { 0.5 / gsf[k] };
            if s[k] > 0.0 { for _ in 0..25 { let (mut f, mut df) = (0.0, 0.0);
                for j in 0..m { if grp[j] != k { continue } let d = 1.0 + a * sf[j] * qk; f += (y[g][j] - sf[j] * qk) / d; df -= sf[j] * (1.0 + a * y[g][j]) / (d * d) }
                let nq = qk - f / df; let nq = if nq > 0.0 && nq.is_finite() { nq } else { qk / 2.0 }; let c = (nq / qk).ln().abs(); qk = nq; if c < 1e-8 { break } } }
            q[k] = qk; w[k] = (0..m).filter(|&j| grp[j] == k).map(|j| sf[j] * qk / (1.0 + a * sf[j] * qk)).sum() }
        let (mut lf, mut se, mut p) = (vec![], vec![], vec![]);
        for k in 1..ngrp { let l = (q[k] / q[0]).ln(); let e = (1.0 / w[k] + 1.0 / w[0]).sqrt();
            lf.push(l / std::f64::consts::LN_2); se.push(e / std::f64::consts::LN_2); p.push(2.0 * norm_sf((l / e).abs())) }
        (lf, se, p)
    }).collect();
    // Benjamini-Hochberg with independent filtering on the mean of normalised counts
    let th: Vec<f64> = (0..48).map(|i| quantile(&base, i as f64 * 0.02)).collect();
    let per: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..ngrp - 1).into_par_iter().map(|k| {
        let p: Vec<f64> = out.iter().map(|o| o.2[k]).collect(); let (mut best, mut bn) = (vec![], -1i64);
        for t in &th { let pp: Vec<f64> = (0..p.len()).map(|i| if base[i] >= *t { p[i] } else { f64::NAN }).collect(); let ad = bh(&pp);
            let nr = ad.iter().filter(|v| **v < 0.05).count() as i64; if nr > bn { bn = nr; best = ad } }
        (out.iter().map(|o| o.0[k]).collect(), out.iter().map(|o| o.1[k]).collect(), best)
    }).collect();
    let mut f = Fit { genes, base, a0, a1, lfc: vec![], se: vec![], padj: vec![] };
    for (a, b, c) in per { f.lfc.push(a); f.se.push(b); f.padj.push(c) }
    f
}

fn read_bin(path: &str, rows: usize, cols: usize) -> Result<Vec<Vec<f32>>> {
    let mut b = vec![];
    if std::path::Path::new(path).exists() { std::fs::File::open(path)?.read_to_end(&mut b)?; } else { flate2::read::MultiGzDecoder::new(std::fs::File::open(format!("{path}.gz"))?).read_to_end(&mut b)?; } if b.len() != rows * cols * 4 { return Err(anyhow!("{path}: unexpected size")) }
    Ok((0..rows).map(|r| (0..cols).map(|c| { let o = (r * cols + c) * 4; f32::from_le_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]]) }).collect()).collect())
}

/// par: dispersion trend a0, a1, median cells per well of the fit, vehicle wells, median cells per vehicle well; ncell: mean cells per treated well.
struct Con { drug: String, dose: usize, line: String, lfc: Vec<f64>, se: Vec<f64>, padj: Vec<f64>, genes: Vec<usize>, base: Vec<f64>, par: [f64; 5], ncell: f64 }
impl Con {
    fn de(&self, cut: f64) -> Vec<usize> { (0..self.lfc.len()).filter(|&i| self.padj[i] < 0.05 && self.lfc[i].abs() >= cut).map(|i| self.genes[i]).collect() }
    fn l2(&self) -> f64 { (0..self.lfc.len()).map(|i| { let (b, s) = (self.lfc[i], self.se[i]); let v = b * b * b / (b * b + s * s); v * v }).sum::<f64>().sqrt() }
}

fn de(wd: &str, inp: &str, od: &str, rj: &str) -> Result<()> {
    std::fs::create_dir_all(od)?; let t0 = std::time::Instant::now();
    let pk = load_pk(inp)?;
    let meta: Vec<Vec<String>> = std::fs::read_to_string(format!("{wd}/meta.tsv"))?.lines().skip(1).map(|l| l.split('\t').map(String::from).collect()).collect();
    let ng = std::fs::read_to_string(format!("{wd}/genes.txt"))?.lines().count();
    let (cnt, cntd) = (read_bin(&format!("{wd}/counts.bin"), meta.len(), ng)?, read_bin(&format!("{wd}/counts_ds.bin"), meta.len(), ng)?);
    let mut lines: Vec<String> = meta.iter().map(|m| m[0].clone()).collect::<HashSet<_>>().into_iter().collect(); lines.sort();
    let ncell = |m: &Vec<String>| m[6].parse::<f64>().unwrap(); let dnm = |m: &Vec<String>| m[2].parse::<f64>().unwrap();
    let run = |line: &str, c: &Vec<Vec<f32>>| -> Vec<Con> {
        let ws: Vec<usize> = (0..meta.len()).filter(|&w| meta[w][0] == line && ncell(&meta[w]) >= 20.0).collect();
        let mut conds: Vec<(String, usize)> = ws.iter().filter(|&&w| meta[w][1] != "control").map(|&w| (meta[w][1].clone(), dose_idx(dnm(&meta[w])))).collect::<HashSet<_>>().into_iter().collect(); conds.sort();
        let grp: Vec<usize> = ws.iter().map(|&w| if meta[w][1] == "control" { 0 } else { 1 + conds.iter().position(|c| c.0 == meta[w][1] && c.1 == dose_idx(dnm(&meta[w]))).unwrap() }).collect();
        let y: Vec<Vec<f64>> = (0..ng).into_par_iter().map(|g| ws.iter().map(|&w| c[w][g] as f64).collect()).collect();
        let mut f = fit(&y, &grp, conds.len() + 1);
        let veh: Vec<f64> = ws.iter().filter(|&&w| meta[w][1] == "control").map(|&w| ncell(&meta[w])).collect();
        let par = [f.a0, f.a1, med(&ws.iter().map(|&w| ncell(&meta[w])).collect::<Vec<_>>()), veh.len() as f64, med(&veh)];
        conds.into_iter().enumerate().map(|(k, (d, ds))| { let nc: Vec<f64> = ws.iter().zip(&grp).filter(|(_, g)| **g == k + 1).map(|(&w, _)| ncell(&meta[w])).collect();
            Con { drug: d, dose: ds, line: line.to_string(), lfc: std::mem::take(&mut f.lfc[k]), se: std::mem::take(&mut f.se[k]), padj: std::mem::take(&mut f.padj[k]), genes: f.genes.clone(), base: f.base.clone(), par, ncell: nc.iter().sum::<f64>() / nc.len() as f64 } }).collect()
    };
    let cons: Vec<Con> = lines.iter().flat_map(|l| run(l, &cnt)).collect();
    eprintln!("contrasts {} [{:.0}s]", cons.len(), t0.elapsed().as_secs_f64());
    let ds_lines: HashSet<&str> = meta.iter().filter(|m| m[6] != m[7]).map(|m| m[0].as_str()).collect();
    let cons_ds: Vec<Con> = lines.iter().filter(|l| ds_lines.contains(l.as_str())).flat_map(|l| run(l, &cntd)).collect();
    // vehicle-split noise: two vehicle wells against the remaining vehicle wells, ten draws per line
    let mut rng = Rng(1); let mut noise: Vec<Con> = vec![];
    for line in &lines {
        let v: Vec<usize> = (0..meta.len()).filter(|&w| meta[w][0] == *line && meta[w][1] == "control" && ncell(&meta[w]) >= 20.0).collect();
        for _ in 0..10 { let fake: HashSet<usize> = rng.sample(v.len(), 2).into_iter().collect(); let grp: Vec<usize> = (0..v.len()).map(|j| fake.contains(&j) as usize).collect();
            let y: Vec<Vec<f64>> = (0..ng).map(|g| v.iter().map(|&w| cnt[w][g] as f64).collect()).collect(); let mut f = fit(&y, &grp, 2);
            noise.push(Con { drug: "vehicle_split".into(), dose: 9, line: line.clone(), lfc: std::mem::take(&mut f.lfc[0]), se: std::mem::take(&mut f.se[0]), padj: std::mem::take(&mut f.padj[0]), genes: f.genes, base: vec![], par: [0.0; 5], ncell: 0.0 }) }
    }
    eprintln!("noise and downsampled fits done [{:.0}s]", t0.elapsed().as_secs_f64());
    let mut res = Map::new();
    // cell numbers
    { let mut cj = Map::new();
      for d in 0..4 { let v: Vec<f64> = meta.iter().filter(|m| m[1] != "control" && dose_idx(dnm(m)) == d).map(ncell).collect(); cj.insert(format!("{}", DOSES[d]), summ(&v)); }
      cj.insert("vehicle".into(), summ(&meta.iter().filter(|m| m[1] == "control").map(ncell).collect::<Vec<_>>()));
      cj.insert("wells_10uM_above_10nM".into(), json!(meta.iter().filter(|m| m[6] != m[7]).count()));
      cj.insert("wells_10uM".into(), json!(meta.iter().filter(|m| m[1] != "control" && dose_idx(dnm(m)) == 3).count()));
      res.insert("cells_per_well".into(), Value::Object(cj)); }
    res.insert("data".into(), json!({"wells": meta.len(), "genes": ng, "contrasts": cons.len(), "compounds": cons.iter().map(|c| c.drug.as_str()).collect::<HashSet<_>>().len(), "cell_lines": lines.len()}));
    // magnitude: top 500 genes by variance of log2FC across contrasts
    let mut acc: HashMap<usize, (f64, f64, f64)> = HashMap::new();
    for c in &cons { for i in 0..c.lfc.len() { let e = acc.entry(c.genes[i]).or_default(); e.0 += 1.0; e.1 += c.lfc[i]; e.2 += c.lfc[i] * c.lfc[i] } }
    let mut gv: Vec<(f64, usize)> = acc.iter().filter(|(_, e)| e.0 >= 0.8 * cons.len() as f64 / lines.len() as f64).map(|(g, e)| ((e.2 - e.1 * e.1 / e.0) / (e.0 - 1.0), *g)).collect();
    gv.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap()); let top: HashSet<usize> = gv.iter().take(500).map(|x| x.1).collect();
    let mag1 = |c: &Con| med(&(0..c.lfc.len()).filter(|&i| top.contains(&c.genes[i])).map(|i| c.lfc[i].abs()).collect::<Vec<_>>());
    let mut metrics: Vec<Value> = vec![]; let mut w = BufWriter::new(std::fs::File::create(format!("{od}/sciplex_contrasts.tsv"))?);
    writeln!(w, "drug\tdose_uM\tline\tcmax_uM\tcmax_unbound_uM\tn_de_0.5\tn_de_1.0\tmedian_abs_lfc_top500\tl2_shrunk")?;
    for c in &cons { let p = &pk[&c.drug]; writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", c.drug, DOSES[c.dose], c.line, p.cmax, p.cmax_u, c.de(0.5).len(), c.de(1.0).len(), mag1(c), c.l2())?; }
    drop(w);
    let mut drugs: Vec<&String> = pk.keys().filter(|d| cons.iter().any(|c| &c.drug == *d)).collect(); drugs.sort();
    let mut brng = Rng(20261116);
    for cut in [0.5, 1.0] {
        let floor: HashMap<&str, f64> = lines.iter().map(|l| (l.as_str(), quantile(&noise.iter().filter(|n| n.line == *l).map(|n| n.de(cut).len() as f64).collect::<Vec<_>>(), 0.95))).collect();
        let get = |d: &str, k: usize, l: &str| cons.iter().find(|c| c.drug == d && c.dose == k && c.line == l);
        for b in 0..3 {
            let frac = |f: &dyn Fn(&Pk) -> bool| drugs.iter().filter(|d| f(&pk[d.as_str()])).count() as f64 / drugs.len() as f64;
            let gt = |dose: f64, mult: f64| frac(&|p: &Pk| { let (t, s) = p.thr_scale(b); dose * s > mult * t });
            // per-drug totals
            let tot: Vec<(f64, f64, f64, f64)> = drugs.iter().map(|d| { let p = &pk[d.as_str()]; let (mut un, mut tn, mut ue, mut te) = (0.0, 0.0, 0.0, 0.0);
                for c in cons.iter().filter(|c| &c.drug == *d) { let n = c.de(cut).len() as f64; let e = (n - floor[c.line.as_str()]).max(0.0); tn += n; te += e; if !p.reach(b, c.dose) { un += n; ue += e } } (un, tn, ue, te) }).collect();
            let share = |ex: bool, rng: &mut Rng| -> (f64, f64, f64) { let f = |v: &[(f64, f64, f64, f64)]| if ex { v.iter().map(|x| x.2).sum::<f64>() / v.iter().map(|x| x.3).sum::<f64>() } else { v.iter().map(|x| x.0).sum::<f64>() / v.iter().map(|x| x.1).sum::<f64>() };
                let bs: Vec<f64> = (0..2000).map(|_| { let s: Vec<_> = (0..tot.len()).map(|_| tot[(rng.next() % tot.len() as u64) as usize]).collect(); f(&s) }).collect(); (f(&tot), quantile(&bs, 0.025), quantile(&bs, 0.975)) };
            let (qr, qe) = (share(false, &mut brng), share(true, &mut brng));
            let cs: Vec<&Con> = cons.iter().collect(); let un: Vec<&&Con> = cs.iter().filter(|c| !pk[&c.drug].reach(b, c.dose)).collect();
            let cond_share = un.len() as f64 / cs.len() as f64;
            let above = |v: &Vec<&&Con>| v.iter().filter(|c| c.de(cut).len() as f64 > floor[c.line.as_str()]).count() as f64 / v.len().max(1) as f64;
            let re: Vec<&&Con> = cs.iter().filter(|c| pk[&c.drug].reach(b, c.dose)).collect();
            let mean = |v: &Vec<&&Con>| v.iter().map(|c| c.de(cut).len() as f64).sum::<f64>() / v.len().max(1) as f64;
            let (mut rat, mut rec, mut rlow) = (vec![], vec![], vec![]); let mut nstr = [0usize; 3];
            for d in &drugs { let p = &pk[d.as_str()]; let st = if !p.reach(b, 0) { 2 } else if !p.reach(b, 3) { 1 } else { 0 }; nstr[st] += 1; if st != 1 { continue }
                let top = (0..3).filter(|&k| p.reach(b, k)).max().unwrap(); let (mut r, mut rc, mut rl) = (vec![], vec![], vec![]);
                for l in &lines { let (Some(hi), Some(lo)) = (get(d, 3, l), get(d, top, l)) else { continue }; let s10: HashSet<usize> = hi.de(cut).into_iter().collect();
                    let e10 = (s10.len() as f64 - floor[l.as_str()]).max(0.0); if e10 <= 0.0 { continue }
                    r.push((lo.de(cut).len() as f64 - floor[l.as_str()]).max(0.0) / e10);
                    rc.push(lo.de(cut).iter().filter(|g| s10.contains(g)).count() as f64 / s10.len() as f64);
                    let low: HashSet<usize> = (0..2).filter_map(|k| get(d, k, l)).flat_map(|c| c.de(cut)).collect(); rl.push(low.iter().filter(|g| s10.contains(g)).count() as f64 / s10.len() as f64) }
                if !r.is_empty() { rat.push(med(&r)); rec.push(med(&rc)); rlow.push(med(&rl)) } }
            metrics.push(json!({"cut": cut, "basis": BASIS[b], "n_compounds": drugs.len(), "n_straddling": nstr[1], "n_always_above": nstr[2], "n_all_reachable": nstr[0],
                "share_1uM_above_cmax": gt(1.0, 1.0), "share_10uM_above_cmax": gt(10.0, 1.0), "share_1uM_above_10x": gt(1.0, 10.0), "share_10uM_above_10x": gt(10.0, 10.0), "share_10nM_above_cmax": gt(0.01, 1.0),
                "share_contrasts_unreachable": cond_share, "share_above_floor_unreachable": above(&un), "share_above_floor_reachable": above(&re), "mean_n_de_unreachable": mean(&un), "mean_n_de_reachable": mean(&re),
                "de_share_above": jnum(qr.0), "de_share_above_ci": [jnum(qr.1), jnum(qr.2)], "excess_share_above": jnum(qe.0), "excess_share_above_ci": [jnum(qe.1), jnum(qe.2)],
                "recovery_n_compounds": rat.len(), "recovery_median_ratio": jnum(med(&rat)), "recall_10uM_at_top_reachable": jnum(med(&rec)), "recall_10uM_at_10_100nM": jnum(med(&rlow)),
                "vehicle_floor": floor.iter().map(|(k, v)| (k.to_string(), jnum(*v))).collect::<Map<_, _>>()}));
        }
    }
    res.insert("metrics".into(), Value::Array(metrics));
    // n_de by dose
    res.insert("n_de_by_dose".into(), (0..4).map(|k| { let v: Vec<f64> = cons.iter().filter(|c| c.dose == k).map(|c| c.de(1.0).len() as f64).collect();
        (format!("{}", DOSES[k]), json!({"mean": v.iter().sum::<f64>() / v.len() as f64, "median": jnum(med(&v)), "q3": jnum(quantile(&v, 0.75))})) }).collect::<Map<_, _>>().into());
    // threshold-free magnitude
    { let m: Vec<(usize, f64, f64, f64, f64)> = cons.iter().map(|c| (c.dose, (DOSES[c.dose] / pk[&c.drug].cmax_u).log10(), mag1(c), c.l2(), c.de(1.0).len() as f64)).collect();
      let vm: Vec<(f64, f64)> = noise.iter().map(|c| (mag1(c), c.l2())).collect(); let mut mj = Map::new();
      mj.insert("vehicle_split".into(), json!({"median_abs_lfc_top500": jnum(med(&vm.iter().map(|x| x.0).collect::<Vec<_>>())), "l2_shrunk": jnum(med(&vm.iter().map(|x| x.1).collect::<Vec<_>>()))}));
      mj.insert("by_dose".into(), (0..4).map(|k| (format!("{}", DOSES[k]), json!({"median_abs_lfc_top500": jnum(med(&m.iter().filter(|x| x.0 == k).map(|x| x.2).collect::<Vec<_>>())), "l2_shrunk": jnum(med(&m.iter().filter(|x| x.0 == k).map(|x| x.3).collect::<Vec<_>>()))}))).collect::<Map<_, _>>().into());
      let bins = [("<=0.1x", f64::NEG_INFINITY, -1.0), ("0.1-1x", -1.0, 0.0), ("1-10x", 0.0, 1.0), ("10-100x", 1.0, 2.0), (">100x", 2.0, f64::INFINITY)];
      mj.insert("by_exposure".into(), bins.iter().map(|(n, lo, hi)| { let s: Vec<&(usize, f64, f64, f64, f64)> = m.iter().filter(|x| x.1 > *lo && x.1 <= *hi).collect();
          (n.to_string(), json!({"n": s.len(), "median_n_de": jnum(med(&s.iter().map(|x| x.4).collect::<Vec<_>>())), "share_n_de_ge_10": s.iter().filter(|x| x.4 >= 10.0).count() as f64 / s.len().max(1) as f64, "median_abs_lfc_top500": jnum(med(&s.iter().map(|x| x.2).collect::<Vec<_>>())), "l2_shrunk": jnum(med(&s.iter().map(|x| x.3).collect::<Vec<_>>())),
              "l2_q1": jnum(quantile(&s.iter().map(|x| x.3).collect::<Vec<_>>(), 0.25)), "l2_q3": jnum(quantile(&s.iter().map(|x| x.3).collect::<Vec<_>>(), 0.75))})) }).collect::<Map<_, _>>().into());
      // within straddling drugs: lowest unreachable over highest reachable dose
      let (mut r1, mut r2, mut h1, mut h2, mut lo2, mut up2) = (vec![], vec![], vec![], vec![], vec![], vec![]);
      for d in &drugs { let p = &pk[d.as_str()]; if !p.reach(0, 0) || p.reach(0, 3) { continue } let top = (0..3).filter(|&k| p.reach(0, k)).max().unwrap();
          for l in &lines { let f = |k: usize| cons.iter().find(|c| &c.drug == *d && c.dose == k && c.line == *l); let (Some(a), Some(b)) = (f(top), f(top + 1)) else { continue };
              r1.push(mag1(b) / mag1(a)); r2.push(b.l2() / a.l2()); h1.push(mag1(a)); h2.push(a.l2()); lo2.push(a.l2()); up2.push(b.l2()) } }
      let wt = wilcoxon(&up2, &lo2);
      mj.insert("straddling".into(), json!({"n_drug_line": r1.len(), "ratio_median_abs_lfc": jnum(med(&r1)), "ratio_l2": jnum(med(&r2)), "ratio_l2_q1": jnum(quantile(&r2, 0.25)), "ratio_l2_q3": jnum(quantile(&r2, 0.75)),
          "highest_reachable_median_abs_lfc": jnum(med(&h1)), "highest_reachable_l2": jnum(med(&h2)), "wilcoxon_p_two_sided": jnum(wt.0)}));
      // above 100x against at or below Cmax, all contrasts
      let (hi, lw): (Vec<f64>, Vec<f64>) = (m.iter().filter(|x| x.1 > 2.0).map(|x| x.3).collect(), m.iter().filter(|x| x.1 <= 0.0).map(|x| x.3).collect());
      let mw = mann_whitney(&hi, &lw); mj.insert("above_100x_vs_reachable".into(), json!({"ratio_of_medians_l2": jnum(med(&hi) / med(&lw)), "mann_whitney_p": jnum(mw.1), "auroc": jnum(mw.0)}));
      res.insert("magnitude".into(), Value::Object(mj)); }
    // downsampling
    { let mut keys: HashSet<(String, String)> = HashSet::new(); for m in meta.iter().filter(|m| m[6] != m[7]) { keys.insert((m[1].clone(), m[0].clone())); }
      let f = |v: &Vec<Con>, d: &str, l: &str| v.iter().find(|c| c.drug == d && c.dose == 3 && c.line == l).map(|c| c.de(1.0).len() as f64);
      let rr: Vec<f64> = keys.iter().filter_map(|(d, l)| match (f(&cons, d, l), f(&cons_ds, d, l)) { (Some(a), Some(b)) if a > 0.0 => Some(b / a), _ => None }).collect();
      let (mut un, mut tn) = (0.0, 0.0);
      for l in &lines { let src = if ds_lines.contains(l.as_str()) { &cons_ds } else { &cons };
          for c in src.iter().filter(|c| c.line == *l) { let n = c.de(1.0).len() as f64; tn += n; if !pk[&c.drug].reach(0, c.dose) { un += n } } }
      res.insert("downsampling".into(), json!({"n_drug_line_pairs": keys.len(), "n_with_de": rr.len(), "median_retained": jnum(med(&rr)), "de_share_above_unbound": jnum(un / tn)})); }
    res.insert("power".into(), power(&cons, &pk, &lines, od)?);
    std::fs::write(rj, serde_json::to_string_pretty(&Value::Object(res))?)?;
    eprintln!("done [{:.0}s]", t0.elapsed().as_secs_f64());
    Ok(())
}

// ------------------------------------------------------------------ power at clinical exposure
/// One compound and cell line with at least 10 DE genes at 10 uM, and one lower dose of it.
struct Unit { drug: usize, band: usize, n: f64, sb: f64, sb2: f64, ss2: f64, ctrl: f64, top: f64, ztop: f64, inst: Vec<(f64, [f64; 5], f64)> } // inst: per-cell rate, line parameters, z threshold
const LN2: f64 = std::f64::consts::LN_2;
fn phi(x: f64) -> f64 { 1.0 - norm_sf(x) }
fn zq(p: f64) -> f64 { use statrs::distribution::{ContinuousCDF, Normal}; Normal::new(0.0, 1.0).unwrap().inverse_cdf(p) }
/// Variance of the log count of a well of n cells: Poisson and count-dependent dispersion scale with 1/(n m), a0 does not.
fn vwell(n: f64, m: f64, p: &[f64; 5]) -> f64 { if n.is_infinite() { p[0] } else { (1.0 + p[1]) / (n * m) + p[0] } }
/// Standard error of the log2 fold change for r treated wells of n cells against the vehicle wells of the line.
fn se_model(n: f64, r: f64, m: f64, p: &[f64; 5]) -> f64 { (vwell(n, m, p) / r + vwell(p[4], m, p) / p[3].max(r)).sqrt() / LN2 }
fn med_power(inst: &[&(f64, [f64; 5], f64)], beta: f64, n: f64, r: f64) -> f64 { med(&inst.iter().map(|(m, p, z)| phi(beta / se_model(n, r, *m, p) - z)).collect::<Vec<_>>()) }
fn med_se(inst: &[&(f64, [f64; 5], f64)], n: f64, r: f64) -> f64 { med(&inst.iter().map(|(m, p, _)| se_model(n, r, *m, p)).collect::<Vec<_>>()) }
/// Smallest x in [lo, hi] (log scale) with f(x) true, for monotone f; None if f(hi) is false.
fn solve<F: Fn(f64) -> bool>(f: F, lo: f64, hi: f64) -> Option<f64> { if !f(hi) { return None } if f(lo) { return Some(lo) } let (mut a, mut b) = (lo.ln(), hi.ln()); for _ in 0..32 { let c = 0.5 * (a + b); if f(c.exp()) { b = c } else { a = c } } Some(b.exp()) }
struct Req { effect: f64, tau: f64, b90: f64, n90: Option<f64>, r90: Option<f64>, ceil90: f64, nmean: Option<f64>, rmean: Option<f64>, ceilmean: f64, nz: Option<f64>, rz: Option<f64>, ztop: f64, top: f64 }
fn require(us: &[&Unit]) -> Req {
    let n: f64 = us.iter().map(|u| u.n).sum(); let mb = us.iter().map(|u| u.sb).sum::<f64>() / n; let effect = mb - us.iter().map(|u| u.ctrl * u.n).sum::<f64>() / n;
    let tau = (us.iter().map(|u| u.sb2).sum::<f64>() / n - mb * mb - us.iter().map(|u| u.ss2).sum::<f64>() / n).max(0.0).sqrt(); let b90 = effect + 1.2816 * tau;
    let inst: Vec<&(f64, [f64; 5], f64)> = us.iter().flat_map(|u| u.inst.iter()).collect();
    let (ztop, top) = (us.iter().map(|u| u.ztop * u.n).sum::<f64>() / n, us.iter().map(|u| u.top * u.n).sum::<f64>() / n);
    let need = |beta: f64| -> (Option<f64>, Option<f64>, f64) { if !(beta > 0.0) { return (None, None, 0.0) }
        (solve(|x| med_power(&inst, beta, x, 2.0) >= 0.8, 1.0, 1e9), solve(|x| med_power(&inst, beta, 150.0, x.ceil()) >= 0.8, 1.0, 1e7).map(|x| x.ceil()), med_power(&inst, beta, f64::INFINITY, 2.0)) };
    let (n90, r90, ceil90) = need(b90); let (nmean, rmean, ceilmean) = need(effect);
    let (nz, rz) = if effect > 0.0 { let t = effect / ztop; (solve(|x| med_se(&inst, x, 2.0) <= t, 1.0, 1e9), solve(|x| med_se(&inst, 150.0, x.ceil()) <= t, 1.0, 1e7).map(|x| x.ceil())) } else { (None, None) };
    Req { effect, tau, b90, n90, r90, ceil90, nmean, rmean, ceilmean, nz, rz, ztop, top }
}
fn power(cons: &[Con], pk: &HashMap<String, Pk>, lines: &[String], od: &str) -> Result<Value> {
    let mut drugs: Vec<&str> = cons.iter().map(|c| c.drug.as_str()).collect::<HashSet<_>>().into_iter().collect(); drugs.sort();
    let band = |c: &Con| { let x = DOSES[c.dose] / pk[&c.drug].cmax_u; if x <= 1.0 { 0 } else if x <= 10.0 { 1 } else if x <= 100.0 { 2 } else { 3 } };
    let mut units: Vec<Unit> = vec![]; let mut nunit = 0;
    for l in lines { let lc: Vec<&Con> = cons.iter().filter(|c| c.line == *l).collect(); let ctrl: Vec<&&Con> = lc.iter().filter(|c| pk[&c.drug].reach(0, c.dose)).collect();
        for (di, d) in drugs.iter().enumerate() { let Some(hi) = lc.iter().find(|c| c.drug == *d && c.dose == 3) else { continue };
            let r: Vec<usize> = (0..hi.lfc.len()).filter(|&i| hi.padj[i] < 0.05 && hi.lfc[i].abs() >= 1.0).collect(); if r.len() < 10 { continue } nunit += 1;
            let sg: Vec<f64> = r.iter().map(|&i| hi.lfc[i].signum()).collect(); let n = r.len() as f64;
            let proj = |c: &Con| r.iter().zip(&sg).map(|(&i, s)| c.lfc[i] * s).sum::<f64>() / n;
            let oc: Vec<f64> = ctrl.iter().filter(|c| c.drug != *d).map(|c| proj(c)).collect(); let cm = oc.iter().sum::<f64>() / oc.len().max(1) as f64;
            let zthr = zq(1.0 - 0.5 * 0.05 * 0.8 * n / hi.lfc.len() as f64);
            let inst: Vec<(f64, [f64; 5], f64)> = r.iter().map(|&i| (hi.base[i] / hi.par[2], hi.par, zthr)).collect();
            let (top, ztop) = (r.iter().map(|&i| hi.lfc[i].abs()).sum::<f64>() / n, r.iter().map(|&i| (hi.lfc[i] / hi.se[i]).abs()).sum::<f64>() / n);
            for c in lc.iter().filter(|c| c.drug == *d && c.dose < 3) { let b: Vec<f64> = r.iter().zip(&sg).map(|(&i, s)| c.lfc[i] * s).collect();
                units.push(Unit { drug: di, band: band(c), n, sb: b.iter().sum(), sb2: b.iter().map(|x| x * x).sum(), ss2: r.iter().map(|&i| c.se[i] * c.se[i]).sum(), ctrl: cm, top, ztop, inst: inst.clone() }) } } }
    let names = ["<=1x", "1-10x", "10-100x", ">100x"]; let mut rng = Rng(11); let mut out = Map::new();
    out.insert("n_units".into(), json!(nunit)); out.insert("n_compounds".into(), json!(units.iter().map(|u| u.drug).collect::<HashSet<_>>().len()));
    let opt = |x: Option<f64>| x.map(jnum).unwrap_or(Value::Null);
    let ci = |v: &[f64]| json!([jnum(quantile(v, 0.025)), jnum(quantile(v, 0.975))]);
    let cio = |v: &[Option<f64>]| { let f: Vec<f64> = v.iter().map(|x| x.unwrap_or(f64::INFINITY)).collect(); let q = |p: f64| { let mut s = f.clone(); s.sort_by(|a, b| a.partial_cmp(b).unwrap()); s[((s.len() - 1) as f64 * p).round() as usize] }; json!([jnum(q(0.025)), jnum(q(0.975))]) };
    for (bi, bn) in names.iter().enumerate() {
        let us: Vec<&Unit> = units.iter().filter(|u| u.band == bi).collect(); if us.is_empty() { continue }
        let ds: Vec<usize> = us.iter().map(|u| u.drug).collect::<HashSet<_>>().into_iter().collect::<std::collections::BTreeSet<_>>().into_iter().collect();
        let q = require(&us);
        let seeds: Vec<Vec<usize>> = (0..1000).map(|_| (0..ds.len()).map(|_| ds[(rng.next() % ds.len() as u64) as usize]).collect()).collect();
        let bs: Vec<Req> = seeds.par_iter().map(|s| { let v: Vec<&Unit> = s.iter().flat_map(|d| us.iter().filter(move |u| u.drug == *d).cloned()).collect(); require(&v) }).collect();
        let col = |f: &dyn Fn(&Req) -> f64| ci(&bs.iter().map(|r| f(r)).collect::<Vec<_>>()); let colo = |f: &dyn Fn(&Req) -> Option<f64>| cio(&bs.iter().map(|r| f(r)).collect::<Vec<_>>());
        let inst: Vec<&(f64, [f64; 5], f64)> = us.iter().flat_map(|u| u.inst.iter()).collect();
        let mde = |n: f64, r: f64| med(&inst.iter().map(|(m, p, z)| (z + zq(0.8)) * se_model(n, r, *m, p)).collect::<Vec<_>>());
        out.insert(bn.to_string(), json!({"n_units": us.len(), "n_compounds": ds.len(), "n_gene_instances": inst.len(),
            "effect_log2": jnum(q.effect), "effect_ci": col(&|r| r.effect), "tau": jnum(q.tau), "effect_p90": jnum(q.b90), "effect_p90_ci": col(&|r| r.b90),
            "effect_10uM": jnum(q.top), "retained_fraction": jnum(q.effect / q.top), "retained_ci": col(&|r| r.effect / r.top), "mean_abs_z_10uM": jnum(q.ztop),
            "control_projection": jnum(us.iter().map(|u| u.ctrl * u.n).sum::<f64>() / us.iter().map(|u| u.n).sum::<f64>()),
            "median_cell_rate": jnum(med(&inst.iter().map(|x| x.0).collect::<Vec<_>>())), "median_alpha_star": jnum(med(&inst.iter().map(|x| 2.0 * norm_sf(x.2)).collect::<Vec<_>>())),
            "mde_150_cells_2_wells": jnum(mde(150.0, 2.0)), "mde_infinite_cells_2_wells": jnum(mde(f64::INFINITY, 2.0)), "mde_150_cells_8_wells": jnum(mde(150.0, 8.0)),
            "p90": {"cells_per_well_2_wells": opt(q.n90), "cells_ci": colo(&|r| r.n90), "wells_150_cells": opt(q.r90), "wells_ci": colo(&|r| r.r90), "power_infinite_cells_2_wells": jnum(q.ceil90), "ceiling": q.ceil90 < 0.8},
            "mean": {"cells_per_well_2_wells": opt(q.nmean), "cells_ci": colo(&|r| r.nmean), "wells_150_cells": opt(q.rmean), "wells_ci": colo(&|r| r.rmean), "power_infinite_cells_2_wells": jnum(q.ceilmean), "ceiling": q.ceilmean < 0.8},
            "match_z_10uM": {"cells_per_well_2_wells": opt(q.nz), "wells_150_cells": opt(q.rz), "wells_ci": colo(&|r| r.rz), "cells_ci": colo(&|r| r.nz)}}));
        if bi == 0 { let mut w = BufWriter::new(std::fs::File::create(format!("{od}/power_curve.tsv"))?); writeln!(w, "effect\twells\tcells_per_well\tmedian_power")?;
            for (en, e) in [("p90", q.b90), ("mean", q.effect)] { for r in [2.0, 4.0, 8.0, 16.0, 32.0] { for k in 0..=40 { let n = 10f64.powf(1.0 + 0.1 * k as f64);
                writeln!(w, "{}\t{}\t{:.1}\t{:.5}", en, r, n, if e > 0.0 { med_power(&inst, e, n, r) } else { 0.0 })?; } } } }
    }
    // trend parameters, depth and agreement of the variance model with the fitted standard errors
    let mut lj = Map::new();
    for l in lines { let c = cons.iter().find(|c| c.line == *l).unwrap(); let m: Vec<f64> = c.base.iter().map(|b| b / c.par[2]).collect();
        let ratio: Vec<f64> = cons.iter().filter(|c| c.line == *l).map(|c| med(&(0..c.se.len()).map(|i| se_model(c.ncell, 2.0, c.base[i] / c.par[2], &c.par) / c.se[i]).collect::<Vec<_>>())).collect();
        lj.insert(l.clone(), json!({"a0": c.par[0], "a1": c.par[1], "median_cells_per_well": c.par[2], "vehicle_wells": c.par[3], "median_cells_per_vehicle_well": c.par[4], "counts_per_cell_tested_genes": m.iter().sum::<f64>(), "median_gene_rate_per_cell": jnum(med(&m)),
            "well_cv_floor": c.par[0].sqrt(), "model_over_fitted_se_median": jnum(med(&ratio))})); }
    out.insert("lines".into(), Value::Object(lj));
    Ok(Value::Object(out))
}
fn summ(v: &[f64]) -> Value { json!({"n": v.len(), "median": jnum(med(v)), "q1": jnum(quantile(v, 0.25)), "q3": jnum(quantile(v, 0.75))}) }

fn main() -> Result<()> {
    let a: Vec<String> = std::env::args().collect();
    match a.get(1).map(|s| s.as_str()) {
        Some("pseudobulk") => pseudobulk(&a[2], &a[3], &a[4]),
        Some("de") => de(&a[2], &a[3], &a[4], &a[5]),
        _ => Err(anyhow!("usage: sciplex_de_power pseudobulk <h5ad> <data dir> <work dir> | sciplex_de_power de <work dir> <data dir> <out dir> <results json>")),
    }
}
