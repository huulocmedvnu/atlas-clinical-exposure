//! Merges the result files of the three analyses and writes the values of Tables 1 to 3 as CSV.
//! Usage: export_tables <results dir>
//! Reads <results dir>/metrics/{tahoe,sciplex,patients}.json. Writes final_metrics.json and table1 to table3 CSV files.
//!
//! Author: Thai Huu Loc Tran
use anyhow::{anyhow, Result};
use serde_json::{json, Value};
use std::io::Write;

fn pct(v: &Value) -> String { format!("{:.1}%", 100.0 * v.as_f64().unwrap_or(f64::NAN)) }
fn pct_ci(v: &Value, ci: &Value) -> String { format!("{:.1}% ({:.1} to {:.1})", 100.0 * v.as_f64().unwrap(), 100.0 * ci[0].as_f64().unwrap(), 100.0 * ci[1].as_f64().unwrap()) }
fn frac(a: &Value, n: &Value) -> String { let (a, n) = (a.as_f64().unwrap(), n.as_f64().unwrap()); format!("{}/{} ({:.1}%)", a, n, 100.0 * a / n) }
fn q(s: &str) -> String { if s.contains(',') { format!("\"{s}\"") } else { s.to_string() } }

fn main() -> Result<()> {
    let dir = std::env::args().nth(1).ok_or(anyhow!("usage: export_tables <results dir>"))?;
    let rd = |n: &str| -> Result<Value> { Ok(serde_json::from_str(&std::fs::read_to_string(format!("{dir}/metrics/{n}.json"))?)?) };
    let (t, s, b) = (rd("tahoe")?, rd("sciplex")?, rd("patients")?);
    std::fs::write(format!("{dir}/final_metrics.json"), serde_json::to_string_pretty(&json!({"tahoe": t, "sciplex": s, "biopsy": b}))?)?;

    // Table 1
    let tr = |set: &str, basis: &str, stratum: &str| t["sweep"].as_array().unwrap().iter().find(|x| x["cut"] == "1.0" && x["pk_set"] == set && x["basis"] == basis && x["stratum"] == stratum).cloned().unwrap_or(Value::Null);
    let sr = |basis: &str| s["metrics"].as_array().unwrap().iter().find(|x| x["cut"] == 1.0 && x["basis"] == basis).cloned().unwrap();
    let mut w = std::fs::File::create(format!("{dir}/table1_exposure_summary.csv"))?;
    writeln!(w, "quantity,tahoe_100m_compilation,tahoe_100m_labels,sciplex_3")?;
    let mut row = |name: &str, a: String, b: String, c: String| writeln!(w, "{},{},{},{}", q(name), q(&a), q(&b), q(&c));
    row("Cell lines", "50".into(), "50".into(), s["data"]["cell_lines"].to_string())?;
    row("Nominal doses (uM)", "0.05, 0.5, 5".into(), "0.05, 0.5, 5".into(), "0.01, 0.1, 1, 10".into())?;
    let n = |v: Value| v["n_compounds"].to_string();
    row("Annotated compounds", n(tr("liston", "unbound", "all")), n(tr("label", "unbound", "all")), sr("unbound")["n_compounds"].to_string())?;
    for (name, st, key) in [("Doses on both sides of unbound Cmax", "straddling", "n_straddling"), ("All doses above unbound Cmax", "always_above", "n_always_above"), ("All doses reachable", "all_reachable", "n_all_reachable")] {
        row(name, n(tr("liston", "unbound", st)), n(tr("label", "unbound", st)), sr("unbound")[key].to_string())?; }
    for (basis, name) in [("unbound", "unbound Cmax"), ("fbs", "serum-adjusted reference"), ("total", "total Cmax")] {
        let (a, l, c) = (tr("liston", basis, "all"), tr("label", basis, "all"), sr(basis));
        row(&format!("Conditions above {name}"), pct(&a["share_contrasts_unreachable"]), pct(&l["share_contrasts_unreachable"]), pct(&c["share_contrasts_unreachable"]))?;
        row(&format!("DE calls from conditions above {name} (95% CI)"), pct_ci(&a["de_share_above"], &a["de_share_above_ci"]), pct_ci(&l["de_share_above"], &l["de_share_above_ci"]), pct_ci(&c["de_share_above"], &c["de_share_above_ci"]))?; }

    // Tables 2 and 3
    let order = [("everolimus", "Everolimus", "GSE119262"), ("afatinib", "Afatinib", "GSE192464"), ("palbociclib", "Palbociclib", "GSE93204"), ("anastrozole", "Anastrozole", "GSE153470"), ("sunitinib", "Sunitinib", "GSE58837"),
                 ("anastrozole_b", "Anastrozole", "GSE93204"), ("olaparib", "Olaparib", "GSE341672")];
    let (mut w2, mut w3) = (std::fs::File::create(format!("{dir}/table2_patient_gene_agreement.csv"))?, std::fs::File::create(format!("{dir}/table3_top50_and_hallmark_sets.csv"))?);
    writeln!(w2, "compound,cohort,set,pairs,lowest_dose_over_unbound_cmax,signature_0.05uM,signature_5uM_only")?;
    writeln!(w3, "compound,cohort,set,top50_0.05uM,top50_5uM_only,gene_sets_0.05uM,gene_sets_0.05uM_rho,gene_sets_5uM,gene_sets_5uM_rho")?;
    let (p, a) = (&b["pooled"], &b["aggregated"]);
    for (i, (k, nm, gse)) in order.iter().enumerate() {
        if i == 5 {
            writeln!(w2, "Pooled,,primary,,,{},{}", frac(&p["low"][0], &p["low"][1]), frac(&p["high_only"][0], &p["high_only"][1]))?;
            let (t5, pw) = (&a["top50_primary"], &a["pathway_primary"]);
            writeln!(w3, "Pooled,,primary,{},{},{},,{},", frac(&t5["low"][0], &t5["low"][1]), frac(&t5["high"][0], &t5["high"][1]), frac(&pw["low"][0], &pw["low"][1]), frac(&pw["high"][0], &pw["high"][1]))?;
        }
        let c = &b["cohorts"][*k]; let set = if c["primary"] == true { "primary" } else { "secondary" };
        let r = c["lowest_dose_over_unbound_cmax"].as_f64().unwrap(); let rs = if r >= 10.0 { format!("{r:.0}") } else if r >= 1.0 { format!("{r:.1}") } else { format!("{r:.2}") };
        writeln!(w2, "{nm},{gse},{set},{},{rs},{},{}", c["n_pairs"], frac(&c["low"]["agree"], &c["low"]["n_measured"]), frac(&c["high_only"]["agree"], &c["high_only"]["n_measured"]))?;
        let (t5, pw) = (&c["top"]["50"], &c["pathway"]); let f2 = |v: &Value| format!("{}/{}", v["agree"], v.get("n").or(v.get("n_sets")).unwrap());
        writeln!(w3, "{nm},{gse},{set},{},{},{},{:.2},{},{:.2}", f2(&t5["low"]), f2(&t5["high_only"]), f2(&pw["0.05uM"]), pw["0.05uM"]["rho"].as_f64().unwrap(), f2(&pw["5uM"]), pw["5uM"]["rho"].as_f64().unwrap())?;
    }
    Ok(())
}
