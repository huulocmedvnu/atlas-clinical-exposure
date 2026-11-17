"""Figures 1 to 7 from the result tables in results/, with the source data of every panel.

Usage: python3 scripts/plot_figures.py
Writes figures/fig*.pdf and .png and results/figure_source_data/*.csv. Figure 1 additionally needs rsvg-convert.
Style: plotly_white, Plotly colours, Arial, black text, no gridlines.

Author: Thai Huu Loc Tran
"""
import json, os, subprocess, numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
R, F, SD = os.path.join(ROOT, "results/"), os.path.join(ROOT, "figures/"), os.path.join(ROOT, "results/figure_source_data/")
os.makedirs(F, exist_ok=True); os.makedirs(SD, exist_ok=True)
J = json.load(open(R + "final_metrics.json")); T, S, B = J["tahoe"], J["sciplex"], J["biopsy"]
NAMES = {"Figure2": "fig2_exposure_ratios", "Figure3": "fig3_response_by_exposure", "Figure4": "fig4_specificity_moa", "Figure5": "fig5_patient_gene_agreement", "Figure6": "fig6_power_and_hallmarks", "Figure7": "fig7_response_model"}
def source(panel, df): pd.DataFrame(df).to_csv(SD + panel + ".csv", index=False)
C = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#19D3F3", "#FF6692"]; K = "#111111"
def style(fig, w, h):
    fig.update_layout(template="plotly_white", width=w, height=h, font=dict(family="Arial", size=13, color=K), margin=dict(l=70, r=20, t=40, b=60), legend=dict(font=dict(color=K, size=12)))
    fig.update_xaxes(showgrid=False, showline=True, linecolor=K, ticks="outside", tickfont=dict(color=K), title_font=dict(color=K), zeroline=False)
    fig.update_yaxes(showgrid=False, showline=True, linecolor=K, ticks="outside", tickfont=dict(color=K), title_font=dict(color=K), zeroline=False)
    for a in fig.layout.annotations:
        if a.text[:2] in ("a ", "b ", "c ", "d "): a.font = dict(family="Arial", size=14, color=K)
    return fig
def save(fig, name): fig.write_image(F + NAMES[name] + ".png", scale=4); fig.write_image(F + NAMES[name] + ".pdf")
BINS = ["<=0.1x", "0.1-1x", "1-10x", "10-100x", ">100x"]; BL = ["≤0.1", "0.1 to 1", "1 to 10", "10 to 100", ">100"]

# ------------------------------------------------ Figure 2: doses against unbound Cmax; share of conditions and of DE calls
st = pd.read_csv(R + "supporting/strata.tsv", sep="\t"); u = st[(st.basis == "unbound") & (st.pk_source == "liston")].sort_values("cmax_unbound_uM").reset_index(drop=True)
sc = pd.read_csv(R + "supporting/sciplex_contrasts.tsv", sep="\t"); su = sc.drop_duplicates("drug").sort_values("cmax_unbound_uM").reset_index(drop=True)
fig = make_subplots(rows=1, cols=3, column_widths=[0.36, 0.3, 0.34], subplot_titles=["a  Tahoe-100M, 65 compounds", "b  sci-Plex 3, 27 compounds", "c  Above unbound C<sub>max</sub>"], horizontal_spacing=0.09)
for k, (d, c) in enumerate(zip([0.05, 0.5, 5], C)):
    fig.add_trace(go.Scatter(x=np.log10(d / u.cmax_unbound_uM), y=u.index, mode="markers", marker=dict(color=c, size=5), name=f"{d} µM", legend="legend"), 1, 1)
for k, (d, c) in enumerate(zip([0.01, 0.1, 1, 10], [C[5], C[0], C[3], C[1]])):
    fig.add_trace(go.Scatter(x=np.log10(d / su.cmax_unbound_uM), y=su.index, mode="markers", marker=dict(color=c, size=6), name=f"{d} µM", legend="legend2"), 1, 2)
for col in (1, 2):
    fig.add_vline(x=0, line=dict(color="#888888", dash="dash", width=1), row=1, col=col)
    fig.update_xaxes(title_text="Dose / unbound C<sub>max</sub>", tickvals=[-4, -2, 0, 2, 4], ticktext=["10<sup>−4</sup>", "10<sup>−2</sup>", "1", "10<sup>2</sup>", "10<sup>4</sup>"], range=[-4.6, 4.6], row=1, col=col)
    fig.update_yaxes(title_text="Compounds, ordered by unbound C<sub>max</sub>" if col == 1 else None, showticklabels=False, ticks="", row=1, col=col)
row = lambda b, key: next(x for x in T["sweep"] if x["cut"] == "1.0" and x["basis"] == b and x["pk_set"] == "liston" and x["stratum"] == "all")
srow = lambda b: next(x for x in S["metrics"] if x["cut"] == 1.0 and x["basis"] == b)
labs = ["Tahoe<br>unbound", "Tahoe<br>serum-adj.", "Tahoe<br>total", "sci-Plex<br>unbound", "sci-Plex<br>serum-adj.", "sci-Plex<br>total"]
rows = [row("unbound", 0), row("fbs", 0), row("total", 0), srow("unbound"), srow("fbs"), srow("total")]
fig.add_trace(go.Bar(x=labs, y=[100 * r["share_contrasts_unreachable"] for r in rows], marker_color=C[2], name="Conditions", legend="legend3"), 1, 3)
fig.add_trace(go.Bar(x=labs, y=[100 * r["de_share_above"] for r in rows], marker_color=C[3], name="DE calls", legend="legend3",
                     error_y=dict(type="data", symmetric=False, array=[100 * (r["de_share_above_ci"][1] - r["de_share_above"]) for r in rows], arrayminus=[100 * (r["de_share_above"] - r["de_share_above_ci"][0]) for r in rows], color=K, thickness=1, width=3)), 1, 3)
fig.update_yaxes(title_text="Share (%)", range=[0, 118], tickvals=[0, 20, 40, 60, 80, 100], row=1, col=3); fig.update_xaxes(tickfont=dict(size=10), row=1, col=3)
style(fig, 1250, 520).update_layout(barmode="group", legend=dict(x=0.24, y=0.99, bgcolor="rgba(255,255,255,0.7)"), legend2=dict(x=0.59, y=0.99, font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"),
                                    legend3=dict(x=0.74, y=1.0, orientation="h", font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"))
save(fig, "Figure2")
source("fig2a_tahoe_dose_ratios", pd.DataFrame({"compound": u.drug, "cmax_unbound_uM": u.cmax_unbound_uM, **{f"ratio_{d}uM": d / u.cmax_unbound_uM for d in [0.05, 0.5, 5]}}))
source("fig2b_sciplex_dose_ratios", pd.DataFrame({"compound": su.drug.str.strip(), "cmax_unbound_uM": su.cmax_unbound_uM, **{f"ratio_{d}uM": d / su.cmax_unbound_uM for d in [0.01, 0.1, 1, 10]}}))
source("fig2c_shares", pd.DataFrame({"atlas_and_reference": [l.replace("<br>", " ") for l in labs], "conditions_above_pct": [100 * r["share_contrasts_unreachable"] for r in rows], "de_calls_above_pct": [100 * r["de_share_above"] for r in rows],
                                     "de_calls_ci_low_pct": [100 * r["de_share_above_ci"][0] for r in rows], "de_calls_ci_high_pct": [100 * r["de_share_above_ci"][1] for r in rows]}))

# ------------------------------------------------ Figure 3: response against exposure
c = pd.read_csv(R + "supporting/contrasts.tsv", sep="\t"); cu = u.set_index("drug").cmax_unbound_uM
inert = c[c.drug.isin(pd.read_csv(ROOT + "/data/pk/compound_moa_annotations.csv").query("inert_reference == 1").compound)]; t = c[c.pk_source == "liston"].copy(); t["bin"] = pd.cut(np.log10(t.dose / t.drug.map(cu)), [-9, -1, 0, 1, 2, 9], labels=BINS)
sc["bin"] = pd.cut(np.log10(sc.dose_uM / sc.cmax_unbound_uM), [-9, -1, 0, 1, 2, 9], labels=BINS)
fig = make_subplots(rows=2, cols=2, subplot_titles=["a  Tahoe-100M", "b  sci-Plex 3", "c  sci-Plex 3", "d  Tahoe-100M"], horizontal_spacing=0.13, vertical_spacing=0.2)
for b, l in zip(BINS, BL):
    fig.add_trace(go.Box(y=t[t.bin == b]["n_de_1.0"] + 1, name=l, marker_color=C[0], boxpoints=False, showlegend=False), 1, 1)
    fig.add_trace(go.Box(y=sc[sc.bin == b]["n_de_1.0"] + 1, name=l, marker_color=C[1], boxpoints="all", jitter=0.5, pointpos=0, marker=dict(size=3), showlegend=False), 1, 2)
    fig.add_trace(go.Box(y=sc[sc.bin == b]["l2_shrunk"], name=l, marker_color=C[1], boxpoints="all", jitter=0.5, pointpos=0, marker=dict(size=3), showlegend=False), 2, 1)
im = inert["n_de_1.0"].median()  # shapes on a log axis take the data value, annotations its log10
fig.add_shape(type="line", xref="x domain", x0=0, x1=1, yref="y", y0=im + 1, y1=im + 1, line=dict(color=C[2], dash="dash"))
fig.add_annotation(xref="x", yref="y", x=6.0, y=np.log10(im + 1), text=f"inert compounds,<br>median ({im:.0f})", showarrow=False, xanchor="right", yanchor="bottom", align="right", font=dict(color=K, size=11))
fig.update_xaxes(range=[-0.6, 6.05], row=1, col=1)
fig.add_hline(y=S["magnitude"]["vehicle_split"]["l2_shrunk"], line=dict(color=C[2], dash="dash"), row=2, col=1, annotation_text="vehicle against vehicle, median", annotation_position="bottom right", annotation_font=dict(color=K, size=11))
for (r_, c_) in [(1, 1), (1, 2)]: fig.update_yaxes(type="log", dtick=1, title_text="DE genes per contrast + 1", row=r_, col=c_)
fig.update_yaxes(title_text="Response magnitude (shrunken log<sub>2</sub>FC, L2 norm)", row=2, col=1)
for (r_, c_) in [(1, 1), (1, 2), (2, 1)]: fig.update_xaxes(title_text="Dose / unbound C<sub>max</sub>", row=r_, col=c_)
sw = pd.DataFrame(T["sweep"]); sw = sw[(sw.pk_set == "liston") & (sw.stratum == "all")]
for b, nm, col in [("unbound", "Unbound C<sub>max</sub>", C[0]), ("fbs", "Serum-adjusted", C[4]), ("total", "Total C<sub>max</sub>", C[3])]:
    x = sw[sw.basis == b]; ci = np.array(x.de_share_above_ci.tolist())
    fig.add_trace(go.Scatter(x=x.cut.astype(float), y=100 * x.de_share_above, mode="lines+markers", marker=dict(color=col, symbol="diamond", size=8), line=dict(color=col), name=nm,
                             error_y=dict(type="data", symmetric=False, array=100 * (ci[:, 1] - x.de_share_above), arrayminus=100 * (x.de_share_above - ci[:, 0]), thickness=1, width=3)), 2, 2)
    fig.add_trace(go.Scatter(x=x.cut.astype(float), y=100 * x.share_contrasts_unreachable, mode="lines", line=dict(color=col, dash="dot"), showlegend=False), 2, 2)
fig.update_xaxes(title_text="|log<sub>2</sub>FC| threshold", tickvals=[0.25, 0.5, 0.75, 1, 1.5], row=2, col=2); fig.update_yaxes(title_text="DE calls from doses above C<sub>max</sub> (%)", range=[0, 100], row=2, col=2)
style(fig, 1150, 900).update_layout(legend=dict(x=0.6, y=0.02, bgcolor="rgba(255,255,255,0.7)"))
save(fig, "Figure3")
source("fig3a_tahoe_de_by_exposure", t[["drug", "dose", "line", "n_de_1.0", "bin"]].rename(columns={"drug": "compound", "dose": "dose_uM", "line": "cell_line", "n_de_1.0": "n_de_genes", "bin": "dose_over_unbound_cmax"}).assign(inert_median=im))
source("fig3bc_sciplex_by_exposure", sc[["drug", "dose_uM", "line", "n_de_1.0", "l2_shrunk", "bin"]].rename(columns={"drug": "compound", "line": "cell_line", "n_de_1.0": "n_de_genes", "l2_shrunk": "response_magnitude", "bin": "dose_over_unbound_cmax"}).assign(vehicle_magnitude_median=S["magnitude"]["vehicle_split"]["l2_shrunk"]))
source("fig3d_threshold_sweep", pd.DataFrame({"reference": sw.basis, "abs_log2fc_threshold": sw.cut, "de_calls_above_pct": 100 * sw.de_share_above, "ci_low_pct": [100 * c_[0] for c_ in sw.de_share_above_ci], "ci_high_pct": [100 * c_[1] for c_ in sw.de_share_above_ci], "conditions_above_pct": 100 * sw.share_contrasts_unreachable}))

# ------------------------------------------------ Figure 4: specificity of the high-dose response
p = pd.read_csv(R + "supporting/overlap_pairs.tsv", sep="\t"); p = p[p.basis == "unbound"]
fig = make_subplots(rows=1, cols=3, subplot_titles=["a  Overlap between mechanisms", "b  Stress and apoptosis gene sets", "c  Mechanism similarity"], horizontal_spacing=0.1)
for stt, lab in [(1, "Lowest dose reachable"), (2, "All doses above C<sub>max</sub>")]:
    x = p[p.stratum == stt]
    fig.add_trace(go.Box(y=x.j_low, x=[lab] * len(x), name="Low dose", marker_color=C[0], boxpoints=False, legendgroup="l", showlegend=stt == 1, legend="legend"), 1, 1)
    fig.add_trace(go.Box(y=x.j_high_only, x=[lab] * len(x), name="5 µM only", marker_color=C[1], boxpoints=False, legendgroup="h", showlegend=stt == 1, legend="legend"), 1, 1)
fig.update_yaxes(title_text="Jaccard index, compound pairs with different mechanisms", range=[0, 0.12], row=1, col=1); fig.update_xaxes(tickfont=dict(size=11), row=1, col=1)
o = T["overlap"]["unbound_stratum1"]["stress"]; gs = [("HALLMARK_APOPTOSIS", "Apoptosis"), ("HALLMARK_P53_PATHWAY", "p53 pathway"), ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "Unfolded protein<br>response"), ("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "Reactive oxygen<br>species")]
for key, nm, col in [("share_low", "Reachable dose", C[0]), ("share_shared", "5 µM, shared with reachable", C[2]), ("share_high_only", "5 µM only", C[1])]:
    fig.add_trace(go.Bar(x=[g[1] for g in gs], y=[100 * o[g[0]][key] for g in gs], name=nm, marker_color=col, legend="legend2"), 1, 2)
fig.update_yaxes(title_text="Share of DE genes in the gene set (%)", range=[0, 3.1], row=1, col=2); fig.update_xaxes(tickfont=dict(size=10), row=1, col=2)
be = T["mechanism"]["by_exposure"]; names = ["<=1x", "1-10x", "10-100x", ">100x"]; lab = ["≤1", "1 to 10", "10 to 100", ">100"]
fig.add_trace(go.Scatter(x=lab, y=[be[n]["median_line_auroc"] for n in names], mode="markers", marker=dict(color=C[3], symbol="diamond", size=11), showlegend=False,
                         error_y=dict(type="data", symmetric=False, array=[be[n]["line_auroc_q3"] - be[n]["median_line_auroc"] for n in names], arrayminus=[be[n]["median_line_auroc"] - be[n]["line_auroc_q1"] for n in names], color=K, thickness=1, width=4)), 1, 3)
fig.add_hline(y=0.5, line=dict(color="#888888", dash="dash", width=1), row=1, col=3)
fig.update_yaxes(title_text="AUROC, same against different mechanism", range=[0.3, 0.75], row=1, col=3); fig.update_xaxes(title_text="Dose / unbound C<sub>max</sub> (both compounds)", type="category", row=1, col=3)
style(fig, 1300, 520).update_layout(boxmode="group", barmode="group", legend=dict(x=0.115, y=0.99, bgcolor="rgba(255,255,255,0.7)"), legend2=dict(x=0.37, y=0.99, font=dict(color=K, size=11), bgcolor="rgba(255,255,255,0.7)"))
save(fig, "Figure4")
source("fig4a_overlap_between_mechanisms", p[["stratum", "j_low", "j_high_only"]].rename(columns={"stratum": "group_1_lowest_dose_reachable_2_all_above", "j_low": "jaccard_low_dose", "j_high_only": "jaccard_5uM_only"}))
source("fig4b_stress_gene_sets", pd.DataFrame([{"gene_set": g[0], "reachable_dose_pct": 100 * o[g[0]]["share_low"], "shared_5uM_pct": 100 * o[g[0]]["share_shared"], "only_5uM_pct": 100 * o[g[0]]["share_high_only"]} for g in gs]))
source("fig4c_mechanism_similarity", pd.DataFrame([{"dose_over_unbound_cmax": n, "median_auroc": be[n]["median_line_auroc"], "q1": be[n]["line_auroc_q1"], "q3": be[n]["line_auroc_q3"]} for n in names]))

# ------------------------------------------------ Figure 5: patient tumours
order = ["everolimus", "afatinib", "palbociclib", "anastrozole", "sunitinib", "anastrozole_b", "olaparib"]
nm = {"everolimus": "Everolimus", "afatinib": "Afatinib", "palbociclib": "Palbociclib", "anastrozole": "Anastrozole<br>(GSE153470)", "sunitinib": "Sunitinib", "anastrozole_b": "Anastrozole<br>(GSE93204)", "olaparib": "Olaparib"}
fig = go.Figure()
for key, lab, col, off in [("low", "0.05 µM signature", C[0], -0.15), ("high_only", "5 µM-only signature", C[1], 0.15)]:
    y = [B["cohorts"][k][key]["rate"] for k in order]; ci = [B["cohorts"][k][key]["ci"] for k in order]
    fig.add_trace(go.Scatter(x=[i + off for i in range(len(order))], y=[100 * v for v in y], mode="markers", marker=dict(color=col, symbol="diamond", size=11), name=lab,
                             error_y=dict(type="data", symmetric=False, array=[100 * (c_[1] - v) for c_, v in zip(ci, y)], arrayminus=[100 * (v - c_[0]) for c_, v in zip(ci, y)], color=K, thickness=1, width=4)))
fig.add_hline(y=50, line=dict(color="#888888", dash="dash", width=1)); fig.add_vline(x=4.5, line=dict(color="#bbbbbb", width=1))
fig.update_xaxes(tickvals=list(range(len(order))), ticktext=[f"{nm[k]}<br>n = {B['cohorts'][k]['n_pairs']}" for k in order]); fig.update_yaxes(title_text="Genes changing in the same direction in patient tumours (%)", range=[25, 75])
style(fig, 1000, 520).update_layout(legend=dict(x=0.01, y=0.99, bgcolor="rgba(255,255,255,0.7)"), margin=dict(b=90))
save(fig, "Figure5")
source("fig5_patient_gene_agreement", pd.DataFrame([{"cohort": k, "pairs": B["cohorts"][k]["n_pairs"], **{f"{key}_{q}": v for key in ("low", "high_only") for q, v in (("pct", 100 * B["cohorts"][k][key]["rate"]), ("ci_low_pct", 100 * B["cohorts"][k][key]["ci"][0]), ("ci_high_pct", 100 * B["cohorts"][k][key]["ci"][1]))}} for k in order]))

# ------------------------------------------------ Figure 6: power at clinical exposure and pathway-level agreement
pc = pd.read_csv(R + "supporting/power_curve.tsv", sep="\t"); P = S["power"]["<=1x"]
fig = make_subplots(rows=1, cols=3, subplot_titles=["a  90th percentile effect", "b  Mean effect", "c  Hallmark gene sets"], horizontal_spacing=0.1)
for j, e in enumerate(["p90", "mean"]):
    for r_, col in zip([2, 4, 8, 16, 32], [C[0], C[5], C[2], C[4], C[1]]):
        x = pc[(pc.effect == e) & (pc.wells == r_)]
        fig.add_trace(go.Scatter(x=x.cells_per_well, y=100 * x.median_power, mode="lines", line=dict(color=col, width=2), name=f"{r_} wells", showlegend=j == 0, legend="legend"), 1, j + 1)
    fig.add_hline(y=80, line=dict(color="#888888", dash="dash", width=1), row=1, col=j + 1); fig.add_vline(x=150, line=dict(color="#bbbbbb", width=1, dash="dot"), row=1, col=j + 1)
    fig.update_xaxes(type="log", title_text="Cells per well", tickvals=[10, 100, 1000, 10000, 100000], ticktext=["10", "10<sup>2</sup>", "10<sup>3</sup>", "10<sup>4</sup>", "10<sup>5</sup>"], row=1, col=j + 1)
    fig.update_yaxes(title_text="Median power (%)" if j == 0 else None, range=[0, 102], row=1, col=j + 1)
for key, lab, col, off in [("0.05uM", "0.05 µM profile", C[0], -0.15), ("5uM", "5 µM profile", C[1], 0.15)]:
    fig.add_trace(go.Scatter(x=[i + off for i in range(len(order))], y=[B["cohorts"][k]["pathway"][key]["rho"] for k in order], mode="markers", marker=dict(color=col, symbol="diamond", size=10), name=lab, legend="legend2"), 1, 3)
fig.add_hline(y=0, line=dict(color="#888888", dash="dash", width=1), row=1, col=3); fig.add_vline(x=4.5, line=dict(color="#bbbbbb", width=1), row=1, col=3)
short = {"everolimus": "Everolimus", "afatinib": "Afatinib", "palbociclib": "Palbociclib", "anastrozole": "Anastrozole", "sunitinib": "Sunitinib", "anastrozole_b": "Anastrozole (2)", "olaparib": "Olaparib"}
fig.update_xaxes(tickvals=list(range(len(order))), ticktext=[short[k] for k in order], tickangle=-40, tickfont=dict(size=11), row=1, col=3)
fig.update_yaxes(title_text="Rank correlation with patient change, 50 gene sets", range=[-0.7, 0.9], row=1, col=3)
style(fig, 1300, 520).update_layout(legend=dict(x=0.385, y=0.99, bgcolor="rgba(255,255,255,0.7)"), legend2=dict(x=0.86, y=0.99, font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"), margin=dict(b=100))
save(fig, "Figure6")
source("fig6ab_power_curves", pc)
source("fig6c_gene_set_correlation", pd.DataFrame([{"cohort": k, "rho_0.05uM": B["cohorts"][k]["pathway"]["0.05uM"]["rho"], "rho_5uM": B["cohorts"][k]["pathway"]["5uM"]["rho"]} for k in order]))


# ------------------------------------------------ Figure 7: response model by exposure tier and after removal of tiers from training
M, MI = json.load(open(R + "metrics/response_model.json")), json.load(open(R + "metrics/response_model_inert.json"))
fig = make_subplots(rows=1, cols=2, subplot_titles=["a  Prediction by exposure tier", "b  Removal of tiers from training"], horizontal_spacing=0.12)
gl = [("inert", "Inert"), ("tier_1", "Tier 1"), ("tier_2", "Tier 2"), ("tier_3", "Tier 3"), ("tier_4", "Tier 4")]; rows7 = []
for key, lab, col, off in [("r_model", "Full profile", C[2], -0.22), ("r_spec_all", "Minus context mean", C[0], 0.0), ("r_spec_plate_all", "Minus plate mean", C[1], 0.22)]:
    v = [MI["by_tier"][g][key] for g, _ in gl]
    fig.add_trace(go.Scatter(x=[i + off for i in range(len(gl))], y=[x["mean"] for x in v], mode="markers", marker=dict(color=col, symbol="diamond", size=10), name=lab, legend="legend",
        error_y=dict(type="data", symmetric=False, array=[x["hi"] - x["mean"] for x in v], arrayminus=[x["mean"] - x["lo"] for x in v], color=K, thickness=1.2, width=4)), 1, 1)
    rows7 += [{"panel": "a", "group": g, "measure": lab, "mean": x["mean"], "ci_low": x["lo"], "ci_high": x["hi"], "compounds": x["compounds"]} for (g, _), x in zip(gl, v)]
fig.update_xaxes(tickvals=list(range(len(gl))), ticktext=[l for _, l in gl], row=1, col=1); fig.update_yaxes(title_text="Correlation of predicted and observed profile", range=[0, 0.75], row=1, col=1)
hl = [("tier_1_2", "Held-out tier 1 and 2"), ("tier_3_4", "Held-out tier 3 and 4")]
for key, lab, col, off in [("r_spec_all", "All conditions", C[0], -0.22), ("r_spec_without_tier_1_2", "Without tier 1 and 2", C[4], 0.0), ("r_spec_without_tier_3_4", "Without tier 3 and 4", C[3], 0.22)]:
    v = [M["by_tier"][g][key] for g, _ in hl]
    fig.add_trace(go.Scatter(x=[i + off for i in range(len(hl))], y=[x["mean"] for x in v], mode="markers", marker=dict(color=col, symbol="diamond", size=10), name=lab, legend="legend2",
        error_y=dict(type="data", symmetric=False, array=[x["hi"] - x["mean"] for x in v], arrayminus=[x["mean"] - x["lo"] for x in v], color=K, thickness=1.2, width=4)), 1, 2)
    rows7 += [{"panel": "b", "group": g, "measure": "trained on: " + lab, "mean": x["mean"], "ci_low": x["lo"], "ci_high": x["hi"], "compounds": x["compounds"]} for (g, _), x in zip(hl, v)]
fig.update_xaxes(tickvals=[0, 1], ticktext=[l for _, l in hl], range=[-0.6, 1.6], row=1, col=2); fig.update_yaxes(title_text="Correlation minus context mean", range=[0, 0.75], row=1, col=2)
style(fig, 1100, 500).update_layout(legend=dict(x=0.01, y=0.02, yanchor="bottom", bgcolor="rgba(255,255,255,0.7)"), legend2=dict(x=0.57, y=0.02, yanchor="bottom", font=dict(color=K, size=12), bgcolor="rgba(255,255,255,0.7)"))
save(fig, "Figure7")
source("fig7_response_model", pd.DataFrame(rows7))


# ------------------------------------------------ Figure 1: study overview (hand-written SVG, readout bars from the results)
tr = next(x for x in T["sweep"] if x["cut"] == "1.0" and x["basis"] == "unbound" and x["pk_set"] == "liston" and x["stratum"] == "all")
sr = next(x for x in S["metrics"] if x["cut"] == 1.0 and x["basis"] == "unbound")
be = T["mechanism"]["by_exposure"]; pl = B["pooled"]
W, RH, X0, X1, X2 = 1300, 190, 20, 330, 800
o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{4 * RH + 50}" viewBox="0 0 {W} {4 * RH + 50}" font-family="Arial" font-size="14" fill="{K}">',
     f'<rect width="{W}" height="{4 * RH + 50}" fill="white"/>']
def t(x, y, s, size=14, w="normal", fill=K, anchor="start", style=""): o.append(f'<text x="{x}" y="{y}" font-size="{size}" font-weight="{w}" fill="{fill}" text-anchor="{anchor}" {style}>{s}</text>')
for x, s in [(X0 + 12, "Data"), (X1 + 12, "Comparison"), (X2 + 12, "Readout")]: t(x, 28, s, 16, "bold")
def bars(x, y, vals, labels, cols, maxv=100, h=100, bw=46, unit="%"):
    o.append(f'<line x1="{x - 6}" y1="{y + h}" x2="{x + len(vals) * (bw + 26)}" y2="{y + h}" stroke="{K}"/>')
    for i, (v, l, c) in enumerate(zip(vals, labels, cols)):
        bh = h * v / maxv; bx = x + i * (bw + 26)
        o.append(f'<rect x="{bx}" y="{y + h - bh}" width="{bw}" height="{bh}" fill="{c}"/>'); t(bx + bw / 2, y + h - bh - 5, f"{v:.0f}{unit}" if unit == "%" else f"{v:.2f}", 12, anchor="middle"); t(bx + bw / 2, y + h + 15, l, 11, anchor="middle")
def plate(x, y, cols):
    o.append(f'<rect x="{x}" y="{y}" width="112" height="64" rx="6" fill="none" stroke="{K}"/>')
    for r in range(3):
        for c in range(6): o.append(f'<circle cx="{x + 14 + c * 17}" cy="{y + 14 + r * 18}" r="6" fill="{cols[(r + c) % len(cols)]}" opacity="{0.35 + 0.3 * r}"/>')
def ladder(x, y, doses, cm, cols):
    o.append(f'<line x1="{x}" y1="{y}" x2="{x}" y2="{y + 120}" stroke="{K}"/>'); t(x - 8, y - 8, "dose", 11, anchor="middle")
    for i, (d, c) in enumerate(zip(doses, cols)): yy = y + 110 - i * (100 / (len(doses) - 1)); o.append(f'<circle cx="{x}" cy="{yy}" r="7" fill="{c}"/>'); t(x + 14, yy + 4, d, 12)
    o.append(f'<line x1="{x - 30}" y1="{cm}" x2="{x + 150}" y2="{cm}" stroke="{C[2]}" stroke-width="2" stroke-dasharray="6 4"/>'); t(x + 154, cm + 4, "unbound C<tspan baseline-shift='sub' font-size='9'>max</tspan> in patients", 12, fill=C[2])
rows = [(C[0], "Exposure calibration", "Fig. 2"), (C[1], "Response against exposure", "Fig. 3"), (C[3], "Specificity of the high-dose response", "Fig. 4"), (C[4], "Agreement with patient tumours", "Fig. 5")]
for i, (c, title, fg) in enumerate(rows):
    y = 44 + i * RH
    o.append(f'<rect x="{X0}" y="{y}" width="6" height="{RH - 14}" fill="{c}"/>'); t(X0 + 14, y + 16, title, 15, "bold", c); t(W - 24, y + 16, fg, 13, "bold", anchor="end")
    for x in (X1, X2): o.append(f'<line x1="{x}" y1="{y + 6}" x2="{x}" y2="{y + RH - 20}" stroke="#bbbbbb"/>')
    if i < 3: o.append(f'<line x1="{X0}" y1="{y + RH - 8}" x2="{W - 20}" y2="{y + RH - 8}" stroke="#dddddd"/>')
y = 44
plate(X0 + 16, y + 34, [C[0], C[1], C[2]]); t(X0 + 140, y + 50, "Tahoe-100M", 13, "bold"); t(X0 + 140, y + 68, "50 lines, 3 doses", 12); t(X0 + 140, y + 86, f"{tr['n_compounds']} compounds with PK", 12)
t(X0 + 16, y + 124, "sci-Plex 3", 13, "bold"); t(X0 + 100, y + 124, f"3 lines, 4 doses, {sr['n_compounds']} compounds", 12); t(X0 + 16, y + 146, "C<tspan baseline-shift='sub' font-size='9'>max</tspan> and protein binding: Liston and Davis 2017", 12)
ladder(X1 + 60, y + 36, ["0.05 µM", "0.5 µM", "5 µM"], y + 112, [C[0], C[1], C[2]])
bars(X2 + 40, y + 34, [100 * tr["share_contrasts_unreachable"], 100 * sr["share_contrasts_unreachable"]], ["Tahoe", "sci-Plex"], [C[0], C[1]]); t(X2 + 200, y + 70, "conditions above", 12); t(X2 + 200, y + 88, "unbound C<tspan baseline-shift='sub' font-size='9'>max</tspan>", 12)
y += RH
plate(X0 + 16, y + 34, [C[1], C[4]]); t(X0 + 140, y + 50, "DE genes per", 12); t(X0 + 140, y + 68, "compound, dose", 12); t(X0 + 140, y + 86, "and cell line", 12); t(X0 + 16, y + 124, "Background: inert compounds (Tahoe),", 12); t(X0 + 16, y + 142, "vehicle against vehicle (sci-Plex)", 12)
o.append(f'<path d="M {X1 + 40} {y + 140} C {X1 + 200} {y + 138}, {X1 + 280} {y + 120}, {X1 + 400} {y + 40}" fill="none" stroke="{C[1]}" stroke-width="3"/>')
o.append(f'<line x1="{X1 + 40}" y1="{y + 150}" x2="{X1 + 420}" y2="{y + 150}" stroke="{K}"/>'); o.append(f'<line x1="{X1 + 190}" y1="{y + 30}" x2="{X1 + 190}" y2="{y + 150}" stroke="{C[2]}" stroke-width="2" stroke-dasharray="6 4"/>')
t(X1 + 110, y + 166, "below", 12, anchor="middle"); t(X1 + 300, y + 166, "above unbound C<tspan baseline-shift='sub' font-size='9'>max</tspan>", 12, anchor="middle"); t(X1 + 40, y + 44, "response", 12)
bars(X2 + 40, y + 34, [100 * tr["de_share_above"], 100 * sr["de_share_above"]], ["Tahoe", "sci-Plex"], [C[0], C[1]]); t(X2 + 200, y + 70, "DE calls from doses", 12); t(X2 + 200, y + 88, "above unbound C<tspan baseline-shift='sub' font-size='9'>max</tspan>", 12)
y += RH
t(X0 + 16, y + 50, "Tahoe-100M DE gene sets", 12); t(X0 + 16, y + 68, "and log<tspan baseline-shift='sub' font-size='9'>2</tspan>FC profiles", 12); t(X0 + 16, y + 100, "Hallmark stress gene sets", 12); t(X0 + 16, y + 132, f"{T['mechanism']['n_drugs']} compounds, {len(T['mechanism']['classes'])} mechanism classes", 12)
for k, (cx, lab) in enumerate([(X1 + 90, "MEK inhibitor A"), (X1 + 250, "MEK inhibitor B"), (X1 + 410, "other mechanism")]):
    for j in range(8): hh = [30, 12, 44, 20, 8, 36, 16, 26] if k < 2 else [10, 40, 14, 34, 28, 8, 42, 18]; o.append(f'<rect x="{cx - 44 + j * 11}" y="{y + 100 - hh[j] - (4 if k == 1 and j % 2 else 0)}" width="8" height="{hh[j] + (4 if k == 1 and j % 2 else 0)}" fill="{C[3] if k < 2 else C[4]}"/>')
    o.append(f'<line x1="{cx - 50}" y1="{y + 100}" x2="{cx + 50}" y2="{y + 100}" stroke="{K}"/>'); t(cx, y + 118, lab, 12, anchor="middle")
t(X1 + 170, y + 146, "same mechanism", 12, anchor="middle", fill=C[3]); t(X1 + 330, y + 146, "different mechanism", 12, anchor="middle", fill=C[4])
bars(X2 + 40, y + 34, [be["<=1x"]["median_line_auroc"], be["10-100x"]["median_line_auroc"], be[">100x"]["median_line_auroc"]], ["≤1", "10 to 100", ">100"], [C[0], C[3], C[1]], maxv=1.0, unit="")
o.append(f'<line x1="{X2 + 30}" y1="{y + 84}" x2="{X2 + 250}" y2="{y + 84}" stroke="#888888" stroke-dasharray="5 4"/>'); t(X2 + 262, y + 70, "AUROC by dose /", 12); t(X2 + 262, y + 88, "unbound C<tspan baseline-shift='sub' font-size='9'>max</tspan>", 12); t(X2 + 262, y + 106, "(dashed: 0.5)", 12)
y += RH
def person(x, yy, c): o.append(f'<circle cx="{x}" cy="{yy}" r="9" fill="{c}"/><path d="M {x - 14} {yy + 42} Q {x} {yy + 2} {x + 14} {yy + 42} Z" fill="{c}"/>')
person(X0 + 40, y + 50, "#999999"); person(X0 + 130, y + 50, C[4]); o.append(f'<path d="M {X0 + 62} {y + 70} L {X0 + 104} {y + 70}" stroke="{K}" marker-end="url(#a)"/>'); t(X0 + 83, y + 62, "drug", 11, anchor="middle")
o.insert(2, f'<defs><marker id="a" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0,0 L8,4 L0,8 z" fill="{K}"/></marker></defs>')
t(X0 + 16, y + 118, "Paired biopsies before and on", 12); t(X0 + 16, y + 136, "treatment, 7 cohorts from GEO", 12); t(X0 + 16, y + 154, "(5 primary)", 12)
for k, (lab, c, ar) in enumerate([("in vitro", C[0], [1, -1, 1, 1, -1, 1]), ("patients", C[4], [1, 1, -1, 1, -1, -1])]):
    t(X1 + 40, y + 62 + k * 50, lab, 12)
    for j, a in enumerate(ar): xx = X1 + 130 + j * 48; yy = y + 58 + k * 50; o.append(f'<path d="M {xx} {yy + 12 * a} L {xx} {yy - 12 * a} M {xx - 6} {yy - 5 * a} L {xx} {yy - 12 * a} L {xx + 6} {yy - 5 * a}" stroke="{c}" stroke-width="3" fill="none"/>')
t(X1 + 250, y + 150, "same direction of change, gene by gene", 12, anchor="middle")
bars(X2 + 40, y + 34, [100 * pl["low_rate"], 100 * pl["high_only_rate"]], ["0.05 µM", "5 µM only"], [C[0], C[1]]); o.append(f'<line x1="{X2 + 30}" y1="{y + 84}" x2="{X2 + 180}" y2="{y + 84}" stroke="#888888" stroke-dasharray="5 4"/>')
t(X2 + 200, y + 70, "genes agreeing in", 12); t(X2 + 200, y + 88, "direction (dashed: 50%)", 12)
o.append("</svg>"); open(F + "fig1_overview.svg", "w").write("\n".join(o))
subprocess.run(["rsvg-convert", "-z", "3", "-o", F + "fig1_overview.png", F + "fig1_overview.svg"], check=True)
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", F + "fig1_overview.pdf", F + "fig1_overview.svg"], check=True)
source("fig1_readout_values", pd.DataFrame([{"row": "conditions above unbound Cmax (%)", "tahoe": 100 * tr["share_contrasts_unreachable"], "sciplex": 100 * sr["share_contrasts_unreachable"]}, {"row": "DE calls from doses above unbound Cmax (%)", "tahoe": 100 * tr["de_share_above"], "sciplex": 100 * sr["de_share_above"]},
    {"row": "mechanism AUROC at <=1x, 10-100x, >100x", "tahoe": f'{be["<=1x"]["median_line_auroc"]:.2f}, {be["10-100x"]["median_line_auroc"]:.2f}, {be[">100x"]["median_line_auroc"]:.2f}', "sciplex": ""},
    {"row": "genes agreeing with patients (%), 0.05 uM and 5 uM only", "tahoe": f'{100 * pl["low_rate"]:.1f}, {100 * pl["high_only_rate"]:.1f}', "sciplex": ""}]))
print("figures written")
