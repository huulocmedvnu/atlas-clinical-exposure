#!/usr/bin/env python3
"""Response model and exposure (analysis plan, section 10).

Trains a small network that predicts the Tahoe-100M log2 fold change profile of a held-out
compound and cell line pair, and reports the compound-specific prediction by exposure tier and
after removal of tiers from training.

Usage: response_model.py <profiles dir> <tier table csv> <out json> <out per-contrast tsv> [compound annotation csv]
With the annotation table, the inert reference compounds are added as tier 0 and only the model trained on
all conditions is run (exploratory).

Author: Thai Huu Loc Tran
"""
import json
import sys

import numpy as np
import pandas as pd
import torch

PROF, TIERS, OUT_JSON, OUT_TSV = sys.argv[1:5]
INERT = sys.argv[5] if len(sys.argv) > 5 else None
HIDDEN, EPOCHS, BATCH, LR, FOLDS, BOOT = 512, 60, 256, 1e-3, 5, 2000
dev = "cuda" if torch.cuda.is_available() else "cpu"
torch.backends.cudnn.deterministic = True

# ---------------------------------------------------------------- data
genes = [g for g in open(f"{PROF}/genes.txt").read().split("\n") if g]
meta = pd.read_csv(f"{PROF}/meta.tsv", sep="\t")
lfc = np.fromfile(f"{PROF}/lfc.f32", dtype="<f4").reshape(len(meta), len(genes))
tt = pd.read_csv(TIERS)
tt = tt[(tt.atlas == "Tahoe-100M") & (tt.in_de_results.astype(str).str.lower() == "true")]
tier = {(c, round(float(d), 4)): int(t.split("_")[1]) for c, d, t in zip(tt.compound, tt.nominal_dose_uM, tt.exposure_tier)}
meta["dose"] = meta.dose.round(4)
if INERT:
    an = pd.read_csv(INERT)
    for c in an.compound[an.inert_reference == 1]:
        for d in (0.05, 0.5, 5.0):
            tier.setdefault((c, d), 0)
keep = np.array([(c, d) in tier for c, d in zip(meta.drug, meta.dose)])
meta_all, lfc_all = meta, np.clip(lfc, -5, 5)
meta, lfc = meta[keep].reset_index(drop=True), lfc_all[keep]
ok = np.isnan(lfc).mean(0) <= 0.20
lfc, lfc_all, genes = lfc[:, ok], lfc_all[:, ok], [g for g, k in zip(genes, ok) if k]
plates = meta.groupby(["drug", "dose", "cell_line"], sort=True).plate.agg(lambda x: tuple(sorted(set(x))))
grp = meta.groupby(["drug", "dose", "cell_line"], sort=True).indices
con = pd.DataFrame(list(grp.keys()), columns=["drug", "dose", "cell_line"])
with np.errstate(all="ignore"):
    Y = np.stack([np.nanmean(lfc[i], 0) for i in grp.values()]).astype(np.float32)
con["tier"] = [tier[(c, d)] for c, d in zip(con.drug, con.dose)]
drugs, lines, doses = sorted(con.drug.unique()), sorted(con.cell_line.unique()), sorted(con.dose.unique())
ci = con.drug.map({d: i for i, d in enumerate(drugs)}).to_numpy()
li = con.cell_line.map({d: i for i, d in enumerate(lines)}).to_numpy()
di = con.dose.map({d: i for i, d in enumerate(doses)}).to_numpy()
pair = ci * len(lines) + li
rng = np.random.default_rng(0)
up = np.unique(pair)
fold_of = dict(zip(rng.permutation(up), np.arange(len(up)) % FOLDS))
fold = np.array([fold_of[p] for p in pair])
supra = con.tier.to_numpy() >= 3
print(f"{len(con)} contrasts, {len(drugs)} compounds, {len(lines)} cell lines, {len(genes)} genes, device {dev}", flush=True)


# ---------------------------------------------------------------- model
class Net(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.c, self.l, self.d = (torch.nn.Embedding(n, HIDDEN) for n in (len(drugs), len(lines), len(doses)))
        self.f = torch.nn.Sequential(torch.nn.ReLU(), torch.nn.Linear(HIDDEN, HIDDEN), torch.nn.ReLU(), torch.nn.Linear(HIDDEN, HIDDEN),
                                     torch.nn.ReLU(), torch.nn.Linear(HIDDEN, len(genes)))

    def forward(self, c, l, d):
        return self.f(self.c(c) + self.l(l) + self.d(d))


Yt = torch.tensor(Y, device=dev)
Mt = ~torch.isnan(Yt)
Y0 = torch.nan_to_num(Yt)
C, L, D = (torch.tensor(x, device=dev, dtype=torch.long) for x in (ci, li, di))


def fit_predict(train, test, seed):
    torch.manual_seed(seed)
    net = Net().to(dev)
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    tr = torch.tensor(np.flatnonzero(train), device=dev)
    for _ in range(EPOCHS):
        for b in tr[torch.randperm(len(tr), device=dev)].split(BATCH):
            err = (net(C[b], L[b], D[b]) - Y0[b]) ** 2 * Mt[b]
            loss = err.sum() / Mt[b].sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
    net.eval()
    te = torch.tensor(np.flatnonzero(test), device=dev)
    with torch.no_grad():
        return torch.cat([net(C[b], L[b], D[b]) for b in te.split(4096)]).cpu().numpy()


def corr(a, b, m):
    a, b = np.where(m, a, 0.0), np.where(m, b, 0.0)
    n = m.sum(1)
    a, b = a - (a.sum(1) / n)[:, None] * m, b - (b.sum(1) / n)[:, None] * m
    return (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))


# Exploratory (not in the plan): mean profile of the other compounds on the same plate, cell line and dose.
pgrp = meta_all.groupby(["plate", "cell_line", "dose"]).indices
drug_all = meta_all.drug.to_numpy()


def plate_mean(rows, test_drugs_by_line):
    out = np.zeros((len(rows), Y.shape[1]), np.float32)
    for j, i in enumerate(rows):
        c, d, l = con.drug[i], con.dose[i], con.cell_line[i]
        ix = np.concatenate([pgrp[(p, l, d)] for p in plates.iloc[i]])
        ix = ix[~np.isin(drug_all[ix], test_drugs_by_line[l])]
        if len(ix):
            with np.errstate(all="ignore"):
                out[j] = np.nan_to_num(np.nanmean(lfc_all[ix], 0))
    return out


ARMS = {"all": np.ones(len(con), bool), "without_tier_3_4": ~supra, "without_tier_1_2": supra}
if INERT:
    ARMS = {"all": ARMS["all"]}
res = {a: np.full(len(con), np.nan) for a in ARMS}
r_plain, r2, r_ctx, r_plate = (np.full(len(con), np.nan) for _ in range(4))
res_p = {a: np.full(len(con), np.nan) for a in ARMS}
for k in range(FOLDS):
    te = fold == k
    ctx = np.zeros_like(Y[te])
    with np.errstate(all="ignore"):
        for l_ in range(len(lines)):
            for d_ in range(len(doses)):
                ctx[(li[te] == l_) & (di[te] == d_)] = np.nan_to_num(np.nanmean(Y[~te & (li == l_) & (di == d_)], 0))
    obs, m = np.nan_to_num(Y[te]), ~np.isnan(Y[te])
    r_ctx[te] = corr(ctx, obs, m)
    rows = np.flatnonzero(te)
    by_line = {l_: set(con.drug[rows][con.cell_line[rows] == l_]) for l_ in lines}
    pm = plate_mean(rows, {l_: np.array(sorted(v)) for l_, v in by_line.items()})
    r_plate[te] = corr(pm, obs, m)
    for a, sel in ARMS.items():
        pred = fit_predict(~te & sel, te, 100 * k + list(ARMS).index(a))
        res[a][te] = corr(pred - ctx, obs - ctx, m)
        res_p[a][te] = corr(pred - pm, obs - pm, m)
        if a == "all":
            r_plain[te] = corr(pred, obs, m)
            r2[te] = 1 - (((pred - obs) ** 2) * m).sum(1) / ((obs ** 2) * m).sum(1)
        print(f"fold {k} {a}: r_spec {np.nanmean(res[a][te]):.3f}", flush=True)

out = con.copy()
out["fold"], out["r_context_mean"], out["r_model"], out["r2_model"], out["r_plate_mean"] = fold, r_ctx, r_plain, r2, r_plate
for a in ARMS:
    out[f"r_spec_{a}"] = res[a]
    out[f"r_spec_plate_{a}"] = res_p[a]
out["sum_sq"] = np.nansum(Y ** 2, 1)
out.round(5).to_csv(OUT_TSV, sep="\t", index=False)

# ---------------------------------------------------------------- summaries
GROUPS = {"inert": [0], "tier_1": [1], "tier_2": [2], "tier_3": [3], "tier_4": [4], "tier_1_2": [1, 2], "tier_3_4": [3, 4], "all": [1, 2, 3, 4]}
brng = np.random.default_rng(1)
idx = brng.integers(0, len(drugs), (BOOT, len(drugs)))


def per_compound(col, tiers):
    s = out[out.tier.isin(tiers)].groupby("drug")[col].mean()
    return s.reindex(drugs).to_numpy()


def ci_mean(v):
    b = np.nanmean(v[idx], 1)
    return {"mean": float(np.nanmean(v)), "lo": float(np.nanpercentile(b, 2.5)), "hi": float(np.nanpercentile(b, 97.5)), "compounds": int(np.isfinite(v).sum())}


def ci_diff(a, b):
    d = np.nanmean(a[idx], 1) - np.nanmean(b[idx], 1)
    return {"difference": float(np.nanmean(a) - np.nanmean(b)), "lo": float(np.nanpercentile(d, 2.5)), "hi": float(np.nanpercentile(d, 97.5))}


S = {"contrasts": len(con), "compounds": len(drugs), "cell_lines": len(lines), "genes": len(genes),
     "contrasts_by_tier": {str(t): int((con.tier == t).sum()) for t in (1, 2, 3, 4)}, "by_tier": {}}
for g, ts in GROUPS.items():
    S["by_tier"][g] = {c: ci_mean(per_compound(c, ts)) for c in ["r_context_mean", "r_plate_mean", "r_model", "r2_model"] + [f"r_spec_{a}" for a in ARMS] + [f"r_spec_plate_{a}" for a in ARMS]}
    S["by_tier"][g]["contrasts"] = int(out.tier.isin(ts).sum())

ss = out.groupby(["drug", out.tier >= 3])["sum_sq"].sum().unstack(fill_value=0).reindex(drugs).to_numpy()
share = lambda m: m[:, 1].sum() / m.sum()
bs = np.array([share(ss[i]) for i in idx])
cond = float(supra.mean())
S["R1_training_signal"] = {"share_of_conditions_tier_3_4": cond, "share_of_squared_change_tier_3_4": float(share(ss)),
                           "lo": float(np.percentile(bs, 2.5)), "hi": float(np.percentile(bs, 97.5)), "supported": bool(np.percentile(bs, 2.5) > cond)}
d41 = ci_diff(per_compound("r_spec_all", [4]), per_compound("r_spec_all", [1]))
t1, t4 = S["by_tier"]["tier_1"]["r_spec_all"], S["by_tier"]["tier_4"]["r_spec_all"]
S["R2_skill_by_exposure"] = {"tier_4_minus_tier_1": d41, "present_in_tier_1": bool(t1["lo"] > 0), "present_in_tier_4": bool(t4["lo"] > 0),
                             "confined_to_supraphysiological_supported": bool(d41["lo"] > 0)}


def paired(arm, ts):
    if INERT:
        return {}
    v = per_compound("r_spec_all", ts) - per_compound(f"r_spec_{arm}", ts)
    r = ci_mean(v)
    r["supported"] = bool(r["lo"] > 0)
    return r


S["exploratory_inert"] = {g + "_minus_inert": ci_diff(per_compound("r_spec_plate_all", GROUPS[g]), per_compound("r_spec_plate_all", [0])) for g in ("tier_1", "tier_4")} if INERT else {}
S["R3_ablation"] = {} if INERT else {"clinical_removed_effect_on_tier_1_2": paired("without_tier_1_2", [1, 2]),
                    "supraphysiological_removed_effect_on_tier_3_4": paired("without_tier_3_4", [3, 4]),
                    "supraphysiological_removed_effect_on_tier_1_2": paired("without_tier_3_4", [1, 2]),
                    "clinical_removed_effect_on_tier_3_4": paired("without_tier_1_2", [3, 4])}
S["exploratory_plate_adjusted"] = {"tier_4_minus_tier_1": ci_diff(per_compound("r_spec_plate_all", [4]), per_compound("r_spec_plate_all", [1])),
    "tier_3_4_minus_tier_1_2": ci_diff(per_compound("r_spec_plate_all", [3, 4]), per_compound("r_spec_plate_all", [1, 2]))}
json.dump(S, open(OUT_JSON, "w"), indent=1)
print("done", flush=True)
