"""Curation of the tables in data/pk and data/gene_sets from the source tables (run once, not part of run_all.sh).

Usage: python3 scripts/curate_pk_tables.py <source dir> <data dir>
The source dir holds: drug_metadata.tsv (Tahoe-100M compound metadata), cmax_matched.tsv (output of
match_liston_davis.py for Tahoe-100M), cmax_matched_sciplex.tsv (the same for sci-Plex 3), pk_expansion.tsv
(values read from product labels, with source links) and h.all.gmt (MSigDB Hallmark gene sets).
Unbound Cmax = Cmax x (1 - bound fraction). A range of binding is replaced by its mean, a lower bound by the bound,
and the unbound fraction is not allowed below 0.005.

Author: Thai Huu Loc Tran
"""
import sys, csv
SRC, DATA = sys.argv[1], sys.argv[2]
import re, numpy as np, pandas as pd
D = pd.read_csv(f"{SRC}/drug_metadata.tsv", sep="\t")
M = pd.read_csv(f"{SRC}/cmax_matched.tsv", sep="\t")[["drug", "Generic Name", "cmax_uM", "cmax_unbound_uM"]]; M["pk_source"] = "liston"
def fu(pb):
    s = str(pb).lower(); nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", s)]
    if not nums: return np.nan
    return max(1 - (max(nums) if ">" in s else np.mean(nums)) / 100, 0.005)
E = pd.read_csv(f"{SRC}/pk_expansion.tsv", sep="\t")
DROP = {"Rapamycin (sirolimus)": "whole blood", "Cyclosporin A": "whole blood", "Futibatinib": "unit misprint in label",
        "crizotinib": "paediatric dose", "Mebendazole": "upper limit only", "Piroxicam": "range only",
        "Fostamatinib (R406 active metabolite)": "dose not stated", "Clopidogrel": "active metabolite, no binding"}
E = E[~E.drug.isin(DROP)].copy(); E["v"] = pd.to_numeric(E.cmax_value, errors="coerce")
conv = {"ng/mL": lambda v, mw: v / mw, "ug/mL": lambda v, mw: v * 1000 / mw, "nM": lambda v, mw: v / 1000, "uM": lambda v, mw: v}
E["cmax_uM"] = [conv[u](v, mw) if u in conv and np.isfinite(v) else np.nan for v, u, mw in zip(E.v, E.cmax_unit, E.molecular_weight_g_per_mol)]
E["fu"] = E.protein_binding_percent.map(fu); E = E.dropna(subset=["cmax_uM", "fu"]); E["cmax_unbound_uM"] = E.cmax_uM * E.fu
key = lambda s: re.sub(r"[^a-z0-9]", "", re.sub(r"\(.*?\)", "", str(s)).split()[0].lower()) if str(s).strip() else ""
E["key"] = E.drug.map(lambda s: "monomethylfumarate" if "fumarate" in s.lower() else ("azacytidine" if "Azacyt" in s else key(s)))
T0 = D[~D.drug.isin(M.drug) & ~D.drug.str.startswith("(R)-")][["drug"]].copy()
T0["key"] = T0.drug.map(lambda s: "monomethylfumarate" if s.startswith("Dimethyl fumarate") else ("azacytidine" if "Azacyt" in s else key(s)))
X = T0.merge(E[["key", "generic_name", "cmax_uM", "cmax_unbound_uM"]], on="key").rename(columns={"generic_name": "Generic Name"}).drop(columns="key"); X["pk_source"] = "label"
cz = pd.DataFrame({"drug": ["crizotinib", "Crizotinib (hydrochloride)"], "Generic Name": "Crizotinib", "cmax_uM": 0.913, "cmax_unbound_uM": 0.913 * 0.09, "pk_source": "label"})
PK = pd.concat([M, X, cz]).drop_duplicates("drug")
INERT = {"Talc", "Allantoin", "Folic acid", "Fumaric acid", "Arbutin", "Gallic acid", "Sodium Salicylate", "Lactate (calcium)"}
FULL = set(open(f"{DATA}/pk/compounds_full_profiles.txt").read().split("\n"))
T = D[["drug", "moa-fine"]].merge(PK, on="drug", how="left").rename(columns={"moa-fine": "moa", "Generic Name": "generic"})
T["generic"] = T.generic.fillna(T.drug).map(lambda s: re.sub(r"[\t*]", "", str(s)).strip()); T["pk_source"] = T.pk_source.fillna("none")
T["inert"] = T.drug.isin(INERT).astype(int); T["full"] = T.drug.isin(FULL).astype(int)
u = T[(T.full == 1) & (T.moa != "unclear") & (T.drug != "Trametinib (DMSO_TF solvate)")]
T["moa_use"] = T.drug.isin(u[u.moa.map(u.moa.value_counts()) >= 2].drug).astype(int)
S = pd.read_csv(f"{SRC}/cmax_matched_sciplex.tsv", sep="\t").dropna(subset=["cmax_unbound_uM"])
# compound annotations
A = T[["drug", "generic", "moa", "inert", "full", "moa_use"]].copy(); A.columns = ["compound", "generic_name", "mechanism_class", "inert_reference", "full_profile", "mechanism_analysis"]
A.to_csv(f"{DATA}/pk/compound_moa_annotations.csv", index=False)
# compilation of Liston and Davis 2017, both atlases
cols = ["atlas", "compound", "generic_name", "cmax_total_uM", "cmax_unbound_uM"]
a = T[T.pk_source == "liston"][["drug", "generic", "cmax_uM", "cmax_unbound_uM"]].copy(); a.insert(0, "atlas", "Tahoe-100M"); a.columns = cols
b = S[["drug", "Generic Name", "cmax_uM", "cmax_unbound_uM"]].copy(); b.insert(0, "atlas", "sci-Plex 3"); b.columns = cols
pd.concat([a, b]).to_csv(f"{DATA}/pk/liston_davis_2017_pk.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)
# product label extension, one row per atlas entry
e = T[T.pk_source == "label"][["drug", "generic", "cmax_uM", "cmax_unbound_uM"]].copy(); e["k"] = e.generic.str.lower(); E["k"] = E.generic_name.str.lower()
m = e.merge(E[["k", "v", "cmax_unit", "cmax_context", "protein_binding_percent", "molecular_weight_g_per_mol", "source_url"]], on="k", how="left").drop(columns="k")
m.columns = ["compound", "generic_name", "cmax_total_uM", "cmax_unbound_uM", "cmax_reported", "cmax_unit", "dosing_context", "protein_binding_percent", "molecular_weight_g_per_mol", "source_url"]
m["note"] = ""; i = m.generic_name == "Crizotinib"
m.loc[i, ["cmax_reported", "cmax_unit", "dosing_context", "protein_binding_percent", "molecular_weight_g_per_mol"]] = [411.0, "ng/mL", "250 mg twice daily, adults, steady state", "91%", 450.34]
m.loc[i, "source_url"] = "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=2a51b0de-47d6-455e-a94c-d2c737b04ff7"
m.loc[i, "note"] = "Protein binding from the label. The label reports Cmax for paediatric dosing only. The adult steady-state Cmax is the value of the regulatory clinical pharmacology review."
m.to_csv(f"{DATA}/pk/fda_ema_label_extension.csv", index=False, quoting=csv.QUOTE_NONNUMERIC)
G = ["HALLMARK_APOPTOSIS", "HALLMARK_P53_PATHWAY", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"]
with open(f"{DATA}/gene_sets/stress.tsv", "w") as f, open(f"{DATA}/gene_sets/hallmark.tsv", "w") as h:
    for l in open(f"{SRC}/h.all.gmt"):
        p = l.strip().split("\t")
        for g in p[2:]:
            h.write(f"{p[0]}\t{g}\n")
            if p[0] in G: f.write(f"{p[0]}\t{g}\n")
print(T.pk_source.value_counts().to_dict(), "sciplex", len(S)); print(sorted(T[T.pk_source == "label"].drug))
