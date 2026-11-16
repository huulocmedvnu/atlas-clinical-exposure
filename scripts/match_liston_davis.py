"""Matches atlas compounds to the clinical Cmax and protein binding tables of Liston and Davis 2017 by generic name.

Usage: python3 scripts/match_liston_davis.py <dir>
The dir holds liston_pmc.html (the article page, PMC5511563) and drug_metadata.tsv (compound names of the atlas).
Writes cmax_matched.tsv to the same dir.

Author: Thai Huu Loc Tran
"""
import re, io, sys, pandas as pd, numpy as np
d = sys.argv[1]
t = pd.read_html(io.StringIO(open(f"{d}/liston_pmc.html").read()))
L = pd.concat([x for x in t if "Protein binding" in x.columns and "Cmax (uM)" in x.columns], ignore_index=True)
def norm(s):
    s = str(s).lower().replace("*", "")
    s = re.sub(r"\(.*?\)", "", s)
    s = re.sub(r"\b(hydrochloride|dihydrochloride|mesylate|dimesylate|tosylate|maleate|malate|citrate|sodium|acetate|phosphate|sulfate|tartrate|hydrate|monohydrate|besylate|fumarate|hemifumarate|salt|free base|di|hcl)\b", "", s)
    return re.sub(r"[^a-z0-9]", "", s)
def fu(pb):
    s = str(pb).lower()
    if "negligible" in s or "low" == s.strip(): return 1.0
    nums = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", s)]
    if not nums: return np.nan
    b = max(nums) if ">" in s else np.mean(nums)
    return max(1 - b / 100, 0.005)
L["key"] = L["Generic Name"].map(norm)
L["cmax_uM"] = pd.to_numeric(L["Cmax (uM)"], errors="coerce")
L["fu"] = L["Protein binding"].map(fu)
L = L.dropna(subset=["cmax_uM"]).drop_duplicates("key")
D = pd.read_csv(f"{d}/drug_metadata.tsv", sep="\t")
D["key"] = D["drug"].map(norm)
M = D.merge(L[["key", "Generic Name", "cmax_uM", "fu", "Protein binding"]], on="key", how="inner")
M["cmax_unbound_uM"] = M["cmax_uM"] * M["fu"]
M.to_csv(f"{d}/cmax_matched.tsv", sep="\t", index=False)
print(len(D), "Tahoe drugs;", len(L), "Liston rows;", len(M), "matched;", M["fu"].notna().sum(), "with binding")
print(M[["drug", "cmax_uM", "fu", "cmax_unbound_uM"]].to_string())
