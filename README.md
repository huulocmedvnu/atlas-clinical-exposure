# The detectable transcriptional response in single-cell chemical perturbation atlases arises above clinical drug exposure

Code, analysis plan, curated pharmacokinetic tables and results for the comparison of the doses used in
Tahoe-100M and sci-Plex 3 with the plasma concentrations reached in patients.

Author: Thai Huu Loc Tran, Department of Obstetrics and Gynecology, Hong Hung Hospital, Tay Ninh City 80116,
Viet Nam (tranthaihuuloc2@gmail.com, ORCID 0009-0009-1810-010X).

Version 1.2.0 is archived at Zenodo: https://doi.org/10.5281/zenodo.24299962

## Exposure tier table

`results/atlas_clinical_exposure_tiers.csv` (and `.parquet`) annotates every compound and dose with clinical
pharmacokinetic data: 348 conditions of 116 Tahoe-100M compound entries and 108 conditions of 27 sci-Plex 3
compounds. Two of the Tahoe-100M entries (selinexor and erdafitinib) have no differential expression results in the
atlas release and are flagged in the column `in_de_results`. Different salt forms of one substance are separate entries.

| Column | Content |
|---|---|
| `atlas`, `compound`, `generic_name` | atlas, compound entry as named in the atlas, generic name |
| `pk_source` | `Liston_Davis_2017` or `product_label` |
| `nominal_dose_uM` | dose applied in the atlas |
| `cmax_total_uM`, `f_unbound_plasma`, `cmax_unbound_uM` | clinical maximum plasma concentration, unbound fraction, unbound concentration |
| `f_unbound_medium_10pct_fbs` | unbound fraction in medium with 10% serum |
| `dose_to_cmax_unbound_ratio`, `dose_to_cmax_serum_adj_ratio` | dose over unbound Cmax, without and with the serum adjustment |
| `exposure_tier` | tier 1: dose at or below unbound Cmax. Tier 2: reachable after serum adjustment only. Tier 3: up to 10 times above. Tier 4: more than 10 times above |
| `w_pk` | optional weight min(1, (unbound Cmax / dose)^2). It has not been evaluated in model training |
| `in_de_results` | whether the compound has differential expression results in the atlas release |

## Licence

The code is released under the MIT licence (`LICENSE`). The data tables created in this study (`data/pk`, the exposure tier table,
the result tables and the figure source data in `results/`) are released under CC0 1.0 (`LICENSE-DATA`). Files derived from
third-party resources keep the terms of their sources, which are listed in `data/README_DATA.md`.

## Layout

- `ANALYSIS_PLAN.md` definitions, thresholds and decision rules, with the outcome of every rule.
- `src/lib.rs` shared statistics and input readers. `src/bin/` one program per step:
  - `01_tahoe_summary.rs` reads the Tahoe-100M differential expression table from Hugging Face and keeps per-contrast summaries.
  - `02_exposure_analysis.rs` exposure calibration, serum adjustment, background of inert compounds, overlap between mechanisms, mechanism similarity.
  - `03_sciplex_de_power.rs` sci-Plex 3 per-well sums, negative binomial differential expression, cell-number control, power at clinical exposure.
  - `04_patient_biopsies.rs` paired patient cohorts: agreement of genes, of the most changed genes and of Hallmark gene set scores, and sensitivity of the patient profiles.
  - `05_export_tiers.rs` exposure tier table. `06_export_tables.rs` values of Tables 1 to 3 and `final_metrics.json`.
  - `07_tahoe_profiles.rs` (additional analysis, not reported in the article) log2 fold changes of the Hallmark genes for every Tahoe-100M contrast, read from Hugging Face.
  - `08_single_cell_cohorts.rs` (additional analysis, not reported in the article) single-cell patient cohorts: counts summed over malignant cells per sample and agreement with the atlas signatures.
- `scripts/response_model.py` (additional analysis, not reported in the article) response model by exposure tier and after removal of tiers from training (PyTorch).
- `scripts/plot_figures.py` Figures 1 to 6 and their source data. `scripts/download_data.sh` public input files.
  `scripts/match_liston_davis.py` and `scripts/curate_pk_tables.py` document how the tables in `data/pk` were derived.
- `data/pk/` curated Cmax and protein binding, mechanism classes and inert reference compounds. `data/gene_sets/` Hallmark gene sets.
  `data/processed/` compact inputs. `data/README_DATA.md` sources and checksums.
- `results/` tier table, `table1` to `table3` CSV files, `figure_source_data/`, `metrics/` (all statistics per analysis),
  `final_metrics.json` (the three together) and `supporting/` (per-compound and per-contrast tables).
- `figures/` Figures 1 to 6 as PDF.

## Running

Requirements: a recent stable Rust toolchain (tested with 1.99) and `cmake` (HDF5 is compiled from source). Figures need Python with `plotly`, `kaleido`,
`pandas`, and `rsvg-convert`. The code runs on macOS (Apple silicon) and Linux. Only the additional response model needs PyTorch.

```
./run_all.sh
python3 scripts/plot_figures.py
```

`run_all.sh` uses what is present and keeps the committed results for steps whose inputs are absent:

| Step | Input | In the repository |
|---|---|---|
| sci-Plex 3 differential expression and power | per-well sums, `data/processed/sciplex3_pseudobulk` | yes |
| exposure tiers, Tables 1 to 3 | `data/pk`, `results/` | yes |
| patient cohorts | `raw/geo` (about 120 MB) and `data/processed/tahoe_signatures.tsv.gz` | after `scripts/download_data.sh` |
| Tahoe-100M exposure analysis | `work/tahoe_summary` (about 1.1 GB) | after `./run_all.sh stream` |
| additional analysis: malignant cells of single-cell patient cohorts | sums per sample, `data/processed/single_cell` | yes (rebuilt from raw files by `./run_all.sh cells`) |
| additional analysis: response model | `work/tahoe_profiles` (about 1.4 GB) and PyTorch | after `./run_all.sh model` |
| sci-Plex 3 per-well sums from single cells | `raw/sciplex3.h5ad` (2.5 GB) | after `scripts/download_data.sh h5ad` |

`./run_all.sh stream` reads the 1,026 parts of the Tahoe-100M table over HTTP and stores only the summaries.
Confidence intervals use fixed random streams, and repeated runs give the same values on both platforms up to
floating-point rounding. The order of rows in some supporting tables can differ between runs.

## Licence

MIT for the code (see `LICENSE`). The curated tables are derived from the cited sources (see `data/README_DATA.md`).
