#!/bin/sh
# Runs every analysis and writes results/.
# Usage: ./run_all.sh [stream]
#   stream   also re-reads the Tahoe-100M differential expression table from Hugging Face (step 1)
# Inputs outside the repository (see data/README_DATA.md): work/tahoe_summary (step 1 output),
# raw/sciplex3.h5ad (optional, the per-well sums are in data/processed) and raw/geo (scripts/download_data.sh).
# Author: Thai Huu Loc Tran
set -e
cd "$(dirname "$0")"
B=target/release; W=work; R=results
cargo build --release
mkdir -p $W/tables $R/metrics $R/supporting
if [ "$1" = "stream" ]; then $B/tahoe_summary $W/tahoe_summary 0 1026 data/pk/compounds_de_rows.txt data/pk/compounds_full_profiles.txt; fi
PB=data/processed/sciplex3_pseudobulk
if [ -f raw/sciplex3.h5ad ]; then PB=$W/sciplex3_pseudobulk; [ -f $PB/counts.bin ] || $B/sciplex_de_power pseudobulk raw/sciplex3.h5ad data $PB; fi
SIG=data/processed/tahoe_signatures.tsv.gz
if [ -d $W/tahoe_summary ]; then
  $B/exposure_analysis $W/tahoe_summary data $W/tables $R/metrics/tahoe.json; SIG=$W/tables/tahoe_signatures.tsv
  for f in sweep strata recovery_per_compound drug_summary moa contrasts overlap_pairs; do cp $W/tables/$f.tsv $R/supporting/$f.tsv; done
else echo "work/tahoe_summary not found: Tahoe-100M exposure analysis skipped, committed results kept"; fi
$B/sciplex_de_power de $PB data $W/tables $R/metrics/sciplex.json
for f in sciplex_contrasts power_curve; do cp $W/tables/$f.tsv $R/supporting/$f.tsv; done
if [ -d raw/geo ]; then
  $B/patient_biopsies raw/geo $SIG data $W/tables $R/metrics/patients.json; cp $W/tables/pathway_scores.tsv $R/supporting/pathway_scores.tsv
else echo "raw/geo not found: patient cohort analysis skipped, committed results kept"; fi
$B/export_tiers data $R/supporting/contrasts.tsv $R/atlas_clinical_exposure_tiers
$B/export_tables $R
