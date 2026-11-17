#!/bin/sh
# Runs every analysis and writes results/.
# Usage: ./run_all.sh [stream|model|cells]
#   stream   also re-reads the Tahoe-100M differential expression table from Hugging Face (step 1)
#   cells    rebuilds the sums over malignant cells of the single-cell patient cohorts (raw/feline, raw/osimertinib)
#            and the ribociclib signature from Hugging Face (step 8 inputs)
#   model    reads the gene panel profiles from Hugging Face (step 7) and trains the response model (needs PyTorch)
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
SC=data/processed/single_cell
if [ "$1" = "cells" ]; then
  $B/tahoe_summary $W/tahoe_ribociclib 0 1026 $SC/compounds.txt $SC/compounds.txt
  $B/sc_cohorts signature $W/tahoe_ribociclib Ribociclib $W/ribociclib_signature.tsv && gzip -c $W/ribociclib_signature.tsv > $SC/ribociclib_signature.tsv.gz
  $B/sc_cohorts feline raw/feline/GSE158724_10x.counts.txt.gz $W/feline_sums.tsv && gzip -c $W/feline_sums.tsv > $SC/feline_sums.tsv.gz
  $B/sc_cohorts tenx raw/osimertinib/GSE292469_counts.h5 raw/osimertinib/GSE292469_cell_metadata.tsv.gz $W/osimertinib_sums.tsv && gzip -c $W/osimertinib_sums.tsv > $SC/osimertinib_sums.tsv.gz
fi
OSI=-; [ -f $SC/osimertinib_sums.tsv.gz ] && OSI=$SC/osimertinib_sums.tsv.gz
$B/sc_cohorts compare data $SC/ribociclib_signature.tsv.gz $SC/feline_sums.tsv.gz $OSI $R/metrics/single_cell_cohorts.json
$B/export_tiers data $R/supporting/contrasts.tsv $R/atlas_clinical_exposure_tiers
$B/export_tables $R
if [ "$1" = "model" ]; then
  [ -f $W/tahoe_profiles/lfc.f32 ] || $B/tahoe_profiles $W/tahoe_profiles 0 1026 data/gene_sets/hallmark.tsv
  python3 scripts/response_model.py $W/tahoe_profiles $R/atlas_clinical_exposure_tiers.csv $R/metrics/response_model.json $R/supporting/response_model_by_contrast.tsv
  python3 scripts/response_model.py $W/tahoe_profiles $R/atlas_clinical_exposure_tiers.csv $R/metrics/response_model_inert.json $R/supporting/response_model_inert_by_contrast.tsv data/pk/compound_moa_annotations.csv
fi
