#!/bin/sh
# Downloads the public input files that are not part of the repository into raw/.
# Usage: scripts/download_data.sh [h5ad]   ("h5ad" also fetches the sci-Plex 3 file, 2.5 GB)
# Author: Thai Huu Loc Tran
set -e
cd "$(dirname "$0")/.."
mkdir -p raw/geo
G=https://ftp.ncbi.nlm.nih.gov/geo
get() { [ -f "raw/geo/$2" ] || curl -fsSL "$1" -o "raw/geo/$2"; }
get $G/series/GSE119nnn/GSE119262/matrix/GSE119262_series_matrix.txt.gz GSE119262_series_matrix.txt.gz
get $G/series/GSE192nnn/GSE192464/matrix/GSE192464_series_matrix.txt.gz GSE192464_series_matrix.txt.gz
get $G/series/GSE93nnn/GSE93204/matrix/GSE93204_series_matrix.txt.gz GSE93204_series_matrix.txt.gz
get $G/series/GSE58nnn/GSE58837/matrix/GSE58837_series_matrix.txt.gz GSE58837_series_matrix.txt.gz
get $G/series/GSE153nnn/GSE153470/suppl/GSE153470_Matrix_normalized_illu_GEX-FAIMoS81pairs.txt.gz GSE153470_norm.txt.gz
get $G/series/GSE341nnn/GSE341672/suppl/GSE341672_MEDIOLA_gene_counts.tsv.gz GSE341672_counts.tsv.gz
for p in GPL6104 GPL6244 GPL6480; do get $G/platforms/GPL6nnn/$p/annot/$p.annot.gz $p.annot.gz; done
get $G/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz GPL10558.annot.gz
if [ "$1" = "h5ad" ] && [ ! -f raw/sciplex3.h5ad ]; then
  curl -fL https://zenodo.org/api/records/7041849/files/SrivatsanTrapnell2020_sciplex3.h5ad/content -o raw/sciplex3.h5ad
fi
(cd raw && shasum -a 256 -c ../data/raw_checksums.sha256 --ignore-missing)
