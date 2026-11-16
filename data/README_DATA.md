# Data sources

## Curated tables in this repository

| File | Content | Source |
|---|---|---|
| `pk/liston_davis_2017_pk.csv` | total and unbound Cmax for 67 Tahoe-100M entries (65 with protein binding) and 27 sci-Plex 3 compounds | Liston DR, Davis M. Clin Cancer Res 2017;23:3489-98 (PMC5511563), Tables 1 and 2 |
| `pk/fda_ema_label_extension.csv` | 51 Tahoe-100M entries (49 with differential expression results) with Cmax and protein binding from product labels | link to each label in the column `source_url` |
| `pk/compound_moa_annotations.csv` | all 379 Tahoe-100M compound entries with mechanism class, inert reference flag and analysis flags | Tahoe-100M compound metadata (`drug_metadata` of the Hugging Face repository) |
| `pk/compounds_de_rows.txt`, `pk/compounds_full_profiles.txt` | compounds whose DE rows and whose full profiles are kept by step 1 | this study |
| `gene_sets/hallmark.tsv`, `gene_sets/stress.tsv` | Hallmark gene sets, and the four stress sets among them | MSigDB 2024.1.Hs, `h.all.v2024.1.Hs.symbols.gmt` (CC BY 4.0) |
| `processed/sciplex3_pseudobulk/` | counts summed per well (`counts.bin.gz`, 32-bit floats, wells x genes in the order of `meta.tsv` and `genes.txt`) and the same after reduction of 10 uM wells to the cell number of the 10 nM wells (`counts_ds.bin.gz`) | step 3 on the sci-Plex 3 file below |
| `processed/tahoe_signatures.tsv.gz` | per compound, dose and gene: cell lines measured, cell lines with up- and down-regulation, median log2FC, for the 24 compounds with full profiles | step 2 on Tahoe-100M |

Unbound Cmax = Cmax x (1 - bound fraction). The label value of voriconazole could not be confirmed in a second
reading of the label. For crizotinib the label reports Cmax for paediatric dosing only, and the adult steady-state
value was used (see the column `note`).

## Public data outside the repository

| Data | Location | Local path |
|---|---|---|
| Tahoe-100M pseudobulk differential expression | Hugging Face dataset `tahoebio/Tahoe-100M`, folder `pseudobulk_differential_expression` of the Parquet conversion branch (`refs/convert/parquet`), 1,026 parts | read over HTTP by step 1, summaries in `work/tahoe_summary` |
| sci-Plex 3 single-cell counts | scPerturb, Zenodo record 7041849, `SrivatsanTrapnell2020_sciplex3.h5ad` | `raw/sciplex3.h5ad` |
| GSE119262, GSE192464, GSE93204, GSE58837 | GEO series matrix files | `raw/geo/<accession>_series_matrix.txt.gz` |
| GSE153470 | GEO supplementary file `GSE153470_Matrix_normalized_illu_GEX-FAIMoS81pairs.txt.gz` | `raw/geo/GSE153470_norm.txt.gz` |
| GSE341672 | GEO supplementary file `GSE341672_MEDIOLA_gene_counts.tsv.gz` | `raw/geo/GSE341672_counts.tsv.gz` |
| GPL6104, GPL6244, GPL6480, GPL10558 | GEO platform annotation files (`.annot.gz`) | `raw/geo/` |

`scripts/download_data.sh` fetches the GEO files and, with the argument `h5ad`, the sci-Plex 3 file.

## SHA-256 checksums

Files as used in this study (`data/raw_checksums.sha256`, relative to `raw/`). GEO platform annotation files are
rebuilt by GEO from time to time and may then differ.

```
bde2420ce24c8aad00d4b0fcbeb0351334e7cd6945d06030ce47be6bcc63b35f  sciplex3.h5ad
c914fdbe1130906ce3b9c97f5a75280591c89c617c168fb687c641f683e76b45  geo/GPL10558.annot.gz
82ae57d6d9ec26ce2bcff01ccd1db498bb8055ee01471829dec0a5ab5666d518  geo/GPL6104.annot.gz
99ee70f87e02c75ff1c8dc57b6c40c0482c7ef2760e004301315ac34344ac6dc  geo/GPL6244.annot.gz
aacec81b2d95f183f0827b43e769388a92917811b81920c13ecc8fe148c5a3ab  geo/GPL6480.annot.gz
f1f9de2c4d309c950ea391361ad49a5fb62cd4d1f9d665ca442c29e1af18cf78  geo/GSE119262_series_matrix.txt.gz
3c0daa1bdbd13c36519e828a2dce52fc60be09ffa9ba425899b5cc8cc08ad3fb  geo/GSE153470_norm.txt.gz
07e3c0a31f3e5bab327ca9ce265d6f52e851cf09a612f0bc2035f7c6c5c1520c  geo/GSE192464_series_matrix.txt.gz
67dc1d92ac58b469e4fc472fa49adb5dffcedb5d38b06fe2bbc44e794f440442  geo/GSE341672_counts.tsv.gz
254ef44286bb5e8bec7462aea7d1c56a22799cf38fb21d45abb77d74f263e62c  geo/GSE58837_series_matrix.txt.gz
158c07898de101c71b0b3a429721bcc81bbd9d5122d44a12aa2673972bf72aae  geo/GSE93204_series_matrix.txt.gz
792d88709065f0edb49cbf675bef89664588f02d79906821d99e1d87eaaa1b74  ../data/processed/tahoe_signatures.tsv.gz
e0db51935262deae5bae247b297c40357886a5184c9c5f6b58dbdb26bf09039b  ../data/processed/sciplex3_pseudobulk/counts.bin.gz
14769c7745b11fd5e2a42de23db89e3e7d8789b4e546aac0512e10f4593bdeea  ../data/processed/sciplex3_pseudobulk/counts_ds.bin.gz
```

The Tahoe-100M parts are read as a stream and no local copy exists, so no checksum is given for them.
