# Analysis plan

Definitions, thresholds and decision rules of the study. Each section was written before the
corresponding results were produced. Sections 1 to 8 were fixed on 16 November 2026 and
sections 10 to 12 on 17 November 2026. Section 9 lists the outcome of every rule.
Sections 10 and 12 are additional analyses that are not reported in the article. Their code, results
and outcomes are kept in this repository.

## 1. Question and data
- Question: how much of the transcriptional response in single-cell chemical perturbation atlases
  comes from concentrations that patients do not reach.
- Tahoe-100M pseudobulk differential expression (DESeq2, compound x dose x cell line against
  vehicle), doses 0.05, 0.5 and 5 uM.
- sci-Plex 3 (scPerturb), A549, K562 and MCF7, 24 h, doses 0.01, 0.1, 1 and 10 uM.
- Clinical Cmax at the highest approved dose and plasma protein binding: Liston and Davis 2017
  (PMC5511563), matched to atlas compounds by generic name.

## 2. Definitions
- Unbound Cmax = Cmax x (1 - bound fraction). A range of binding is replaced by its mean,
  a lower bound such as "> 99%" by the bound, and the unbound fraction is not allowed below 0.005.
- DE gene: adjusted P < 0.05 and |log2FC| >= 1.
- Reachable dose: dose <= unbound Cmax (primary). Sensitivity: dose <= total Cmax.
- Serum-adjusted reference: unbound fraction in medium with 10% serum
  fu_med = 1 / (1 + 0.1 (1 - fu) / fu). Reachable if dose x fu_med <= unbound Cmax.
- Inert reference compounds (Tahoe-100M): talc, allantoin, folic acid, fumaric acid, arbutin,
  gallic acid, sodium salicylate, calcium lactate. Background of a plate and cell line: 95th
  percentile of the number of DE genes of these compounds.
- Background in sci-Plex 3: two vehicle wells against the remaining vehicle wells of the cell
  line, 95th percentile per cell line.
- Exposure groups: doses on both sides of the reference (at least one dose reachable and the
  top dose not), all doses above the reference, all doses reachable.

## 3. Share of the response above clinical exposure (Tahoe-100M)
Quantities
- Share of annotated compounds whose lowest dose exceeds unbound Cmax, and whose top dose
  exceeds 10 x unbound Cmax.
- Share of all DE calls of annotated compounds that come from doses above the reference.
- For compounds with doses on both sides of the reference: ratio of DE calls at the highest
  reachable dose to those at 5 uM (median over cell lines), and overlap of the gene sets.
Rule
- Supported: top-dose share >= 50%, share of DE calls >= 50%, median ratio <= 0.5 and the
  5 uM count above the background.
- Not supported: share of DE calls < 25% or median ratio > 0.8. Otherwise inconclusive.
Sensitivity
- |log2FC| thresholds 0.25, 0.5, 0.75, 1.0 and 1.5. Stable if the share of DE calls on the
  unbound reference is >= 50% at every threshold.
- Serum-adjusted and total Cmax references, reported.

## 4. Replication in sci-Plex 3
- Per-well pseudobulk, negative binomial model per cell line, same DE definition.
- Replicates if the share of DE calls from doses above unbound Cmax is >= 50% (strong: > 65%)
  and the median ratio of DE calls at the highest reachable dose to 10 uM is <= 0.5.
  Fails if the share is < 25% or the ratio > 0.8. Otherwise inconclusive.
- Cell number: where the 10 uM well has more cells than the 10 nM well of the same compound,
  cell line and replicate, the 10 uM well is reduced to that number of cells and the model is
  refitted. Passes if the reduced 10 uM contrasts keep >= 80% of their DE genes (median) and
  the share of DE calls stays > 65%.
- Threshold-free magnitude per contrast: (i) median |log2FC| over the 500 genes with the highest
  variance of log2FC across contrasts, (ii) L2 norm of beta^3 / (beta^2 + SE^2). A step at the
  clinical concentration is supported if, in compounds with doses on both sides of the reference,
  the median ratio of the lowest unreachable to the highest reachable dose is >= 2 for both
  measures and the vehicle-against-vehicle magnitude is below that of the highest reachable dose.

## 5. Selection of annotated compounds
- 5 uM DE counts and dose slope of annotated against other compounds (Mann-Whitney and
  Kolmogorov-Smirnov), also restricted to compounds of the mechanism classes present among
  the annotated ones.
- Extension set: Cmax and protein binding from product labels with source links, five rows
  checked against the label. Consistent if the share of DE calls in the added set is > 70%.

## 6. Specificity of the response at high doses (Tahoe-100M)
- For compounds with doses on both sides of the reference, DE genes at 5 uM are divided into
  those also DE at a reachable dose and those DE only at 5 uM. For compounds with all doses
  above the reference, 0.05 and 0.5 uM take the place of the reachable doses (descriptive).
- Convergence on a common response is supported if the Jaccard index between compounds of
  different mechanism classes (same cell line, all sets >= 10 genes) is higher for 5 uM-only
  genes than at reachable doses, with Wilcoxon P < 0.001 and a ratio of medians >= 1.5.
- Enrichment of Hallmark apoptosis, p53, unfolded protein response and reactive oxygen species
  genes among 5 uM-only genes against genes shared with reachable doses, Fisher test, reported.
- Mechanism similarity: classes with >= 2 compounds, profiles of the 1,000 genes with the
  highest variance, Spearman correlation between compounds per cell line and dose, AUROC of
  pairs of the same class against pairs of different classes. Loss of mechanism information at
  5 uM is supported if the AUROC at 5 uM is below that at 0.05 uM (Wilcoxon over cell lines,
  P < 0.05) and the similarity between classes is higher at 5 uM (P < 0.05). It is contradicted
  if the AUROC at 5 uM is at or above that at 0.05 uM with P < 0.05.

## 7. Agreement with patient tumours
- Cohorts with paired samples before and during single-agent treatment: GSE119262 everolimus,
  GSE192464 afatinib, GSE93204 palbociclib (day 15 against day 1 of cycle 1) and anastrozole,
  GSE153470 anastrozole, GSE58837 sunitinib, GSE341672 olaparib. Patient change = mean
  within-patient difference. The olaparib cohort (sequencing batch coincides with time point)
  and the second anastrozole comparison are reported separately.
- For all primary compounds the lowest atlas dose exceeds unbound Cmax, so that the comparison
  is between doses. Low-dose signature: DE at 0.05 uM. High-dose signature: DE at 5 uM and at
  neither lower dose. A gene enters a signature if DE with the same sign in >= 3 cell lines.
  In vitro change = median log2FC across cell lines.
- Measures: sign agreement with the patient change (binomial against 0.5), Fisher test between
  signatures, Spearman correlation of full profiles.
- Better agreement of the low-dose signature is supported if its agreement is higher in most
  cohorts and in the pooled comparison (P < 0.05), or if stress genes of the high-dose signature
  are not changed in patients (their share among patient-changed genes below half of that in
  the signature).
- Aggregated comparisons: (i) the four agents with a direct target in cancer cells, (ii) the
  50 and 100 signature genes with the largest in vitro change, (iii) scores of the 50 Hallmark
  gene sets. In vitro score = mean over set genes of the median log2FC (genes measured in >= 10
  cell lines) minus the mean over all genes, at 0.05 and 5 uM. Patient score = mean paired t
  statistic of the set minus the mean over all genes. Sets with >= 10 genes on both sides.
  Target sets: mTORC1 signalling (everolimus), E2F targets and G2M checkpoint (palbociclib),
  early and late oestrogen response (anastrozole), angiogenesis and hypoxia (sunitinib).
- Agreement of gene set scores at a dose: pooled sign agreement over primary cohorts and sets
  above 50% (binomial one-sided P < 0.05) and Spearman > 0 in most primary cohorts.
- Better agreement at 5 uM than at 0.05 uM: higher Spearman in most primary cohorts and higher
  pooled sign agreement with Fisher P < 0.05. Otherwise no difference is claimed.

## 8. Power at clinical exposure (sci-Plex 3) and exposure tiers
- Responsive genes of a compound and cell line: DE genes at 10 uM, at least 10 genes.
- Projected change at a lower dose: b = log2FC x sign(log2FC at 10 uM). Negative control: the
  same projection on contrasts of other compounds at reachable doses in the same cell line.
  Change = mean(b) - mean(control).
- Spread of true changes: tau^2 = var(b) - mean(SE^2), not below 0. 90th percentile = change +
  1.2816 tau. Groups by dose / unbound Cmax: <= 1, 1 to 10, 10 to 100. Confidence intervals by
  resampling compounds.
- Variance of the log count of a gene with per-cell rate m in a well of N cells:
  v(N) = (1 + a1) / (N m) + a0, with a0 and a1 from the fitted dispersion trend.
  SE^2 of the log fold change = v(N) / r + v(N0) / r0 for r treated and r0 vehicle wells.
- Significance level: alpha* = 0.05 x 0.8 x |R| / (genes tested).
- Reported: cells per well at r = 2 and wells at N = 150 that give a median power of 80% over
  responsive genes, for the 90th percentile and for the mean change at reachable doses, and the
  same for matching the mean test statistic at 10 uM.
- A ceiling is declared if the median power at r = 2 stays below 80% for unlimited cells.
- Exposure tiers: tier 1 dose <= unbound Cmax. Tier 2 not tier 1 and dose x fu_med <= unbound
  Cmax. Tier 3 neither and dose / unbound Cmax <= 10. Tier 4 neither and ratio > 10.
  Optional weight w = min(1, (unbound Cmax / dose)^2), not evaluated in model training.

## 10. What a response model learns from the atlas (Tahoe-100M)
Fixed on 17 November 2026, before the profiles were assembled and before any model was trained.
The section was shortened on the same day, also before any result, to the three rules below.
- Question: when a model is trained to predict the Tahoe-100M response, does its prediction skill
  and its training signal come from conditions above clinical exposure.
- Data: log2 fold changes of the released DE table for the genes of the 50 hallmark gene sets.
  Values with baseMean < 5 are missing. Values are clipped to [-5, 5]. Genes missing in more than
  20% of contrasts are dropped. Contrasts on several plates are averaged. Annotated compounds with
  DE results only, so that every condition has a tier.
- Split: 5 folds over (compound, cell line) pairs. The three doses of a pair are held out together.
- Model: compound, cell line and dose embeddings (512 units, summed), two hidden layers of 512
  units with ReLU, linear output for all genes, mean squared error over non-missing genes, Adam,
  learning rate 0.001, batch 256, 60 epochs. No tuning on test folds.
- Reference: context mean, the mean training profile of the same cell line and dose.
- Metric per held-out contrast: r_spec, the Pearson correlation across genes of
  (prediction - context mean) and (observed - context mean). This is the compound-specific part of
  the prediction. The plain correlation r and R2 against no change are reported as well.
- Summaries: mean per compound and tier, then mean over compounds. 95% intervals from 2,000
  resamples of compounds.
- Rules:
  - R1 training signal. Share of the sum of squared log2 fold changes (the loss of a model that
    predicts no change) held by tier 3 and 4 conditions against their share of conditions.
    Supported if the lower interval bound of the loss share exceeds the share of conditions.
  - R2 skill by exposure. Compound-specific prediction is present in a tier if the interval of
    r_spec lies above 0. Skill confined to supraphysiological exposure is supported if r_spec in
    tier 4 exceeds r_spec in tier 1 (interval of the difference above 0). Whether the tier 1
    interval includes 0 is reported.
  - R3 training ablation. Arms: all conditions, without tier 3 and 4 conditions, without tier 1
    and 2 conditions. Learning from clinical exposure is supported if removal of tier 1 and 2
    conditions lowers r_spec on held-out tier 1 and 2 conditions (interval of the difference above
    0). Dependence on supraphysiological conditions is supported if removal of tier 3 and 4
    conditions lowers r_spec on held-out tier 3 and 4 conditions.
- Either outcome of each rule is reported.

## 11. Sensitivity of the patient reference (cohorts of section 7)
Fixed on 17 November 2026, after the results of section 7 and before any result of this section.
- Question: whether gene-level agreement close to 50% in section 7 is explained by measurement
  noise of bulk tumour profiles (five microarray cohorts, one RNA sequencing cohort).
- (a) Agreement between cohorts. The two anastrozole comparisons (GSE153470, Illumina, and
  GSE93204, Agilent). Signature of a cohort = the 500 genes with the smallest paired t-test P
  among genes measured in both cohorts. Agreement = proportion with the same sign of the mean
  change in the other cohort, in both directions, with Wilson intervals.
- (b) Agreement within a cohort. Cohorts with >= 16 pairs. Patients are divided at random into
  two halves, 200 times. Signature = the 500 genes with the smallest P in one half. Agreement =
  proportion with the same sign of the mean change in the other half, averaged over divisions.
- (c) Atlas signatures restricted to genes that changed in patients. The agreement of section 7
  is recomputed for signature genes with a paired t-test P < 0.05 in the cohort, pooled over the
  primary cohorts, for the 0.05 uM and the 5 uM-only signature.
- Rules:
  - The patient reference is sensitive if the agreement in (a) exceeds 50% in both directions
    (lower interval bound above 50%) and the agreement in (b) exceeds 60% in most cohorts.
  - Noise as explanation of the section 7 result is supported if the pooled agreement in (c) lies
    above 50% for either signature (lower interval bound above 50%). It is not supported if the
    reference is sensitive and both intervals in (c) include 50% or lie below it.
  - If the reference is not sensitive, the gene-level comparison of section 7 is declared
    inconclusive and only the gene set result is retained.
- Either outcome is reported.

## 12. Agreement with malignant cells of patient tumours (single-cell cohorts)
Fixed on 17 November 2026, after inspection of sample annotation and cell numbers only and before
any expression value of these cohorts was read.
- Question: whether the gene-level result of section 7 is explained by the non-malignant cells of
  bulk tissue. The comparison is repeated in malignant cells only.
- Cohorts. GSE158724 (FELINE trial, breast cancer, single-nucleus RNA sequencing, 10x samples):
  biopsies at day 0 and day 14 of letrozole (arm A) or letrozole with ribociclib (arms B and C),
  cancer cells as annotated by the authors. GSE292469 (lung adenocarcinoma, single-nucleus RNA
  sequencing): tumour cells as flagged by the authors, treatment-naive samples and samples of
  residual disease under osimertinib. Samples at progression are not used. The decitabine series
  GSE279925 holds promoter methylation only and is not used.
- Profiles: counts summed over the malignant cells of a sample, samples with >= 30 malignant cells,
  counts per million, log2(x + 1). Genes with a mean of >= 5 counts per million.
- Patient change. Ribociclib (primary): mean within-patient change from day 0 to day 14 in arms
  B and C minus that in arm A, Welch t-test between arms. Secondary: within-patient change in
  arms B and C, paired t-test. Osimertinib (secondary, between patients): mean of residual disease
  samples minus mean of treatment-naive samples, Welch t-test.
- Atlas signatures and measures as in section 7 (0.05 uM signature, 5 uM-only signature, sign
  agreement with Wilson interval, Fisher test between signatures, Spearman correlation of full
  profiles), and as in section 11 (c) for signature genes with patient P < 0.05. The 0.5 uM
  signature is added for ribociclib, for which 0.05 and 0.5 uM are reachable doses.
- Sensitivity of the reference: change of the E2F targets and G2M checkpoint gene set scores under
  ribociclib (expected negative), and agreement between random halves of patients as in
  section 11 (b) for the secondary ribociclib comparison.
- Rules:
  - Agreement in malignant cells is present for a signature if the lower interval bound of its
    sign agreement in the primary ribociclib comparison lies above 50%.
  - Non-malignant cells as explanation of the section 7 result is supported if agreement is
    present for at least one signature in the primary ribociclib comparison. It is not supported
    if the intervals of all signatures include 50% or lie below it.
  - A difference between doses is claimed only with Fisher P < 0.05.
  - The osimertinib comparison is reported separately and does not enter the rules.
- Either outcome is reported.

## 9. Outcome of the rules
| Rule | Outcome |
|---|---|
| Share of the response above unbound Cmax, Tahoe-100M (section 3) | Supported: top-dose share 74%, share of DE calls 78%, median ratio 0.00 |
| Same with total Cmax | Inconclusive: share of DE calls 55%, top-dose share 40% |
| Stability across thresholds | Stable: 76.1% to 77.7% |
| Replication in sci-Plex 3 (section 4) | Replicates: share of DE calls 91.1%, median ratio 0.006 |
| Cell number | Passes: 96% of DE genes kept, share 91.0% |
| Step of magnitude at the clinical concentration | Not supported: ratio 1.05 and 1.09 |
| Label extension set (section 5) | Consistent: 71.3% |
| Convergence on a common response at 5 uM (section 6) | Not supported: Jaccard 0.010 for 5 uM-only genes against 0.024 |
| Loss of mechanism information at 5 uM | Neither supported nor contradicted: AUROC 0.48 at both doses |
| Better agreement of the low-dose signature with patients (section 7) | Not supported: 48.0% against 48.6% |
| Agreement of gene set scores | Met at 0.05 uM (56.0%), not met at 5 uM (46.8%) |
| Better agreement of gene set scores at 5 uM | Not met |
| Ceiling for the mean change at reachable doses (section 8) | Declared: power 24% with two wells and unlimited cells |
| Training signal held by tier 3 and 4 conditions (section 10, R1) | Not supported: 61.9% of squared change (54.5 to 69.1) against 56.7% of conditions |
| Prediction confined to supraphysiological exposure (R2) | Not supported: r_spec 0.52 in tier 1 and 0.55 in tier 4, difference 0.024 (-0.001 to 0.048) |
| Learning from clinical exposure (R3) | Supported: removal of tier 1 and 2 lowers r_spec on tier 1 and 2 by 0.31 (0.28 to 0.33) |
| Dependence on supraphysiological conditions (R3) | Supported for tier 3 and 4 only: lowered by 0.37 (0.35 to 0.40), effect on tier 1 and 2 0.02 |
| Sensitivity of the patient reference (section 11) | Sensitive: agreement between the anastrozole cohorts 94.6% and 93.0%, within-cohort agreement above 60% in 5 of 6 cohorts (palbociclib 59.2%) |
| Noise as explanation of the gene-level result | Supported for the 0.05 uM signature only: among genes changed in patients 56.0% (50.1 to 61.8), 5 uM-only signature 48.6% (42.1 to 55.2), Fisher P 0.12 |
| Agreement in malignant cells, ribociclib (section 12) | Not present: 36.4% (23.8 to 51.1) at 0.05 uM, 56.8% (45.4 to 67.4) at 0.5 uM, 58.8% (36.0 to 78.4) for 5 uM-only, 44, 74 and 17 genes |
| Non-malignant cells as explanation of the section 7 result | Not supported. The reference check for the primary comparison was not met (E2F score difference between arms 0.09, within arms -1.71 and -1.30) |
| Osimertinib, residual disease against untreated (reported separately) | 38.7% (34.0 to 43.6) at 0.05 uM and 39.2% (35.8 to 42.7) for 5 uM-only, Fisher P 0.90 |

Exploratory analyses added after the results of section 10: subtraction of the mean profile of the
other compounds on the same plate, cell line and dose (r_spec 0.38 in tier 1, 0.44 in tier 4), and
addition of the eight inert compounds (r 0.62, plate-adjusted 0.43 against 0.44 in tier 4).

Corrections in section 12 after the first run, rules unchanged: genes without variance of the
within-patient change are dropped (two genes), and the time point of each osimertinib sample is
taken from the cell table (six residual disease samples had been grouped with untreated samples).
