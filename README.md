# cohortguard

**Audit, before validation labels exist, whether a transcriptomic signature will transfer to a new cohort.**

A signature is usually judged twice: by a cross-validated AUC in the cohort that produced it, and later by an external AUC if a second labelled cohort appears. Both numbers are retrospective. `cohortguard` measures the quantities that are available *before* the labels arrive — drift and model agreement — and reports how far each falls short of predicting the drop.

This package accompanies the manuscript *"Cross-cohort transferability of transcriptomic signatures is poorly predicted by within-cohort performance and label-free drift scores under a fixed microarray pipeline."* It lets anyone recompute every number in that paper from public GEO series on a common scale.

## What it computes

| Quantity | Meaning |
|---|---|
| Kernel MMD drift | Distance between the two expression matrices (RBF kernel, median-heuristic bandwidth, optionally restricted to the top-N most variable genes) |
| Model-agreement drift | How far independently fitted models stop agreeing on the target cohort |
| Cross-validated AUC | Within-cohort performance (repeated stratified k-fold) |
| External AUC and drop | Performance on the target cohort, and the difference from cross-validation |
| Split-conformal coverage and abstention | Per-sample prediction sets and the nominal-coverage guarantee under exchangeability |

## Requirements

R >= 4.0 and base R only. Figure generation additionally uses `ggplot2`.

```bash
R CMD INSTALL .
```

```r
library(cohortguard)
```

## Getting the data

The raw expression files are **not** committed to this repository (they total several gigabytes and are all public). Download the 25 GEO series matrices and their platform annotations from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) into this layout, which `data/datasets.csv` already points at:

```
data/raw/GSE96804/GSE96804_series_matrix.txt.gz     # plus its GPL*.txt annotation
data/raw/GSE30122/GSE30122_series_matrix.txt.gz
...                                                  # one folder per GSE in datasets.csv
data/raw/_platform/GPL570-55999.txt                  # annotations shared by several series
```

The accession numbers and platforms are listed in `data/datasets.csv`.

## Running the analysis

Regenerate the whole 16-transfer benchmark and the figures:

```bash
Rscript scripts_r/run_benchmark.R      # -> data/processed/benchmark_R.csv
Rscript scripts_r/make_figures.R       # -> figures/*.pdf
```

The revision and sensitivity analyses behind the paper are in `scripts_r/` as well:

```bash
Rscript scripts_r/revision_metrics.R   # exclusion, bandwidth, random-gene controls
Rscript scripts_r/m12_sensitivity.R    # training-statistics standardisation
Rscript scripts_r/m8_m11.R             # continuous agreement + DeLong intervals
Rscript scripts_r/m14_agreement_se.R   # bootstrap SE of the agreement drift
```

## Using the functions directly

```r
library(cohortguard)

# Az, yA : z-scored training matrix (genes x samples) and 0/1 labels
# Bz, yB : z-scored target matrix (same genes) and labels
res <- cg_evaluate(Az, yA, Bz, yB, k = 50, seed = 7)
print(res)   # cv, ext, drop, mmd, mmd_norm, agree_delta, cov90, abstain90
```

Lower-level pieces: `cg_read_series_matrix`, `cg_read_gpl`, `cg_to_gene`, `cg_zscore`,
`cg_auc`, `cg_tstat`, `cg_fit_logit`, `cg_predict`, `cg_mmd`, `cg_mmd_topvar`,
`cg_agreement`, `cg_cv_auc`, `cg_split_ref`, `cg_conformal`, `cg_evaluate`.

## Repository layout

```
R/cohortguard.R          core functions (readers, AUC, ridge logistic, MMD, agreement, conformal)
scripts_r/               driver scripts: benchmark, figures, revision and sensitivity analyses
data/datasets.csv        transfer configuration for all 16 transfers
tests/test_r.R           unit tests (AUC ties, conformal coverage, drift identities)
DESCRIPTION, NAMESPACE   package metadata and exports
```

## What the paper finds

Across 16 cross-cohort transfers under one fixed pipeline, the drop from cross-validated to external AUC ranged from -0.200 to +0.692, and none of the label-free quantities tested predicted it usefully: MMD without labels explained 11% of its variance, model-agreement drift 2%, and external AUC was essentially uncorrelated with cross-validated AUC (R2 = 0.002). Split-conformal coverage held within cohorts but collapsed across them, and abstention tracked in-cohort difficulty rather than shift.

## Scope and limitations

- Drift is characterised by an RBF-kernel MMD and by model-agreement decline; the latter is insensitive to scale differences.
- Coverage guarantees require exchangeability between calibration and target; across cohorts they do not hold, which is the reason for cohort-level judgement.
- Everything here is restricted to microarray data and one classifier family (ridge logistic regression on 50 genes); see the manuscript's Limitations section.

## Citation

Manuscript under review at *BMC Bioinformatics*. If you use this code, please cite the archived release:

> Li Y. cohortguard: label-free audit of transcriptomic signature transferability. Zenodo; 2026. doi:10.5281/zenodo.24284477

## License

MIT. See `LICENSE`.
