# DXA-aging-pipeline

Code for a multi-omics study of **bone-imaging-derived biological age acceleration**: a deep
learning model estimates skeletal age from knee/hip radiographs, the gap between estimated and
chronological age ("age acceleration") is used as an exposure, and its associations with
bone mineral density, lifestyle, disease/survival outcomes, common-variant genetics and
single-cell bone tissue atlases are dissected.

> **This repository contains code only.** No individual-level data, no model weights, no
> summary statistics and no derived result tables are distributed here. All input paths appear
> as `/path/to/...` placeholders — see [Quick start](#quick-start).

---

## Contents

| Stage | Directory | What it does |
|---|---|---|
| 1 | [`01_imaging_model/`](01_imaging_model) | Build the analysis table, 5-fold cross-validated ResNet101 training, external validation |
| 2 | [`02_epi_outcomes/`](02_epi_outcomes) | Age calibration, age acceleration, lifestyle RCS, disease/survival, SMR/TWAS, proteomics, scDRS/scPagwas |
| 3 | [`03_gwas/`](03_gwas) | Genotype QC, sparse GRM, fastGWA association |
| 4 | [`04_downstream/`](04_downstream) | LDSC, MiXeR, PRS-CS, classic PRS, TWAS, LAVA, N-GWAMA, PolyFun fine-mapping |
| 5 | [`05_singlecell/`](05_singlecell) | cell2location spatial deconvolution of bone tissue |

```
radiographs ──► ResNet101 ──► estimated age ──┐
                                              ├─► age acceleration ──┐
chronological age ────────────────────────────┘                      │
                                                                     ├─► DXA BMD
   lifestyle · covariates · disease & death records ─────────────────┤   survival
                                                                     │   proteomics
   GWAS ──► LDSC ──► MiXeR / LAVA / N-GWAMA                           │
        └─► PRS-CS / classic PRS ──► PRS ────────────────────────────┤
        └─► PolyFun fine-mapping ──► credible sets                   │
        └─► TWAS / SMR                                               │
                                                                     │
   scRNA-seq + spatial transcriptomics ──► cell2location ──► cell types
```

---

## Pipeline in detail

### 1. Imaging model — `01_imaging_model/`

| Script | Role |
|---|---|
| `build_metadata.py` | Merge image inventory, phenotype, covariates and outcomes into one modelling table; stratified fold assignment |
| `train_cv_final.py` | 5-fold CV training of ResNet101 for 14 site combinations (left/right knee, left/right hip, single & double), out-of-fold predictions |
| `external_validation.py` | Cross-site (left↔right) and cross-cohort validation, MAE / RMSE / *r* |

Output: per-subject estimated age for each site combination.

### 2. Epidemiology & outcomes — `02_epi_outcomes/`

`F1` is the hinge between stages: it linearly calibrates the predicted age
(`lm(predicted_age ~ age)`) and defines

```
age_gap_adjusted = adjusted_age − chronological age
```

then merges BMD, covariates and outcomes, and writes the six key tables consumed downstream.

| Script | Role |
|---|---|
| `F1_calibration_and_phenotype.R` | Calibration, age acceleration, outcome merge |
| `F2_lifestyle_rcs.R` | Lifestyle exposures vs age acceleration (restricted cubic splines) |
| `F3_genetic_analysis.R` | Druggable targets, GWAS Catalog, FUMA, classic PRS, LDpred2 |
| `F4_singlecell_scdrs_scpagwas.R` | scDRS and scPagwas cell-type association |
| `F5_smr_twas.R` | SMR and TWAS integration |
| `F6_proteomics_enrichment.R` | Protein/metabolite association and enrichment |
| `scpagwas_single_site.R`, `scpagwas_combined_site.R` | scPagwas for single-site and combined phenotypes |
| `plot_density.ipynb` | Age-acceleration density plots |

### 3. GWAS — `03_gwas/`

`01`–`05` are plain PLINK2/GCTA shell steps run per chromosome:

```
01 keep samples → 02 missingness/MAF/HWE/INFO → 03 HapMap3 prune
→ 04 sparse GRM (180 parts) → 05 fastGWA-mlm
```

`06_gwas_covariate_model.R` prepares phenotypes/covariates and checks that the analysed sample
matches the genotype sample.

### 4. Downstream genetics — `04_downstream/`

| Module | Method | Key point |
|---|---|---|
| `ldsc/` | LD-SNP heritability and genetic correlation | `03_build_all_rg.sh` assembles `all.rg` from `*_rg.log` |
| `mixer/` | Univariate & bivariate causal variant sharing | 20 replicates, dockerised `gsa-mixer` |
| `prs_cs/` | PRS-CS with continuous shrinkage priors | 8 steps incl. duplicate-rsID harmonisation between weights and genotype |
| `prs_classic/` | LD clumping + thresholded PRS | `ieugwasr::ld_clump` |
| `twas/` | S-PrediXcan across GTEx tissues | |
| `lava/` | Local genetic correlation | `01_run_pair.R` is the reusable single-pair template |
| `ngwama/` | N-weighted / model-averaging multi-trait GWAMA | LDSC cross-trait intercept matrix corrects for sample overlap |
| `polyfun/` | Functionally-informed fine-mapping (PolyFun Approach 3) | 2 748 × 3 Mb windows, SuSiE with local PLINK LD |

PolyFun is the longest branch: baseline-LF annotations → S-LDSC binned regression → per-SNP
prior variance (SNPVAR) → regional SuSiE → PIP and credible sets → loci.

> **Three rules that matter in `polyfun/`:** always pass `--geno` (local LD, not precomputed
> LD scores), always pass `--allow-swapped-indel-alleles`, never pass `--allow-missing`.
> `03_prepare_genotype_cm.sh` must replace the `.bim` cM column with a real genetic map via
> `plink --cm-map`; with cM = 0 the LD matrix is wrong and PIPs collapse towards zero.

### 5. Single cell — `05_singlecell/`

```
01 gene filtering → 02 reference signature (RegressionModel)
→ 03 spatial mapping (Cell2location) → 04 extract → 05 QC → 06 abundance
```

---

## Quick start

Every script contains `/path/to/...` placeholders. Two ways to instantiate them:

**Option A — automatic**

```bash
bash scripts/instantiate_paths.sh \
    --project /data/my_project \
    --epi     /data/my_project_epi \
    --ref     /data/reference \
    --tools   /opt/bioinformatics
```

This rewrites the placeholders in place (a `.bak` is kept the first time).

**Option B — manual**

```bash
cp config/paths.example.sh config/paths.sh   # edit, then `source config/paths.sh`
cp config/paths.example.R  config/paths.R
```

| Placeholder | Meaning |
|---|---|
| `/path/to/project` | Python-side project root (metadata, GWAS, downstream results) |
| `/path/to/epi` | R-side project root (F1–F6, `data/`, `result/`) |
| `/path/to/ref` | Shared reference data: 1000G EUR, HapMap3, PRS-CS LD blocks, LDSC weights, GTEx models |
| `/path/to/tools` | Third-party executables and package sources |
| `/path/to/shared_project` | Shared upstream tables (IDs, proteomics matrices, scPagwas helpers) |
| `/path/to/external_project` | Code adapted from an earlier project |
| `/path/to/raw_genotype` | Raw imputed genotype input for QC |

All scripts assume the **project root as working directory**; `GWAS/...`, `data/...` are
relative to it. The directory structure each stage expects is described in
[`docs/data_layout.md`](docs/data_layout.md).

---

## Requirements

**R ≥ 4.1** — `data.table`, `dplyr`, `ggplot2`, `survival`, `rms`, `metafor`, `AICcmodavg`,
`LAVA`, `CMplot`, `foreach`, `doParallel`, `ieugwasr`, `genetics.binaRies`, `Seurat`,
`scPagwas`, `clusterProfiler`.

**Python ≥ 3.9** — `torch`, `torchvision`, `numpy`, `pandas`, `scikit-learn`, `scanpy`,
`cell2location`, `scvi-tools`, `matplotlib`, `tqdm`.

**Command-line** — PLINK 1.9 / PLINK 2.0, GCTA (`gcta64`), LDSC, PRS-CS, PolyFun, MiXeR
(docker or singularity image `ghcr.io/precimed/gsa-mixer`), S-PrediXcan, SMR.



## License

MIT — see [LICENSE](LICENSE). Third-party tools keep their own licenses.
