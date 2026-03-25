# GROWell WGBS Comethylation Analysis


 Bioinformatics analysis for identifying co-methylation networks associated with maternal metabolic outcomes in the GROWell clinical trial, focusing on postpartum weight retention (PPWR) and adverse pregnancy outcomes (APOs).

 **Project Overview**

This project applies a co-methylation network framework (comethyl) to whole-genome bisulfite sequencing (WGBS) data derived from participant-collected dried blood spots (DBS) in the GROWell (Goals for Reaching Optimal Wellness) study.

The goal is to:

Identify modules of correlated DNA methylation regions
Associate these modules with:
Metabolic outcomes (PPWR)
Clinical outcomes (APOs)
Behavioral and environmental covariates

# Data Description
**Data type(s):** Whole-Genome Bisulfite Sequencing (WGBS)
**Study/cohort:** GROWell (Goals for Reaching Optimal Wellness)
**Sample type:** Dried Blood Spots (DBS)
**Timepoints:**
 1. Baseline (early pregnancy)
 2. 36–38 weeks gestation
 3. ~3 months postpartum
**Genome build:** hg38
---

## Quick start

### Clone
```bash
git clone https://github.com/dreusebio/wgbs_growell_comethylation_analysis.git
cd [REPO]
```

---

## Repository structure

- `data/` - raw + processed data + codebooks, data dictionaries, provenance notes (**not commited**)
- `scripts/` — entrypoints + SLURM submit scripts
- `analysis/` — downstream statistics/figures +  configuration + sample sheets + generated logs
- `docs/` — methods + workflow documentation
- `results/` — generated outputs
- `test/` —  Test data and vignettes

---

#Pipeline Overview (comethyl)

The workflow follows a structured, reproducible pipeline:

1. Preprocessing
Input: Bismark cytosine reports (CpG_report.txt.gz)
Coverage filtering
CpG clustering (≥3 CpGs within ≤150 bp)
2. Region Filtering
Minimum coverage threshold
Variability filtering (e.g., SD > 0.05–0.08)
3. Methylation Matrix Construction
Region-level methylation aggregation
Output: Region_Methylation.rds
4. Adjustment for Confounders

Three analysis strategies:

v1_all_pcs → adjust for all PCs
v2_exclude_outcome_exposure_pcs → remove biologically relevant PCs
v3_technical_pcs_only → adjust only technical variation
5. Network Construction
Soft-threshold selection (getSoftPower)
Adjacency → TOM → clustering
Module detection (WGCNA-based)
6. Module–Trait Association
Correlation:
bicor (robust, preferred)
pearson (optional)
Outputs:
Heatmaps
Dot plots
Scatter plots
7. Functional Annotation
GREAT (preferred)
Offline annotation (EnsDb / TxDb fallback)
Enrichment:
GO
KEGG
ClinVar / GWAS (optional)

---

## Reproducibility standard

See `docs/lab_reproducibility_standard.md`.

---

## License

MIT License (see `LICENSE`)
