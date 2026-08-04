# Pain Burden, Plasma Proteomics, and Incident Parkinson Disease

This folder contains the analysis code supporting the manuscript `Main text.docx`.
The code has been reorganized to follow the final analysis workflow in the Methods
and Results sections.

Raw participant-level UK Biobank, CHARLS, Olink, and PPMI data are not included in
this repository. Place local input files under `code/data/` or set
`PD_PAIN_DATA_DIR` to another private data directory.

## Main Analysis Methods

- UK Biobank cohort construction, incident PD follow-up time, and 2-year landmark
  exclusion to reduce reverse causation.
- Pain exposure definition by seven body sites, duration, single-site/multisite
  burden, and weighted composite pain score.
- Cox proportional hazards models for site-specific pain, pain status, cumulative
  pain burden, subgroup interaction tests, and sensitivity analyses where source
  variables are available.
- Plasma proteome-wide linear models for composite pain score and Cox models for
  incident PD, with FDR correction.
- Directionally concordant overlap protein selection, IORW mediation analyses,
  and prediction analyses using the 25 shared proteins.
- Pain-related proteomic signature (PPS) development using nested cross-validated
  Elastic Net. The final PPS model is Elastic Net, not pure LASSO.
- Final organized analyses do not use PCA. The former PCA-based prediction branch
  was removed from the main prediction script.

## Directory Structure

```text
code/
  00_setup/
    config.R
    ukb_helpers.R
    model_helpers.R
  01_data_preparation/
    01_build_ukb_followup_outcome.R
    02_define_ukb_pain_exposures.R
    03_define_ukb_covariates.R
    04_prepare_ukb_analysis_dataset.R
    05_baseline_characteristics.R
  02_epidemiology/
    01_site_specific_pain_cox.R
    02_pain_status_multisite_cox.R
    03_interactions_and_sensitivity.R
    04_composite_pain_score_and_rcs.R
  03_proteomics/
    01_pain_score_protein_association.R
    02_protein_pd_cox.R
    03_compare_ppmi_ukb_protein_pd.py
  04_mediation/
    01_overlap25_iorw_joint_mediation.R
    02_overlap25_iorw_per_protein_mediation.R
  05_prediction/
    01_overlap25_prediction_cindex_nri.R
    02_build_stable_nested_elastic_net_pfs.py
    03_stable_elastic_net_pfs_pd_cindex_nri.R
    04_predict_pd_stable_elasticnet374_lgbm_cv.py
  06_figures/
    01_plot_main_epidemiology_results.R
  archive/
    PD_pain_legacy.R
  requirements.txt
```

## Input Files

By default, put private input files in `code/data/`:

- `PD_data.csv`: UK Biobank analytical source file.
- `protein_new_missing60_filtered.csv`: Olink NPX protein matrix after QC.
- `PD_data_matched_protein_covariates.csv`: proteomic technical covariates.
- `PD_data_matched_protein_times.csv`: blood-to-assay timing variables.
- `protein_intersection_FDR05.csv`: FDR-significant pain-PD overlap results.
- `results_PPMI.csv`: PPMI protein-PD comparison data.

Generated intermediate and result files are written to `code/results/`. Figures
are written to `code/results/figures/`.

## Recommended Run Order

Run commands from the `code/` directory.

```bash
Rscript 01_data_preparation/01_build_ukb_followup_outcome.R
Rscript 01_data_preparation/02_define_ukb_pain_exposures.R
Rscript 01_data_preparation/03_define_ukb_covariates.R
Rscript 01_data_preparation/04_prepare_ukb_analysis_dataset.R
Rscript 01_data_preparation/05_baseline_characteristics.R

Rscript 02_epidemiology/01_site_specific_pain_cox.R
Rscript 02_epidemiology/02_pain_status_multisite_cox.R
Rscript 02_epidemiology/03_interactions_and_sensitivity.R
Rscript 02_epidemiology/04_composite_pain_score_and_rcs.R

Rscript 03_proteomics/01_pain_score_protein_association.R
Rscript 03_proteomics/02_protein_pd_cox.R
python 03_proteomics/03_compare_ppmi_ukb_protein_pd.py

Rscript 04_mediation/01_overlap25_iorw_joint_mediation.R
Rscript 04_mediation/02_overlap25_iorw_per_protein_mediation.R

Rscript 05_prediction/01_overlap25_prediction_cindex_nri.R
python 05_prediction/02_build_stable_nested_elastic_net_pfs.py
Rscript 05_prediction/03_stable_elastic_net_pfs_pd_cindex_nri.R
python 05_prediction/04_predict_pd_stable_elasticnet374_lgbm_cv.py

Rscript 06_figures/01_plot_main_epidemiology_results.R
```

## Environment

R packages:

- `data.table`
- `survival`
- `rms`
- `ggplot2`
- `jsonlite`

Additional packages used only by the archived legacy script include
`survminer`, `Gmisc`, `Greg`, `tidyverse`, `forestplot`, `ggthemes`,
`ggcorrplot`, `RColorBrewer`, `interactionR`, and `TwoSampleMR`.

Python packages are listed in `requirements.txt`:

- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `lightgbm`

## Reproducibility Notes

- Elastic Net preprocessing is inside a sklearn `Pipeline`, so median imputation
  and standardization are fitted within each outer training fold.
- Nested Elastic Net uses fixed random seeds: `RANDOM_STATE = 20260525`,
  10 outer folds, and 10 inner folds.
- The final prediction scripts use Cox model linear predictors in their fitted
  risk direction. Continuous NRI is computed from Cox absolute risk at the
  evaluation horizon.
- The organized main code has no hard-coded personal absolute paths. Personal
  absolute paths remain only in `archive/PD_pain_legacy.R` for provenance.

## Analyses That Still Need Manual Confirmation

The manuscript Methods/Results mention several analyses that were not available
as complete, clean, final scripts in the current `code/` folder:

- CHARLS external validation logistic regression code.
- STRING protein-protein interaction and GO/KEGG/Reactome/WikiPathways
  enrichment scripts.
- Formal Schoenfeld residual proportional hazards assumption checks.
- Fully cleaned sensitivity scripts for dementia exclusion, pain medication,
  vascular disease, cancer, diabetes, and sleep duration.

These gaps should be resolved before public GitHub release if the repository is
expected to reproduce every figure and supplementary table in the manuscript.
