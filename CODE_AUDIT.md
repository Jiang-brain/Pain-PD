# Code Reorganization and Audit Notes

## Files Created

- `00_setup/config.R`
- `00_setup/ukb_helpers.R`
- `00_setup/model_helpers.R`
- `01_data_preparation/01_build_ukb_followup_outcome.R`
- `01_data_preparation/02_define_ukb_pain_exposures.R`
- `01_data_preparation/03_define_ukb_covariates.R`
- `01_data_preparation/04_prepare_ukb_analysis_dataset.R`
- `01_data_preparation/05_baseline_characteristics.R`
- `02_epidemiology/01_site_specific_pain_cox.R`
- `02_epidemiology/02_pain_status_multisite_cox.R`
- `02_epidemiology/03_interactions_and_sensitivity.R`
- `02_epidemiology/04_composite_pain_score_and_rcs.R`
- `06_figures/01_plot_main_epidemiology_results.R`
- `05_prediction/04_predict_pd_stable_elasticnet374_lgbm_cv.py`
- `README.md`
- `requirements.txt`
- `.gitignore`

## Files Moved or Renamed

- `PD_pain.R` -> `archive/PD_pain_legacy.R`
- `analysis_pain_score_protein_association.R` -> `03_proteomics/01_pain_score_protein_association.R`
- `analysis_protein_pd_cox.R` -> `03_proteomics/02_protein_pd_cox.R`
- `compare_ppmi_ukb_protein_pd.py` -> `03_proteomics/03_compare_ppmi_ukb_protein_pd.py`
- `analysis_overlap25_iorw_joint_mediation.R` -> `04_mediation/01_overlap25_iorw_joint_mediation.R`
- `analysis_overlap25_iorw_per_protein_mediation.R` -> `04_mediation/02_overlap25_iorw_per_protein_mediation.R`
- `analysis_overlap25_pca_prediction_cindex.R` -> `05_prediction/01_overlap25_prediction_cindex_nri.R`
- `build_stable_nested_elastic_lasso_pfs.py` -> `05_prediction/02_build_stable_nested_elastic_net_pfs.py`
- `analysis_stable_nested_pfs_pd_cindex.R` -> `05_prediction/03_stable_elastic_net_pfs_pd_cindex_nri.R`
- added/standardized `05_prediction/04_predict_pd_stable_elasticnet374_lgbm_cv.py`

## Important Fixes

- Removed PCA models from the overlap-25 prediction script. The organized final
  analysis uses the 25 shared proteins directly, not a PCA signature.
- Confirmed that the PPS model is Elastic Net because it uses `ElasticNetCV`
  with `l1_ratio` values between 0 and 1. This includes both L1 and L2
  regularization.
- Removed the old dependency on nested LASSO coefficient files from the PPS
  builder.
- Replaced the missing `analysis_stable_nested_pfs_pd_cindex_base.R` dependency
  with a self-contained Elastic Net PPS Cox/NRI script.
- Replaced the LightGBM wrapper dependency on missing
  `predict_pd_intersection_lgbm_cv.py` with a self-contained stable Elastic Net
  374-protein LightGBM script.
- Replaced personal absolute input/output paths in the main organized scripts
  with `code/data/` and `code/results/`.
- Moved the original large script to `archive/` rather than deleting it.

## Remaining Code Issues

- `archive/PD_pain_legacy.R` still contains personal absolute paths and
  exploratory MR code. It is retained only for provenance.
- The IORW mediation scripts still contain duplicated helper definitions from
  the original implementation. Their inputs/outputs were standardized, but the
  core algorithm was intentionally left unchanged.
- R scripts could not be executed or parsed in this environment because
  `Rscript` was not available.
- Full data-level checks could not be run because raw UK Biobank, CHARLS, Olink,
  and PPMI data are not included.

## Manuscript Analyses Without Complete Final Code in This Folder

- CHARLS external validation logistic regression.
- STRING protein-protein interaction and GO/KEGG/Reactome/WikiPathways
  enrichment analyses.
- Formal Schoenfeld residual proportional hazards assumption checks.
- Fully cleaned sensitivity analyses for dementia exclusion, pain medication,
  vascular disease, cancer, diabetes, and sleep duration.

## Leakage and Modeling Audit

- Elastic Net median imputation and standardization are inside the sklearn
  pipeline fitted on each training fold.
- Nested Elastic Net uses fixed seeds for outer and inner K-fold splits.
- The organized Elastic Net PPS builder no longer reads LASSO-selected features.
- C-index/NRI scripts use fitted Cox risks; higher predicted absolute risk is
  treated as higher event risk.
- The PPS Cox/NRI evaluation remains an internal cohort evaluation. External
  validation is not implemented in the current code.
