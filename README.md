# Pain Burden, Plasma Proteomics, and Incident Parkinson Disease
## Data Availability and Access

This repository contains analysis code only. Individual-level participant data
and proteomic measurements from UK Biobank, the China Health and Retirement
Longitudinal Study (CHARLS), the Parkinson's Progression Markers Initiative
(PPMI), and the Olink-based proteomic datasets used in this study are not
included or redistributed.

Access to these data is controlled by the respective data providers and is
subject to their application procedures, data-use agreements, and participant
confidentiality requirements. Researchers wishing to reproduce the analyses
must apply independently to the relevant data custodians and obtain the
necessary approvals. The authors of this repository cannot provide, transfer,
or authorize access to these datasets.

After obtaining authorized access, users may place the required local input
files under `code/data/` or set the environment variable `PD_PAIN_DATA_DIR` to
another secure, private data directory. The filenames listed below refer to
locally prepared analysis files and are not distributed with this repository.

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
