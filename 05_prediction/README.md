# Final Prediction Using Non-Proteomic-Derived Pain Scores

Release status: this originally controlled sensitivity branch now supplies the
adopted Figure 5 results (370 stable proteins). Its historical directory name
and original/new comparisons are retained. Use
`06_figures/03_plot_final_figure5.py` as the public plotting entry point. The
older standalone 25/374 workflow is not included in the submission tree.

This branch uses the single fixed-weight score from the independent pain-score workflow, not the mean of
the 1,000 bootstrap scores from the score-stability workflow. All original results remain intact.
The full original prediction scripts and saved results were located in PA5;
these, rather than the previously reorganized but unverified prediction scripts,
are the reference for controlled replication.

## Analyses

1. Reconstruct original Elastic Net stable scores and outer-fold predictions
   from saved coefficients. Reproduce original stable-panel Cox C-index/NRI and
   shared-25 Cox models before comparing new results. Independently rerun the
   original 374-protein LightGBM models and compare AUC/AUPRC/Brier with saved
   metrics. Reproduction failures are reported, not forced to match.
2. Match the new pain scores to the same 48,720 Elastic Net participants, keeping
   their order, 2,920 input proteins, outer/inner folds, seeds, target transform
   `log(pain_score + 0.1)`, alpha grid, and penalties fixed.
3. Refit all 10 outer folds, each with the full original 10-fold inner
   ElasticNetCV, 70 alphas, and L1 ratios 0.5/0.75/0.9. This is Elastic Net with
   both L1 and L2 penalties. Protein selection is rerun from the full original
   protein universe, not restricted to the original 374 or 834 marginal hits.
4. Select proteins retained in all 10 folds. Do not force the new count to 374.
   Build the original mean-coefficient stable signature using the new folds.
   Report genuine outer-fold OOF pain prediction separately from the stable
   mean-coefficient score's apparent fit in the derivation population.
5. Refit the original Cox prediction models for the new stable panel and its
   single proteomic pain score, separately. The manuscript's original C-index
   0.831/0.888 refers to the MULTIVARIABLE protein panel, not the single PPS.
6. Compare the 24 shared proteins with the original 25 on the same original
   complete-case cohorts. The only shared-set change is dropping COMP; original
   participant IDs, covariates, and protein values are held fixed.
7. Run original LightGBM classification models with the new stable protein panel
   and new shared-24 panel, using the original 20-PC covariate specification,
   class weighting, native missing/categorical handling, model settings, and
   stratified 10-fold splits. The scoring cohort remains the original LightGBM
   cohort, not the more restrictive pain-score training cohort.

## Important Interpretive Limits

- VERIFIED ORIGINAL R COX BUG: sequential recoding collapses deprivation and PRS
  to a single observed level; `usable_covariates` then excludes them. The original
  R prediction models do not actually adjust for these two variables. This
  one-factor sensitivity preserves the original behavior, and records it in
  `legacy_Cox_covariate_variation_audit.csv` and `methodological_warnings.csv`.
  Corrected-covariate models require a separately labelled analysis.
- Original Elastic Net preprocessing is fitted within each OUTER training fold,
  but the imputer/scaler precede ElasticNetCV's INNER splitting. Thus this is a
  controlled replication of the legacy implementation, not a newly corrected
  fully nested preprocessing pipeline. Silently changing this would prevent a
  one-factor sensitivity comparison.
- Stable/shared protein panels are selected before PD cross-validation rather
  than reselected inside each PD training fold. Do not call the entire workflow
  leakage-free or external validation.
- Cox C-indices and continuous NRI are evaluated in their fitted cohorts and
  are apparent performance estimates. This especially matters for large panels
  relative to the number of PD events. NRI uses the original analytic formula,
  without a censoring-adjusted reclassification estimator or bootstrap CIs.
- Cox models preserve the original 20-PC and proteomic technical covariates.
  Some original near-term Cox fits warn of potentially infinite coefficients;
  numerical reproduction does not eliminate this model-instability concern.
- Original LightGBM's long-term window is AFTER FIVE YEARS, not limited to
  years 5-10. Original shared-protein Cox outputs use within-5y and within-10y
  horizons. These definitions must not be conflated in manuscript wording.
- LightGBM's Youden thresholds are selected from the evaluated predictions as
  in the original code. Threshold-derived sensitivity/specificity are exploratory,
  not performance of a prospectively frozen threshold. Primary comparisons use
  AUC, AUPRC, and Brier score.
- LightGBM execution defaults to four threads for resource control. The original
  baseline metrics are explicitly checked after rerunning; a thread-related or
  software-related numerical difference cannot be silently treated as reproduced.
- No PCA model is fitted or used. The original shared-analysis derived CSV is
  read only for IDs and its 25 individual standardized protein columns. Its PCA
  score/loadings are never read into a model. Original snapshots with obsolete
  PCA/LASSO branches are retained unchanged solely for provenance.

## Files and Run Order

- `common.py`: isolated paths and read-only loading of original source snapshots.
- `01_prepare_and_validate.py`: matched target/data arrays, original-score and
  fold-prediction reconstruction, original fold assignments, input hashes.
- `02_fit_elastic_net_fold.py FOLD`: one complete outer fold; saves coefficients,
  OOF predictions, metrics, and warnings before marking completion.
- `03_collect_elastic_net.py`: collect all folds, select the new stable panel,
  rebuild scores, and compare panel membership with the original 374.
- `04_cox_prediction.R baseline|new`: original Cox reproduction, shared-25 vs
  shared-24 comparison, and new stable-panel/single-PPS C-index/NRI.
- `05_lightgbm_prediction.py baseline|new`: checkpointed original and sensitivity
  LightGBM runs, preserving original phenotype and covariate construction.
- `06_validate_and_report.py`: output checks and numerical comparison report.
- `07_plot_figure5.py`: reproduce all nine new stable-panel OOF AUCs and redraw
  Figure 5 using the existing gray/teal/orange triptych with inset AUC bars.
  This script does not refit models or change saved predictions. Panel C is
  correctly labelled after five years, not within ten years. Legends identify
  the multivariable 370-protein panel rather than a single scalar PPS.
- `source_reference/`: unmodified original PA5 source scripts; do not execute
  these directly because their relative output paths target their own directory.

From the repository root, configure local private input directories:

```powershell
$env:PD_PAIN_DATA_DIR = 'PATH_TO_PRIVATE_DATA'
$env:PD_PAIN_PREDICTION_SOURCE_DIR = 'PATH_TO_PRIVATE_PREDICTION_REFERENCE'
$env:PD_PAIN_LGBM_THREADS = '4'
python 05_prediction/01_prepare_and_validate.py
Rscript 05_prediction/04_cox_prediction.R baseline
python 05_prediction/05_lightgbm_prediction.py baseline
1..10 | ForEach-Object { python 05_prediction/02_fit_elastic_net_fold.py $_ }
python 05_prediction/03_collect_elastic_net.py
Rscript 05_prediction/04_cox_prediction.R new
python 05_prediction/05_lightgbm_prediction.py new
python 05_prediction/06_validate_and_report.py
python 05_prediction/07_plot_figure5.py
```

R dependencies: `data.table`, `survival`. Python dependencies: numpy, pandas,
scipy, scikit-learn, lightgbm, and threadpoolctl. The new branch was tested using
Python 3.12.4, numpy 2.0.0, pandas 2.2.3, scipy 1.15.2, scikit-learn 1.8.0, and
LightGBM 4.6.0. R 4.4.1 uses the locally installed survival package. Version
information, input hashes, fold seeds, and model warnings are retained.

Original PA5 inputs include the protein matrix and original pain scores; saved
Elastic Net coefficients/fold metrics/PPS; saved original Cox/NRI and LightGBM
metrics; `protein_intersection_FDR05.csv`; and the derived individual shared-
protein columns in `overlap25_directional_pca_signature.csv`. Raw UKB outcomes,
covariates, technical factors, and protein timing come from `PD_PAIN_DATA_DIR`.
New pain scores/shared proteins come from branch-07 output. Originals are read
only; no source data or published result is overwritten.

## Outputs and Resume

All outputs are under
`results/sensitivity_nonproteomic_weights/prediction/`. Participant-level arrays,
OOF predictions, design matrices, and intermediate cohorts must remain private
and excluded from GitHub. Allow several GB of disk and memory. Do not reuse
checkpoints after changing input scores, folds, data, feature sets, or parameters.

Elastic Net folds are independently checkpointed in `elastic_net_folds/fold_XX/`.
Only missing folds need to be launched after an interruption; a completed fold
is never overwritten automatically. All 10 must exist before collection.
At most two outer-fold processes were used together on the local 16-GB machine.
Later folds limit inner execution to two jobs rather than four to bound working
memory; this changes execution resources, not the CV grid, folds, or model.
LightGBM checkpoints every fold and validates participant order/outcomes when
resuming. It does not silently replace failed fits or discard convergence warnings.

Figure 5 plotting additionally requires matplotlib and Arial. The plotting script
defaults to `figures/figure5_new_prediction/` alongside the repository folder, and accepts
`--input-dir`, `--output-dir`, and `--dpi` (default 600). It exports the combined
figure as PNG/TIFF/PDF/SVG, separate A/B/C panels, an English caption, aggregate
ROC coordinates, a nine-row metric table, and validation/QA records. Only
aggregate curves and metrics are exported; participant-level predictions remain
in the private analysis results. SVG/PDF text remains editable. The original
manuscript image and all original analysis outputs are left untouched.

Independent LightGBM windows can be run with `--task incident_pd_2to5y`,
`--task incident_pd_after5y`, and `--task incident_pd_full_followup`; after all
three finish, use `05_lightgbm_prediction.py new --collect-only` to combine them.

Key outputs: `PREDICTION_SENSITIVITY_RESULTS.md`, `prediction_comparison.csv`,
`validation_checks.csv`, original-reproduction check tables, new fold metrics,
`new_elastic_net_oof_predictions.csv`, `new_stable_elastic_net_pfs.csv`,
`new_stable_elastic_net_coefficients.csv`, stable-panel retained/lost/gained
membership, shared/stable C-index and NRI tables, and original/new LightGBM
OOF predictions, fold metrics, and cohort files.
