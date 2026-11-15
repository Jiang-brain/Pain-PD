# Bootstrap Stability of Non-Proteomic-Derived Pain Scores

This analysis follows the independent pain-score workflow without changing its data, scores, protein
associations, or shared-protein results. It quantifies uncertainty in pain-score
weights and the resulting composite scores in a fixed proteomic cohort.

## Design

- Resample the 408,583 non-proteomic derivation participants with replacement,
  drawing 408,583 rows per replicate, for exactly 1,000 replicates.
- Refit the original mutually adjusted Cox model with Efron handling of ties.
  Resampled duplicates are actual repeated rows, not Cox case weights: the two
  are not generally equivalent under Efron ties.
- As in the PA source, truncate negative pain coefficients to zero, normalize by
  the sum of the retained coefficients, and freeze these weights within each
  replicate. The covariate definitions and cohort eligibility are not re-derived
  within each bootstrap; the sensitivity concerns coefficient estimation.
- Apply each set of weights to the same 49,056 proteomic participants. Their PD
  outcomes never enter the bootstrap fitting or scoring calculation.
- Calculate all 499,500 distinct pairwise correlations among the 1,000 score
  vectors, plus each vector's correlation with the unbootstrapped branch-07 score.
- Report Pearson and Spearman correlations in all 49,056 participants and in the
  29,417 whose pain indicators are not all absent. This second analysis limits
  the influence of participants assigned zero under every set of weights.
- Save raw coefficients, truncated coefficients, normalized weights, all scores,
  coefficient/score uncertainty, and complete correlation matrices.

The 2.5th and 97.5th percentiles of correlations describe their empirical
distribution; they are not confidence intervals for an independent correlation
parameter. The pairwise correlations are dependent. High score correlations
alone do not establish coefficient stability, external validity, or clinical
utility. No additional protein associations or shared-set analyses are run here.

## Files

- `00_bootstrap_helpers.R`: matrix Cox fitting and exact correlation helpers.
- `01_prepare_bootstrap_inputs.R`: recreate the original PA design matrix;
  reproduce branch-07 coefficients and scores; validate matrix Cox fitting
  against `coxph` both before and after actual participant resampling; save inputs.
- `02_run_bootstrap.R`: fit/checkpoint independent deterministic bootstrap
  replicates and collect all 1,000 coefficient vectors.
- `03_score_correlations.R`: calculate/save individual scores, full Pearson and
  Spearman matrices, empirical summaries, and coefficient uncertainty.
- `04_validate_and_report.R`: integrity checks and an automatically generated
  `BOOTSTRAP_STABILITY_RESULTS.md` based on the actual completed results.
- `05_plot_spearman_stability.py`: visualize the completed all-proteomic
  Spearman matrix and the distribution of its 499,500 unique off-diagonal
  correlations. This optional Python plotting step does not refit any model.
- `test_correlation_helpers.R`: synthetic mathematical unit tests, not study
  results; verifies weighted-tie and reduced-covariance correlation calculations.
- `run_bootstrap.ps1`: optional Windows launcher using hidden independent R
  processes. It does not require local sockets, clusters, or network access.

## Run

Requirements: completed independent pain-score analysis, R, `data.table`, and `survival`.
No additional package installation or internet connection is needed. From the repository root:

```powershell
$env:PD_PAIN_DATA_DIR = 'PATH_TO_PRIVATE_DATA'
# Use your installed Rscript executable if it is not on PATH.
.\04_proteomics_mediation\02_bootstrap_score_stability\run_bootstrap.ps1 -Rscript 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' -Workers 4
```

Alternatively, run sequentially:

```text
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/test_correlation_helpers.R
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/01_prepare_bootstrap_inputs.R
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/02_run_bootstrap.R
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/02_run_bootstrap.R --collect-only
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/03_score_correlations.R
Rscript 04_proteomics_mediation/02_bootstrap_score_stability/04_validate_and_report.R
python 04_proteomics_mediation/02_bootstrap_score_stability/05_plot_spearman_stability.py
```

For manual independent-process execution, pass `SHARD_INDEX SHARD_COUNT` to step
02, for example `1 4`, `2 4`, `3 4`, and `4 4`, using the same total shard count.
Do not launch overlapping shard assignments. Run `--collect-only` after all
processes finish. Do not run the launcher while another set of shards is active.

The master seed is 20261026. Each replicate has a saved unique seed; execution
order, scheduling, and shard count do not change its sample. Completed atomic
checkpoints can be resumed; source/input hashes must match. Failed or warned
fits are retained and reported, never silently discarded or replaced. Correct
the cause only after reviewing the diagnostic checkpoint; do not claim 1,000
valid replicates if fewer succeeded.

## Output and Privacy

Outputs go to `results/sensitivity_nonproteomic_weights/bootstrap1000/`, or the
same subdirectory under `PD_PAIN_RESULTS_DIR`. The uncompressed R score matrix
uses approximately 392 MB; allow several GB for working memory and outputs.
Four workers were selected for the local 16-GB machine. Reduce workers when other
applications need memory. All source definitions and prior outputs remain intact.

Important outputs:

- `all_model_beta_1000_columns.csv`: all adjustment/pain coefficients, one column
  per bootstrap; `pain_raw_beta_1000_columns.csv` contains only the 14 raw pain
  coefficients; `pain_effective_beta_1000_columns.csv` and
  `pain_normalized_weight_1000_columns.csv` hold the scoring coefficients.
- `proteomic_scores_1000.rds`: a list containing `eid` and a 49,056 x 1,000
  double-precision `scores` matrix. A wide `.csv.gz` version is also saved.
- `correlation_summary.csv`: primary numerical summaries for both populations,
  both correlation methods, and pairwise/reference comparisons.
- `correlation_matrix_*.csv.gz` and `pairwise_correlation_matrices.rds`: four
  complete 1,000 x 1,000 correlation matrices with replicate names.
- `correlations_with_no_bootstrap_reference.csv`: 1,000 correlations per
  population/method versus the original branch-07 fixed-weight score.
- `individual_score_uncertainty.csv`: each participant's mean, SD, and empirical
  score percentiles; these are score-estimation summaries, not PD risk intervals.
- `raw_beta_bootstrap_summary.csv`, `pain_weight_bootstrap_summary.csv`: raw
  coefficient and normalized-weight dispersion, including zero-weight frequency.
- `bootstrap_fit_diagnostics.csv`, `bootstrap_validation_checks.csv`, and
  `correlation_calculation_checks.csv`: model and mathematical integrity checks.
- `BOOTSTRAP_STABILITY_RESULTS.md`: the completed-results report.

All data, IDs, design matrices, and score outputs are private and remain under
the ignored `results/` directory; do not upload them to GitHub.

Pearson correlations are calculated algebraically from the 14-indicator
covariance matrix. Spearman correlations use exact frequency-weighted midranks
for identical pain patterns. Both are checked against direct participant-level
`cor()` on selected replicates, including ties, separately in each population.
These are computational reductions, not sampling approximations.

## Spearman Figure

The optional plotting script requires Python, numpy, pandas, matplotlib, and
Arial. It validates the matrix, all 1,000 successful fits, and reproduction of
the saved summary before exporting. Defaults are the completed bootstrap output
directory for inputs and `figures/bootstrap_spearman_stability/` alongside the repository folder
for outputs; `--input-dir` and `--output-dir` can override them.

Exports include PNG/TIFF (600 dpi), editable-text PDF/SVG, a preview, an English
Results sentence and caption, aggregate matrix/histogram source data, a Spearman
summary, and integrity records. The figure includes all 49,056 scoring
participants, not the pain-only subgroup. Its two panels show the 1,000-by-1,000
correlation matrix (without diagonal self-correlations) and the 499,500 unique
pairwise correlations. The empirical percentile range is not a confidence
interval. No patient IDs, pain records, or participant-level scores are exported.
