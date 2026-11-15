# Mediation With Independently Weighted Pain Scores

Release status: this branch supplies the adopted shared-24 Figure 4 analysis.
Final plotting entry point: `06_figures/02_plot_final_figure4.py`. For the right
side of E use `08_plot_figure4E_forest_compact.py` (large, closely spaced circles),
not the earlier small-marker composite. G uses `07_plot_figure4G_top5.py`.

This branch replaces the original pain score with the fixed non-proteomic-derived
score from the independent pain-score workflow and replaces the original shared25 set with shared24 (COMP
excluded). It does not rerun PPI or enrichment. Original source files, data and
results are read-only, and outputs are isolated in
`results/sensitivity_nonproteomic_weights/mediation24/`.

## Controlled Analysis

- Retain the original 46,591 mediation participants and 707 incident PD events.
- Retain their original aligned, standardized, median-imputed protein values,
  technical factors, 20 ancestry PCs and original covariate recodings.
- Reproduce the saved original observed CHGA, ITGAV and joint25 results before
  fitting new-score models. Stop on disagreement.
- Fit all 24 single-mediator models and the joint24 IORW model using the original
  formula functions. Gaussian stabilized density-ratio weights are trimmed at
  their 99th percentile and normalized to mean one. Cox ties use Efron handling.
- Keep the new pain score fixed during mediation bootstrapping. This is distinct
  from the score-stability workflow, which resamples the non-proteomic score-development cohort.
- Run 1,000 bootstrap resamples per single protein plus 1,000 for the joint set
  (25,000 attempts total). Preserve the original serial sample.int RNG streams,
  seed 20260602, and original protein order with COMP removed from the end.
- Checkpoint every 20 attempts. Resume only after checking input hashes; never
  overwrite or substitute failed resamples. Each error/warning is retained.
- Calculate 95% percentile CIs for proportions mediated and indirect log-HRs.
  Final reporting uses indirect log-HR CI exclusion of zero: 16/24 in the saved
  results. These intervals are not multiplicity-adjusted. Preserve original
  one-sided and two-sided tail/BH columns for historical comparison only;
  the previous 18/24 one-sided-FDR count is not the final criterion.

For speed, bootstrap coefficient-only fits use the same `lm.fit` and
`survival::coxph.fit` engines as the original formula models. They omit repeated
robust-variance, concordance and naive-adjusted calculations, none of which is
needed to obtain percentile CIs. Duplicate censored rows can be collapsed with
equivalent multiplicity weights, but **all duplicated event rows remain separate**
to preserve Efron ties. Both paths are checked numerically against the original
formula fits and actual repeated-row resamples. Early checkpoints computed with
all repeated rows are equally valid. The computational reduction does not alter
the sampled participants, mediators, covariates, score or statistical estimator.

## Important Limits

- Observational mediation cannot establish causality or intervention benefit.
- The original deprivation/PRS recoding collapses both factors to constants;
  the original usable_covariates helper omits them. This branch preserves and
  audits that limitation rather than silently correcting the comparison.
- A bootstrap tail probability of zero is not an exact P value of zero. Legacy
  one-sided and two-sided FDR columns are distinct and are not the final CI
  criterion. For an indirect HR interval the null is one, not zero.
- A finite `valid` effect estimate is not a claim that every nuisance coefficient
  is stable. Bootstrap Cox fits can warn that batch or other coefficients may
  be infinite; these warnings are retained, counted and compared against the
  original formula implementation in `test_warned_resamples.R`. Do not describe
  these runs as uniformly warning-free or fully converged.
- Shared-protein selection is not repeated inside mediation bootstraps. CIs are
  conditional on the selected mediators and the frozen independently derived
  pain-score weights, rather than uncertainty intervals for the full pipeline.
- Source data and all participant-level analysis tables stay private. Do not
  upload the ignored results directory to GitHub.

## Run Order

Run from the repository root. R requires `data.table` and `survival`; source functions and
saved original mediator data come from the configured reference directory.
The launcher is standard-library Python. Plotting additionally uses numpy,
pandas, scipy, matplotlib, Pillow and Arial. The three original R source helpers
are bundled unchanged in `source_reference/`; they are parsed for functions,
not executed as standalone analysis scripts. Private reference cohorts/results
still come from the configured reference directory.

```powershell
$env:PD_PAIN_PREDICTION_SOURCE_DIR = 'PATH_TO_PRIVATE_PREDICTION_REFERENCE'
Rscript 04_proteomics_mediation/03_shared_protein_mediation/01_prepare_and_validate.R
Rscript 04_proteomics_mediation/03_shared_protein_mediation/test_bootstrap_acceleration.R
python 04_proteomics_mediation/03_shared_protein_mediation/run_mediation.py --rscript 'PATH_TO_RSCRIPT' --reference-dir $env:PD_PAIN_PREDICTION_SOURCE_DIR --workers 8
Rscript 04_proteomics_mediation/03_shared_protein_mediation/04_path_coefficients.R
python 04_proteomics_mediation/03_shared_protein_mediation/05_plot_figure4_panels.py
python 04_proteomics_mediation/03_shared_protein_mediation/08_plot_figure4E_forest_compact.py
python 04_proteomics_mediation/03_shared_protein_mediation/07_plot_figure4G_top5.py
```

`run_mediation.py` prepares deterministic job streams if needed, launches hidden
local R workers, waits for all workers, and collects results automatically.
Do not run multiple launchers concurrently. Interrupted atomic checkpoints can
be reused; do not change input files, score definitions, folds or source models
when resuming. Worker logs are appended, not overwritten.

## Outputs

- `per_protein_mediation_results.csv`: 24 observed fits, bootstrap CIs, primary
  `indirect_95CI_excludes_zero`/`proportion_95CI_excludes_zero` flags, and legacy
  one-sided/two-sided FDR columns retained for comparison.
- `joint24_mediation_results.csv`: joint observed fit and percentile CIs.
- `mediation24_bootstrap.csv`, `bootstrap_counts.csv`: all 25,000 attempted
  fits and per-model success/warning/error counts.
- `reproduction_and_implementation_checks.csv`,
  `censored_duplicate_equivalence_checks.csv`,
  `mediation_validation_checks.csv`: original-result and numerical checks.
- `covariates_actually_used.csv`, `input_hashes.csv`, session information and
  source hashes: transparent model/input provenance.
- `per_protein_standardized_path_coefficients.csv`,
  `joint24_standardized_path_coefficients.csv`: descriptive path coefficients
  for the replacement Figure4E diagram, using the original formula functions.
- `MEDIATION_RESULTS.md`: automatically generated numerical report and limits.

Replacement Figure4 A/D/E panels default to
`figures/figure4_independent_weights/` alongside the repository folder, with 600-dpi PNG/TIFF and
editable-text PDF/SVG. A/D can be exported before bootstrap completion with
`--associations-only`. E requires all 1,000 successful resamples per model.
Original B/C and author-regenerated F/G are not combined silently with these
replacement panels.

The supplied G inputs are bundled under `figure_data/figure4/`: five terms per
GO category and per pathway database by ascending FDR, with ties retaining
export order. This does not rerun STRING. Earlier whole-figure assembly scripts
are not included in the submission tree.

`test_ci_reporting.R` recomputes summaries from saved bootstrap effects and
checks that estimates/CIs are unchanged and the final count is 16. It does not
refit models, resample participants or write analysis results. Existing output
CSVs can be read by final plotting even before the new CI flag columns are added.
Running `03_collect_and_report.R` updates the summary/flags from existing
checkpoints; it does not repeat the bootstrap fits.
