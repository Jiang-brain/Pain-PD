# Sensitivity: Derive Pain Weights Outside the Proteomic Cohort

## Scope

Originally developed as a sensitivity analysis, this branch now supplies the
independently weighted pain score adopted for the final Figure 4/5 proteomic
analyses. It preserves original files and reassesses pain-protein associations
and shared-protein selection underlying Figure 4A/D. Downstream mediation is
in the shared-protein mediation workflow and prediction in the prediction workflow; this branch does not perform those
steps itself. The original protein-PD association models are unchanged.

## Methods and Preservation of the Original Analysis

1. Reconstruct the eligible score-development cohort using an unchanged snapshot
   of `PD_Pain_WQ.R` from the final PA analysis. This preserves its variable
   recodings, missing-value handling, generalised-pain exclusions, age restriction,
   follow-up definition, and exclusion of PD events before 730 days.
2. Reproduce the original all-cohort scores and compare participant IDs and every
   score against `composite_pain_scores.csv`. Stop if IDs differ or the maximum
   absolute score difference is at least `1e-12`. This is a validation gate only;
   all-cohort coefficients are not used to compute the sensitivity scores.
3. Exclude all participants listed in the unfiltered protein matrix or supplied
   protein roster, plus all participants listed in the filtered protein matrix,
   from weight development. This excludes measured participants even when they
   failed protein QC. No participant can belong to both groups.
4. Fit the original mutually adjusted Cox formula in the remaining participants.
   Extract the 14 acute/chronic pain coefficients, truncate negative coefficients
   to zero, and freeze these coefficients and their sum. For each eligible
   proteomic participant, compute `sum(indicator * fixed_beta) / sum(fixed_beta)`.
   Any supplementary tertile boundaries are fitted in the derivation group only;
   the protein association uses the continuous score, not tertiles.
5. Rerun the unchanged PA protein-wise linear models, replacing only the score
   file and output paths. Preserve the same protein universe, per-protein timing
   covariate, technical factors, main covariates, and 20 genetic PCs. Apply BH FDR
   and Bonferroni adjustment across all tested proteins. Require protein-specific
   sample sizes to match the original analysis.
6. Reuse the original `protein_PD_cox_results.csv`: the new pain weights do not
   alter these protein-PD models. Report both the full significant intersection
   and its directionally concordant subset: `pain_beta * PD_log_HR > 0`. The latter
   is the shared-protein criterion reported for the manuscript's original 25.
   Compare new membership with a baseline set recalculated from original results.

Important source-specific details are preserved rather than harmonised:

- The PA weight-development Cox model includes age, sex, ethnicity, deprivation,
  income, education, smoking, alcohol, sedentary behaviour, obesity, and PRS;
  genetic PCs are commented out in that particular source formula.
- The source truncates negative pain coefficients to zero, although the current
  manuscript's composite-score paragraph does not explicitly state this step.
- The actual PA association script regresses `pain_score` on each protein, not
  each protein on `pain_score` as described in the current manuscript. The
  sensitivity preserves the implemented model; this discrepancy needs author
  confirmation before manuscript revision. Beta magnitudes must not be described
  as change in protein level per unit pain score under the implemented model.
- Original score-development and protein-association recodings are not identical
  in all respects. Each is preserved from its own source snapshot.

## Scripts

- `00_common.R`: independent path configuration and integrity helpers.
- `01_derive_weights_and_score.R`: original-score reproduction check, disjoint
  derivation cohort, fixed pain weights, and proteomic-participant scores.
- `02_pain_score_protein_association.R`: original PA protein-association models
  using the new scores, with protein-specific sample-size checks.
- `03_shared_proteins.R`: significant intersections, directional selection,
  and retained/lost/gained comparison with the original 25-protein selection.
- `04_validate_outputs.R`: independent raw-field score recalculation, disjoint
  cohorts, unchanged protein samples, multiple-testing correction, and shared-set
  integrity checks.
- `source_reference/`: unchanged copies of the two final PA source scripts.
  Do not execute these snapshots directly; the sensitivity drivers evaluate
  selected parsed R expressions with explicit input/output path overrides.

## Run

R 4.4.1 was used locally. Required packages: `data.table`, `survival`.
From the repository root, with `Rscript` available on PATH:

```powershell
$env:PD_PAIN_DATA_DIR = 'PATH_TO_PRIVATE_DATA'
# Optional: parent directory for sensitivity outputs (defaults to results/).
# $env:PD_PAIN_RESULTS_DIR = 'PATH_TO_PRIVATE_RESULTS'
# If protein_new.csv is unavailable, provide the complete measured roster:
# $env:PD_PAIN_PROTEIN_ROSTER = 'PATH_TO_PRIVATE_ROSTER_CSV'
Rscript 04_proteomics_mediation/01_independent_pain_score/01_derive_weights_and_score.R
Rscript 04_proteomics_mediation/01_independent_pain_score/02_pain_score_protein_association.R
Rscript 04_proteomics_mediation/01_independent_pain_score/03_shared_proteins.R
Rscript 04_proteomics_mediation/01_independent_pain_score/04_validate_outputs.R
```

Required private inputs: `PD_data.csv`, `composite_pain_scores.csv`,
`protein_new_missing60_filtered.csv`, `protein_new.csv` (or complete measured
roster with an `eid` column), `PD_data_matched_protein_covariates.csv`,
`PD_data_matched_protein_times.csv`, `pain_score_protein_association_results.csv`,
and `protein_PD_cox_results.csv`. None is included with the code.

## Outputs

All outputs go to `results/sensitivity_nonproteomic_weights/` (or the same
subdirectory under `PD_PAIN_RESULTS_DIR`), never to the original PA directory.

Main outputs:

- `cohort_counts.csv`, `fixed_weight_checks.csv`: counts and disjointness check.
- `baseline_score_reproduction_check.csv`: original-score validation gate.
- `pain_score_fixed_beta_comparison.csv`, `fixed_weights.rds`: original and
  derivation-only coefficients, frozen denominator, and model formula.
- `composite_pain_scores.csv`: new scores for eligible proteomic participants.
- `proteomic_score_comparison.csv`: paired original and sensitivity scores.
- `pain_score_protein_association_results.csv`: all tested proteins.
- `pain_score_protein_significant_FDR05.csv`: FDR-significant pain proteins.
- `protein_association_counts.csv`, `protein_association_sample_size_check.csv`:
  model counts, failures, and per-protein sample-size checks.
- `shared_proteins_fdr05_concordant.csv`: primary sensitivity shared set.
- `protein_intersection_fdr05_all_directions.csv`: unfiltered-direction overlap.
- `shared_protein_retained_lost_gained.csv`: comparison with the original set.
- `shared_protein_summary.csv`: original and sensitivity discovery counts.
- `validation_checks.csv`: pass/fail results of the independent output checks.
- `*_hash*.csv`, `sessionInfo_*.txt`: source provenance and R environment.

The output directory includes participant IDs and individual scores. Keep it
private and excluded from GitHub, just like all other participant-level outputs.
This sensitivity assesses whether weights derived outside the proteomic cohort
change the discoveries; it does not establish causality or external validity.
