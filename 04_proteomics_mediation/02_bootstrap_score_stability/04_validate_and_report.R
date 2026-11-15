source("04_proteomics_mediation/02_bootstrap_score_stability/00_bootstrap_helpers.R")
prepared <- fread(boot_output("preparation_checks.csv"))
diagnostics <- fread(boot_output("bootstrap_fit_diagnostics.csv"))
coefs <- readRDS(boot_output("bootstrap_coefficients.rds"))
target <- readRDS(boot_output("bootstrap_scoring_input.rds"))
saved_scores <- readRDS(boot_output("proteomic_scores_1000.rds"))
matrices <- readRDS(boot_output("pairwise_correlation_matrices.rds"))
summary <- fread(boot_output("correlation_summary.csv"))
calculation_checks <- fread(boot_output("correlation_calculation_checks.csv"))
checks <- list()
check <- function(name, passed) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(passed))
}
check("exactly_1000_successful_replicates", nrow(diagnostics) == 1000L && all(diagnostics$valid))
check("all_1000_replicate_IDs", identical(diagnostics$replicate, seq_len(1000L)))
check("all_resamples_same_size_as_derivation", all(diagnostics$n_sampled == prepared$n_derivation))
check("sampling_with_replacement_duplicates_present", all(diagnostics$n_unique < diagnostics$n_sampled))
check("ordinary_bootstrap_event_count_varies", length(unique(diagnostics$events)) > 1L)
# The binary diagnostics preserve empty strings; CSV readers can infer an
# all-empty warning/error column as logical NA, which is not a model warning.
check("no_fit_warnings_or_errors", all(coefs$diagnostics$warnings == "" & coefs$diagnostics$error == ""))
check("reference_standard_matrix_fit_agree", prepared$baseline_coxph_matrix_max_beta_difference < 1e-10)
check("bootstrap_standard_matrix_fit_agree", prepared$bootstrap_coxph_matrix_max_beta_difference < 1e-8)
check("derivation_protein_disjoint", prepared$derivation_protein_overlap == 0L)
check("same_target_ids_and_order", identical(saved_scores$eid, target$eid))
check("score_matrix_dimensions", identical(dim(saved_scores$scores), c(length(target$eid), 1000L)))
check("all_scores_finite_and_bounded", all(is.finite(saved_scores$scores)) &&
        min(saved_scores$scores) >= 0 && max(saved_scores$scores) <= 1 + 1e-12)
selected <- c(1L, 2L, 500L, 999L, 1000L)
expected <- target$indicators %*% coefs$normalized_weights[, selected, drop = FALSE]
check("saved_scores_recalculate_from_saved_betas", max(abs(saved_scores$scores[,selected] - expected)) < 1e-12)
check("negative_raw_betas_clipped_exactly", max(abs(coefs$effective_beta -
        pmax(coefs$raw_beta[rownames(coefs$effective_beta), , drop = FALSE], 0))) == 0)
check("normalized_weights_sum_to_one", max(abs(colSums(coefs$normalized_weights) - 1)) < 1e-12)
check("four_correlation_matrices", length(matrices) == 4L)
for (name in names(matrices)) {
  m <- matrices[[name]]
  check(paste0(name, "_valid_matrix"), identical(dim(m), c(1000L,1000L)) && all(is.finite(m)) &&
          max(abs(m-t(m))) < 1e-10 && max(abs(diag(m)-1)) < 1e-12 && max(abs(m)) <= 1 + 1e-10)
}
check("499500_unique_pairs_per_matrix", all(summary[comparison == "bootstrap_vs_bootstrap", n_correlations] == choose(1000L, 2L)))
check("1000_reference_comparisons_per_population_method", all(summary[comparison == "bootstrap_vs_no_bootstrap", n_correlations] == 1000L))
check("direct_correlation_validation_passed", all(calculation_checks$passed))
source_hashes <- fread(boot_output("input_source_hashes.csv"))
check("previous_sensitivity_inputs_unchanged", identical(source_hashes$md5, unname(tools::md5sum(source_hashes$path))))
result <- rbindlist(checks)
fwrite(result, boot_output("bootstrap_validation_checks.csv"))
if (!all(result$passed)) stop("Bootstrap validation failed; inspect bootstrap_validation_checks.csv.")

display_rows <- function(comparison) {
  comparison_name <- comparison
  subset <- summary[get("comparison") == comparison_name]
  vapply(seq_len(nrow(subset)), function(i) sprintf(
    "| %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f |",
    subset$population[i], subset$method[i], subset$median[i], subset$minimum[i],
    subset$p025[i], subset$p975[i], subset$maximum[i]), "")
}
report <- c(
  "# Bootstrap Stability of Composite Pain Scores", "", "## Material Passport",
  "- Status: COMPLETED; 1000/1000 bootstrap fits valid and all implemented checks passed.",
  "- Mode: participant-level nonparametric bootstrap of the non-proteomic derivation cohort.",
  "- Source: unchanged local PA model specifications and prior sensitivity cohorts.",
  "- Main analyses, prior sensitivity results, and manuscript were not overwritten.", "",
  "## Design",
  sprintf("- Derivation participants: %s; reference PD events: %s.",
          format(prepared$n_derivation, big.mark=","), format(prepared$derivation_events, big.mark=",")),
  sprintf("- Fixed scoring cohort: %s, including %s participants with any reported pain.",
          format(prepared$n_scored, big.mark=","), format(prepared$n_scored_with_pain, big.mark=",")),
  "- Each replicate sampled the full derivation cohort with replacement and refitted the original Cox model.",
  "- Negative pain coefficients were set to zero, following the original PA score implementation.",
  "- Each replicate's coefficient sum was used to normalize its weights; these weights were then fixed for scoring.",
  "- Every replicate scored the same proteomic participants without using their PD outcomes.",
  "- The reference was the unbootstrapped non-proteomic-weight score from sensitivity branch 07.",
  "- All pairwise comparisons: 499,500 distinct score-vector pairs for each population and method.",
  "- Random seed: 20261026; saved replicate-specific seeds make results independent of worker scheduling.", "",
  "## Pairwise Bootstrap Score Correlations",
  "| Population | Method | Median | Minimum | 2.5th percentile | 97.5th percentile | Maximum |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
  display_rows("bootstrap_vs_bootstrap"), "",
  "## Correlations With the Unbootstrapped Score",
  "| Population | Method | Median | Minimum | 2.5th percentile | 97.5th percentile | Maximum |",
  "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
  display_rows("bootstrap_vs_no_bootstrap"), "",
  "## Interpretation and Limitations",
  sprintf("Median pairwise Pearson correlations were %.4f in all scored participants and %.4f in those with pain;",
          summary[comparison == "bootstrap_vs_bootstrap" & population == "all_proteomic" & method == "pearson", median],
          summary[comparison == "bootstrap_vs_bootstrap" & population == "with_any_pain" & method == "pearson", median]),
  sprintf("the corresponding Spearman medians were %.4f and %.4f. These results do not support describing the",
          summary[comparison == "bootstrap_vs_bootstrap" & population == "all_proteomic" & method == "spearman", median],
          summary[comparison == "bootstrap_vs_bootstrap" & population == "with_any_pain" & method == "spearman", median]),
  "scores or within-pain rankings as uniformly highly stable. Similarity to the single unbootstrapped",
  "reference is higher than similarity among arbitrary pairs of bootstrap scores. The lower within-pain",
  "correlations show that including consistently zero-scored participants contributes to the overall agreement.",
  "",
  "The correlations describe sensitivity of scores and participant rankings to the estimation of pain weights.",
  "They are conditional on the fixed proteomic cohort and its pain phenotypes. They do not establish external",
  "prediction performance, causal effects, or clinical utility. The with-pain analysis excludes participants",
  "whose 14 pain indicators are all absent, limiting the contribution of universally zero scores.",
  "Percentiles summarize the empirical distribution of correlations; they are not confidence intervals for",
  "an independently sampled correlation parameter. Pairwise correlations share replicates and are dependent.",
  "Individual pain weights can remain uncertain even when the composite scores are strongly correlated;",
  "the beta and weight summaries therefore report coefficient dispersion and the frequency of zero weights.",
  "This bootstrap varies the fitted model coefficients, not the original cohort eligibility or covariate",
  "cutpoints. No protein associations, shared-set selection, mediation, or prediction models were rerun.", "",
  "## Outputs and Validation",
  "- all_model_beta_1000_columns.csv: one coefficient column per bootstrap fit, including adjustment terms.",
  "- pain_raw_beta_1000_columns.csv: the 14 raw pain coefficients, one column per bootstrap fit.",
  "- pain_effective_beta_1000_columns.csv: the 14 pain coefficients after truncation.",
  "- pain_normalized_weight_1000_columns.csv: the 14 normalized weights used for scoring.",
  "- proteomic_scores_1000.rds / proteomic_scores_1000_columns.csv.gz: participant scores, one column per replicate.",
  "- correlation_matrix_*.csv.gz / pairwise_correlation_matrices.rds: four complete 1000 x 1000 matrices.",
  "- correlation_summary.csv / correlations_with_no_bootstrap_reference.csv: aggregate and replicate-level correlations.",
  "- individual_score_uncertainty.csv: mean, SD, and percentile summaries for each scored participant.",
  "- raw_beta_bootstrap_summary.csv / pain_weight_bootstrap_summary.csv: coefficient and weight uncertainty.",
  "- bootstrap_fit_diagnostics.csv: seeds, event counts, unique sampled participants, warnings, and errors.",
  sprintf("- bootstrap_validation_checks.csv: all %d checks passed.", nrow(result)),
  "- correlation_calculation_checks.csv: reduced calculations checked against direct participant-level correlations.",
  "- Execution used four independent R processes. An initial local-socket cluster setup was blocked by the environment",
  "  before any bootstrap fit ran; its log is retained under setup_diagnostics/. No statistical specification changed.",
  "Participant IDs and score matrices are confidential outputs and must remain excluded from GitHub.")
writeLines(report, boot_output("BOOTSTRAP_STABILITY_RESULTS.md"))
writeLines(capture.output(sessionInfo()), boot_output("sessionInfo_validation.txt"))
cat("All", nrow(result), "bootstrap validation checks passed.\n")
print(summary)
