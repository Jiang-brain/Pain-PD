source("04_proteomics_mediation/02_bootstrap_score_stability/00_bootstrap_helpers.R")
coefficient_path <- boot_output("bootstrap_coefficients.rds")
if (!file.exists(coefficient_path)) stop("Run all 1000 bootstrap fits first.")
coefs <- readRDS(coefficient_path)
target <- readRDS(boot_output("bootstrap_scoring_input.rds"))
if (ncol(coefs$normalized_weights) != 1000L || !all(coefs$diagnostics$valid)) stop("Require 1000 valid fits.")
weights <- coefs$normalized_weights[colnames(target$indicators), , drop = FALSE]
score_matrix <- target$indicators %*% weights
colnames(score_matrix) <- colnames(weights)
if (any(!is.finite(score_matrix))) stop("Invalid bootstrap scores.")
saveRDS(list(eid = target$eid, scores = score_matrix), boot_output("proteomic_scores_1000.rds"), compress = FALSE)
fwrite(data.table(eid = target$eid, as.data.table(score_matrix)),
       boot_output("proteomic_scores_1000_columns.csv.gz"), compress = "gzip", nThread = 4L)
cat("Saved", nrow(score_matrix), "x 1000 participant scores.\n")

# Identical pain-indicator patterns have identical scores. Frequency-weighted
# midranks on these patterns give exactly the full participant-level Spearman
# correlations, including ties, without ranking 49 million entries repeatedly.
indicator_table <- as.data.table(target$indicators)
term_names <- names(indicator_table)
patterns <- indicator_table[, .(pattern_id = .GRP, frequency = .N), by = term_names]
pattern_ids <- patterns[indicator_table, on = term_names, pattern_id]
pattern_matrix <- as.matrix(patterns[, ..term_names])
all_weights <- cbind(weights, no_bootstrap_reference = target$reference_weights[colnames(target$indicators)])
pattern_scores <- pattern_matrix %*% all_weights
all_pain_patterns <- rowSums(pattern_matrix) > 0
summaries <- list()
reference_correlations <- list()
validation <- list()
correlation_matrices <- list()
add_check <- function(name, difference, tolerance = 1e-10) {
  validation[[length(validation) + 1L]] <<- data.table(check = name, max_absolute_difference = difference,
                                                    tolerance = tolerance,
                                                    passed = is.finite(difference) && difference < tolerance)
}
for (population in c("all_proteomic", "with_any_pain")) {
  participant_rows <- if (population == "all_proteomic") rep(TRUE, length(target$eid)) else target$any_pain
  pattern_rows <- if (population == "all_proteomic") rep(TRUE, nrow(patterns)) else all_pain_patterns
  if (sum(participant_rows) < 3) stop("Insufficient participants in ", population)
  pearson <- score_covariance_correlation(target$indicators[participant_rows, , drop = FALSE], all_weights)
  ranks <- vapply(seq_len(ncol(pattern_scores)), function(j) {
    weighted_midranks(pattern_scores[pattern_rows, j], patterns$frequency[pattern_rows])
  }, numeric(sum(pattern_rows)))
  colnames(ranks) <- colnames(all_weights)
  spearman <- weighted_column_correlation(ranks, patterns$frequency[pattern_rows])
  test_columns <- c(1L, 2L, 3L, 500L, 1000L)
  test_scores <- cbind(score_matrix[participant_rows, test_columns, drop = FALSE],
                       no_bootstrap_reference = target$reference_scores[participant_rows])
  selected <- c(test_columns, 1001L)
  add_check(paste0(population, "_pearson_vs_direct_cor"),
            max(abs(pearson[selected, selected] - cor(test_scores, method = "pearson"))))
  add_check(paste0(population, "_spearman_vs_direct_cor"),
            max(abs(spearman[selected, selected] - cor(test_scores, method = "spearman"))))

  for (method in c("pearson", "spearman")) {
    correlation <- if (method == "pearson") pearson else spearman
    pairwise <- correlation[seq_len(1000L), seq_len(1000L)]
    values <- pairwise[upper.tri(pairwise)]
    vs_reference <- correlation[seq_len(1000L), 1001L]
    summaries[[length(summaries) + 1L]] <- summarize_correlations(values, "bootstrap_vs_bootstrap", population, method)
    summaries[[length(summaries) + 1L]] <- summarize_correlations(vs_reference, "bootstrap_vs_no_bootstrap", population, method)
    reference_correlations[[length(reference_correlations) + 1L]] <- data.table(
      replicate = seq_len(1000L), population = population, method = method, correlation = vs_reference)
    name <- paste(population, method, sep = "_")
    correlation_matrices[[name]] <- pairwise
    fwrite(data.table(replicate = colnames(weights), as.data.table(pairwise)),
           boot_output(paste0("correlation_matrix_", name, ".csv.gz")), compress = "gzip", nThread = 4L)
  }
  cat("Computed exact Pearson and Spearman correlations:", population, "n =", sum(participant_rows), "\n")
}
saveRDS(correlation_matrices, boot_output("pairwise_correlation_matrices.rds"))
summary <- rbindlist(summaries)
fwrite(summary, boot_output("correlation_summary.csv"))
fwrite(rbindlist(reference_correlations), boot_output("correlations_with_no_bootstrap_reference.csv"))

pattern_summary <- rbindlist(lapply(seq_len(nrow(pattern_scores)), function(i) {
  values <- pattern_scores[i, seq_len(1000L)]
  qs <- quantile(values, c(.025, .5, .975), names = FALSE)
  data.table(mean_score = mean(values), sd_score = sd(values), p025_score = qs[1],
             median_score = qs[2], p975_score = qs[3], minimum_score = min(values),
             maximum_score = max(values), fraction_zero_score = mean(values == 0))
}))
individual_summary <- pattern_summary[pattern_ids]
individual_summary[, eid := target$eid]
individual_summary[, no_bootstrap_score := target$reference_scores]
setcolorder(individual_summary, c("eid", "no_bootstrap_score", setdiff(names(individual_summary), c("eid", "no_bootstrap_score"))))
fwrite(individual_summary, boot_output("individual_score_uncertainty.csv"))
fwrite(patterns, boot_output("pain_pattern_frequencies.csv"))
fwrite(data.table(population = c("all_proteomic", "with_any_pain"),
                  n = c(length(target$eid), sum(target$any_pain)),
                  n_patterns = c(nrow(patterns), sum(all_pain_patterns))), boot_output("scoring_population_counts.csv"))

coefficient_summary <- rbindlist(lapply(rownames(coefs$raw_beta), function(term) {
  values <- coefs$raw_beta[term, ]
  qs <- quantile(values, c(.025, .5, .975), names = FALSE)
  data.table(coefficient = term, reference_beta = target$reference_all_beta[term],
             mean_raw_beta = mean(values), sd_raw_beta = sd(values), p025_raw_beta = qs[1],
             median_raw_beta = qs[2], p975_raw_beta = qs[3])
}))
pain_weight_summary <- rbindlist(lapply(rownames(weights), function(term) {
  qs <- quantile(weights[term, ], c(.025, .5, .975), names = FALSE)
  data.table(coefficient = term, reference_normalized_weight = target$reference_weights[term],
             mean_normalized_weight = mean(weights[term, ]), sd_normalized_weight = sd(weights[term, ]),
             p025_normalized_weight = qs[1], median_normalized_weight = qs[2], p975_normalized_weight = qs[3],
             fraction_effective_weight_zero = mean(coefs$effective_beta[term, ] == 0))
}))
fwrite(coefficient_summary, boot_output("raw_beta_bootstrap_summary.csv"))
fwrite(pain_weight_summary, boot_output("pain_weight_bootstrap_summary.csv"))
add_check("each_replicate_normalized_weight_sum", max(abs(colSums(weights) - 1)), 1e-12)
add_check("pattern_vs_individual_score_calculation",
          max(abs(score_matrix[, c(1, 500, 1000)] - pattern_scores[pattern_ids, c(1, 500, 1000)])), 1e-12)
validation <- rbindlist(validation)
fwrite(validation, boot_output("correlation_calculation_checks.csv"))
if (!all(validation$passed)) stop("Correlation/score checks failed; inspect validation CSV.")
writeLines(capture.output(sessionInfo()), boot_output("sessionInfo_correlations.txt"))
print(summary)
cat("All correlation and scoring calculation checks passed.\n")
