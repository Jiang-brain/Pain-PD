source("04_proteomics_mediation/01_independent_pain_score/00_common.R")

# Evaluate the original PA model code unchanged, replacing only input/output
# paths. It models pain_score ~ protein + covariates, including 20 genetic PCs.
analysis <- new.env(parent = globalenv())
analysis$data_path <- sens_input("PD_data.csv")
analysis$score_path <- sens_output("composite_pain_scores.csv")
analysis$protein_path <- sens_input("protein_new_missing60_filtered.csv")
analysis$protein_cov_path <- sens_input("PD_data_matched_protein_covariates.csv")
analysis$protein_time_path <- sens_input("PD_data_matched_protein_times.csv")
analysis$out_path <- sens_output("pain_score_protein_association_results.csv")
if (!file.exists(analysis$score_path)) stop("Run sensitivity step 01 first.")
path_variables <- c("data_path", "score_path", "protein_path", "protein_cov_path",
                    "protein_time_path", "out_path")
for (expr in parse(reference_path("analysis_pain_score_protein_association.R"))) {
  if (assigned_name(expr) %in% path_variables || is_libpaths(expr)) next
  eval(expr, envir = analysis)
}
result <- analysis$result
fwrite(result[fdr < 0.05], sens_output("pain_score_protein_significant_FDR05.csv"))
fwrite(result[bonferroni < 0.05], sens_output("pain_score_protein_significant_Bonferroni05.csv"))
fwrite(data.table(base_complete_case_n = analysis$base_n, proteins_tested = nrow(result),
                  failed_models = sum(!is.finite(result$p_value)),
                  nominal_significant = sum(result$p_value < 0.05, na.rm = TRUE),
                  fdr_significant = sum(result$fdr < 0.05, na.rm = TRUE),
                  bonferroni_significant = sum(result$bonferroni < 0.05, na.rm = TRUE)),
       sens_output("protein_association_counts.csv"))
original_result <- fread(sens_input("pain_score_protein_association_results.csv"))
checks <- merge(result[, .(protein, n_new = n)], original_result[, .(protein, n_original = n)],
                by = "protein", all = TRUE)
fwrite(checks, sens_output("protein_association_sample_size_check.csv"))
if (anyNA(checks) || any(checks$n_new != checks$n_original)) {
  stop("Protein-specific sample sizes changed; inspect sample-size check before interpretation.")
}
writeLines(capture.output(sessionInfo()), sens_output("sessionInfo_association.txt"))
write_provenance(reference_path("analysis_pain_score_protein_association.R"), "association_source_hashes.csv")
