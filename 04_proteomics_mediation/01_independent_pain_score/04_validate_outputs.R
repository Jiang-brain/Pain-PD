source("04_proteomics_mediation/01_independent_pain_score/00_common.R")

checks <- list()
check <- function(name, condition, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(check = name, passed = isTRUE(condition), detail = detail)
}
scores <- fread(sens_output("composite_pain_scores.csv"))
train_ids <- fread(sens_output("derivation_participant_ids.csv"))
assert_unique_eid(scores, "Scored participants")
assert_unique_eid(train_ids, "Derivation participants")
roster_paths <- fread(sens_output("protein_roster_sources.csv"))$path
roster_ids <- unique(unlist(lapply(roster_paths, function(path) fread(path, select = "eid")$eid)))
check("no_derivation_protein_roster_overlap", !any(train_ids$eid %in% roster_ids))
check("no_derivation_scored_overlap", !any(train_ids$eid %in% scores$eid))
check("all_scored_have_proteins", all(scores$eid %in% roster_ids))
reproduction <- fread(sens_output("baseline_score_reproduction_check.csv"))
check("baseline_scores_reproduced", reproduction$passed && reproduction$max_absolute_difference < 1e-12)
frozen <- readRDS(sens_output("fixed_weights.rds"))
check("fixed_denominator", identical(frozen$denominator, sum(frozen$weights$effective_beta)))
check("negative_weights_clipped", all(frozen$weights$effective_beta == pmax(frozen$weights$raw_beta, 0)))

# Independent scoring check: reconstruct the 14 indicators directly from raw
# duration fields, without using the cohort/model matrix from step 01.
site_fields <- c(headache = "3799-0.0", facial = "4067-0.0", neck_shoulder = "3404-0.0",
                 back = "3571-0.0", stomach = "3741-0.0", hip_pain = "3414-0.0", knee = "3773-0.0")
raw <- fread(sens_input("PD_data.csv"), select = c("eid", "6159-0.0", unname(site_fields)))
raw <- raw[match(scores$eid, raw$eid)]
indicators <- matrix(0, nrow = nrow(raw), ncol = nrow(frozen$weights))
for (j in seq_len(nrow(frozen$weights))) {
  term <- frozen$weights$term[j]
  site <- sub("_(acute|chronic)$", "", term)
  duration_code <- if (grepl("_acute$", term)) 0 else 1
  values <- raw[[site_fields[[site]]]]
  indicators[, j] <- as.integer(!is.na(values) & values == duration_code)
}
indicators[which(raw[["6159-0.0"]] == -7), ] <- 0
independent_scores <- as.numeric(indicators %*% frozen$weights$effective_beta) / frozen$denominator
score_difference <- max(abs(independent_scores - scores$pain_score))
check("independent_fixed_beta_score_recalculation", is.finite(score_difference) && score_difference < 1e-12,
      paste("max_absolute_difference", format(score_difference, digits = 16)))
check("score_range_zero_to_one", all(scores$pain_score >= 0 & scores$pain_score <= 1 + 1e-12))

pain <- fread(sens_output("pain_score_protein_association_results.csv"))
baseline_pain <- fread(sens_input("pain_score_protein_association_results.csv"))
sizes <- merge(pain[, .(protein, n_new = n)], baseline_pain[, .(protein, n_original = n)],
               by = "protein", all = TRUE)
check("same_protein_universe_and_sample_sizes", !anyNA(sizes) && all(sizes$n_new == sizes$n_original))
check("all_protein_models_finite", all(is.finite(pain$p_value)) && all(is.finite(pain$beta)))
check("BH_FDR_recalculation", isTRUE(all.equal(pain$fdr, p.adjust(pain$p_value, method = "BH"), tolerance = 1e-12)))
check("Bonferroni_recalculation", isTRUE(all.equal(pain$bonferroni, p.adjust(pain$p_value, method = "bonferroni"), tolerance = 1e-12)))
shared <- fread(sens_output("shared_proteins_fdr05_concordant.csv"))
check("shared_significance_and_direction", all(shared$fdr_pain < 0.05 & shared$fdr_PD < 0.05 &
                                               shared$beta_pain * shared$beta_PD > 0))
comparison <- fread(sens_output("shared_protein_membership_comparison.csv"))
check("original_25_shared_reproduced", sum(comparison$original_shared) == 25L)
check("shared_membership_output_consistency", setequal(shared$protein, comparison[sensitivity_shared == TRUE, protein]))
hashes <- fread(sens_output("unchanged_PD_result_hash.csv"))
check("PD_results_unchanged", identical(hashes$md5, unname(tools::md5sum(hashes$path))))
result <- rbindlist(checks)
fwrite(result, sens_output("validation_checks.csv"))
writeLines(capture.output(sessionInfo()), sens_output("sessionInfo_validation.txt"))
print(result)
if (!all(result$passed)) stop("One or more validation checks failed.")
cat("All", nrow(result), "validation checks passed.\n")
