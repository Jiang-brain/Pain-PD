source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
inputs <- readRDS(file.path(OUT, "prepared_inputs.rds"))
jobs <- readRDS(file.path(OUT, "bootstrap_jobs.rds"))
paths <- list.files(file.path(OUT, "bootstrap_checkpoints"), pattern = "\\.rds$", full.names = TRUE)
if (length(paths) != length(jobs)) stop("Not all bootstrap jobs have completed")
checkpoints <- lapply(paths, readRDS)
boot <- rbindlist(lapply(checkpoints, `[[`, "rows"), fill = TRUE)
if (anyDuplicated(boot[, .(model, bootstrap)])) stop("Duplicate model/bootstrap rows")
counts <- boot[, .(attempts = .N, valid = sum(valid), warnings = sum(nzchar(warnings)), errors = sum(nzchar(error))), by = model]
if (nrow(counts) != 25L || any(counts$attempts != 1000L)) stop("Unexpected bootstrap counts")
fwrite(boot, file.path(OUT, "mediation24_bootstrap.csv"))
fwrite(counts, file.path(OUT, "bootstrap_counts.csv"))
fwrite(boot[nzchar(warnings), .(n_fits = .N), by = warnings][order(-n_fits)],
       file.path(OUT, "bootstrap_warning_summary.csv"))
observed <- fread(file.path(OUT, "per_protein_observed.csv"))
single_boot <- copy(boot[model != "joint24"])
single_boot[, `:=`(protein = model, mediator = paste0("aligned_", model, "_z"))]
result <- legacy$add_bootstrap_results(observed, single_boot)
for (i in seq_len(nrow(result))) {
  protein_name <- result$protein[i]
  b <- single_boot[protein == protein_name & valid == TRUE]
  result[i, `:=`(
    proportion_lower = quantile(b$proportion_mediated, 0.025, na.rm = TRUE),
    proportion_upper = quantile(b$proportion_mediated, 0.975, na.rm = TRUE),
    n_boot_warnings = sum(nzchar(b$warnings)))]
}
result[, indirect_p_boot_two_sided_fdr := p.adjust(indirect_p_boot_two_sided, "BH")]
# Primary reporting uses the indirect log-HR percentile CI, not legacy tail FDR.
result[, indirect_95CI_excludes_zero := indirect_lower > 0 | indirect_upper < 0]
result[, proportion_95CI_excludes_zero := proportion_lower > 0 | proportion_upper < 0]
setorder(result, proportion_mediated)
fwrite(result, file.path(OUT, "per_protein_mediation_results.csv"))
joint <- fread(file.path(OUT, "joint24_observed.csv"))
joint_boot <- boot[model == "joint24"]
joint <- legacy$add_boot_ci(joint, joint_boot)
joint[, indirect_95CI_excludes_zero := indirect_lower > 0 | indirect_upper < 0]
joint[, proportion_95CI_excludes_zero := proportion_lower > 0 | proportion_upper < 0]
fwrite(joint, file.path(OUT, "joint24_mediation_results.csv"))

if (!all(unname(tools::md5sum(inputs$input_hashes$path)) == inputs$input_hashes$md5)) stop("Original inputs changed")
checks <- rbindlist(list(
  data.table(check = "original_cohort_and_events", passed = nrow(inputs$data) == 46591L && sum(inputs$data$status) == 707L),
  data.table(check = "all_24_single_and_joint_bootstraps_attempted", passed = all(counts$attempts == 1000L)),
  data.table(check = "all_bootstrap_effect_estimates_finite", passed = all(counts$valid == 1000L)),
  data.table(check = "original_inputs_unchanged", passed = TRUE),
  data.table(check = "24_mediators_excluding_COMP", passed = length(inputs$proteins) == 24L && !"COMP" %in% inputs$proteins),
  data.table(check = "all_CIs_finite", passed = all(is.finite(result$proportion_lower)) && all(is.finite(result$proportion_upper)) && is.finite(joint$proportion_lower) && is.finite(joint$proportion_upper))
))
fwrite(checks, file.path(OUT, "mediation_validation_checks.csv"))
summary <- c(
  "# IORW Mediation Using the Independently Weighted Pain Score",
  "",
  sprintf("Fixed analysis cohort: %s participants; %s incident PD events.", nrow(inputs$data), sum(inputs$data$status)),
  "The new non-proteomic-derived pain score replaces the original score. COMP is excluded from the shared set; all other participants, aligned standardized protein values, covariate codings and IORW rules are retained.",
  "",
  sprintf("Primary reporting: %s of 24 single-protein indirect log-HR 95%% bootstrap percentile CIs exclude zero (unadjusted for multiplicity).", sum(result$indirect_95CI_excludes_zero)),
  sprintf("Legacy comparison only, not the final criterion: one-sided bootstrap BH FDR < 0.05, %s of 24; two-sided bootstrap BH FDR < 0.05, %s of 24.", sum(result$indirect_p_boot_fdr < 0.05), sum(result$indirect_p_boot_two_sided_fdr < 0.05)),
  sprintf("Observed single-protein proportions mediated: %.4f%% to %.4f%%.", min(result$proportion_mediated) * 100, max(result$proportion_mediated) * 100),
  sprintf("Joint 24-protein proportion mediated: %.4f%% (95%% bootstrap percentile CI %.4f%% to %.4f%%).", joint$proportion_mediated * 100, joint$proportion_lower * 100, joint$proportion_upper * 100),
  sprintf("Bootstrap indirect-effect two-sided tail proportion for joint set: %g.", joint$indirect_p_boot),
  sprintf("Fits with finite effect estimates: %s/%s. Fits with Cox warnings: %s (retained, not silently discarded).", sum(counts$valid), sum(counts$attempts), sum(counts$warnings)),
  "",
  "## Interpretation and Preserved Limitations",
  "- These are observational mediation estimates, not causal or intervention effects.",
  "- The pain score is fixed during mediation bootstrap; this bootstrap resamples mediation participants, not pain-score derivation participants.",
  "- Each of the 24 single-mediator models and the joint model requests 1,000 participant-level bootstrap resamples, preserving the original serial RNG stream (seed 20260602).",
  "- Original Gaussian density-ratio weights, 99th-percentile trimming, Efron Cox ties, and median-imputed/aligned standardized proteins are retained.",
  "- The original deprivation/PRS recoding produces constant factors that are excluded by usable_covariates. This known limitation is not silently repaired; see covariates_actually_used.csv.",
  "- Final reporting uses whether the indirect log-HR 95% bootstrap percentile CI excludes zero (or indirect HR CI excludes one). These are unadjusted intervals, not FDR-corrected significance. Original one-sided/two-sided tail and BH columns are preserved only for historical comparison.",
  "- A bootstrap tail proportion of zero is finite-resolution evidence, not an exact P value of zero; do not report P=0.",
  "- Many bootstrap fits can warn that nuisance coefficients, especially batch contrasts, may be infinite. Representative warned resamples are checked against the original formula fits, which show the same issue. Finite exposure/indirect estimates are retained under the original implementation, but these warnings are a model-stability limitation, not evidence of fully converged fits.",
  "- Checkpointed coefficient-only fits use the same lm.fit and coxph.fit engines. Every duplicated event row is retained, preserving Efron ties; duplicated censored rows may be collapsed into mathematically equivalent multiplicity weights. Both implementations are checked for numerical agreement. Robust SE/concordance/naive-adjusted fits are not recomputed in each bootstrap because percentile CIs require only the effect estimates.",
  "- This mediation script does not rerun PPI/enrichment. Figure 4G uses the separately supplied 24-protein enrichment exports; see 07_plot_figure4G_top5.py."
)
writeLines(summary, file.path(OUT, "MEDIATION_RESULTS.md"))
print(result[, .(protein, proportion_mediated, proportion_lower, proportion_upper, indirect_lower, indirect_upper, indirect_95CI_excludes_zero)])
print(joint[, .(proportion_mediated, proportion_lower, proportion_upper, indirect_p_boot, n_boot_success)])
if (!all(checks$passed)) stop("Mediation validation failed; inspect result diagnostics")
