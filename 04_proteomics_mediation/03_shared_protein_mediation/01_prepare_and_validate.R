source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
covariates <- c("age", "gender", "ethnicity", "deprivation", "income", "university_education",
                "smoking", "alcohol", "sedentary", "obesity", "prs", paste0("comp_", 1:20),
                "protein_cov_zscore", "season_final", "Batch", "protein_time")
factor_variables <- c("gender", "ethnicity", "deprivation", "income", "university_education",
                      "smoking", "alcohol", "sedentary", "obesity", "prs", "season_final", "Batch")
old_data_path <- file.path(REFERENCE_DIR, "overlap25_iorw_joint_mediation_analysis_data.csv")
dat <- as.data.frame(fread(old_data_path))
for (v in factor_variables) dat[[v]] <- factor(dat[[v]])
if (anyDuplicated(dat$eid) || nrow(dat) != 46591L || sum(dat$status) != 707L) stop("Unexpected original mediation cohort")
old_mediators <- grep("^aligned_.+_z$", names(dat), value = TRUE)
if (length(old_mediators) != 25L) stop("Expected 25 original aligned mediators")
old_proteins <- sub("_z$", "", sub("^aligned_", "", old_mediators))
new_shared_path <- file.path(SENS_DIR, "shared_proteins_fdr05_concordant.csv")
new_shared <- fread(new_shared_path)
if (nrow(new_shared) != 24L || !setequal(new_shared$protein, setdiff(old_proteins, "COMP"))) stop("Unexpected new shared set")
proteins <- old_proteins[old_proteins %in% new_shared$protein]
mediators <- paste0("aligned_", proteins, "_z")
if (any(sign(new_shared$beta_pain) != sign(new_shared$beta_PD))) stop("Shared directions disagree")

baseline_single <- fread(file.path(REFERENCE_DIR, "overlap25_iorw_per_protein_mediation_results.csv"))
baseline_joint <- fread(file.path(REFERENCE_DIR, "overlap25_iorw_joint_mediation_results.csv"))
validation <- list()
cat("Reproducing original observed mediation coefficients...\n"); flush.console()
for (protein in c("CHGA", "ITGAV")) {
  result <- legacy$compute_iorw_single(dat, covariates, paste0("aligned_", protein, "_z"), protein)
  reference_index <- which(baseline_single$protein == protein)
  reference <- baseline_single[reference_index]
  difference <- max(abs(unlist(result[, .(total_logHR, direct_logHR, indirect_logHR, proportion_mediated)]) -
                        unlist(reference[, .(total_logHR, direct_logHR, indirect_logHR, proportion_mediated)])))
  validation[[length(validation) + 1L]] <- data.table(check = paste0("original_single_", protein), max_difference = difference, passed = difference < 1e-8)
}
original_joint <- legacy$compute_iorw(dat, covariates, old_mediators)$result
difference <- max(abs(unlist(original_joint[, .(total_logHR, direct_logHR, indirect_logHR, proportion_mediated)]) -
                      unlist(baseline_joint[, .(total_logHR, direct_logHR, indirect_logHR, proportion_mediated)])))
validation[[length(validation) + 1L]] <- data.table(check = "original_joint25", max_difference = difference, passed = difference < 1e-8)
if (!all(rbindlist(validation)$passed)) stop("Original observed mediation does not reproduce")

new_score_path <- file.path(SENS_DIR, "composite_pain_scores.csv")
scores <- fread(new_score_path)
if (anyDuplicated(scores$eid)) stop("Duplicate new score IDs")
matched <- match(dat$eid, scores$eid)
if (anyNA(matched)) stop("New score missing for original mediation participants")
dat$pain_score_original <- dat$pain_score
dat$pain_score <- scores$pain_score[matched]
dat <- dat[, c("eid", "time", "status", "pain_score", "pain_score_original", covariates, mediators)]
if (anyNA(dat)) stop("Unexpected missingness in frozen new mediation data")
design <- make_design(dat, covariates, mediators)
covariate_audit <- data.table(covariate = covariates,
  n_observed_levels = vapply(covariates, function(v) length(unique(dat[[v]])), integer(1)),
  included_in_model = covariates %in% design$covariates)
fwrite(covariate_audit, file.path(OUT, "covariates_actually_used.csv"))
fwrite(as.data.table(dat), file.path(OUT, "mediation24_analysis_data.csv"))
observed <- list()
cat("Fitting observed 24-protein and joint models with the original functions...\n"); flush.console()
for (i in seq_along(proteins)) {
  fitted <- capture_fit(legacy$compute_iorw_single(dat, covariates, mediators[i], proteins[i]))
  if (nzchar(fitted$error)) stop(fitted$error)
  row <- fitted$value
  row[, warnings := paste(fitted$warnings, collapse = " | ")]
  observed[[i]] <- row
  cat("Observed", i, proteins[i], "complete\n"); flush.console()
}
observed <- rbindlist(observed)
fwrite(observed, file.path(OUT, "per_protein_observed.csv"))
joint <- legacy$compute_iorw(dat, covariates, mediators)$result
joint[, `:=`(method = "IORW_joint_24_direction_aligned_proteins", mediator_set = "24_direction_aligned_overlap_proteins")]
fwrite(joint, file.path(OUT, "joint24_observed.csv"))

cat("Validating coefficient-only bootstrap implementation against original formula fits...\n"); flush.console()
set.seed(20260602)
indices <- list(seq_len(nrow(dat)), sample.int(nrow(dat), nrow(dat), replace = TRUE))
benchmark <- list()
for (kind in c("CHGA", "ITGAV", "joint24")) {
  mi <- if (kind == "joint24") seq_along(mediators) else match(kind, proteins)
  for (j in seq_along(indices)) {
    idx <- indices[[j]]
    started <- proc.time()[3]
    fast <- fast_iorw(design, idx, mi, joint = kind == "joint24")
    elapsed <- proc.time()[3] - started
    reference <- if (kind == "joint24") {
      legacy$compute_iorw(dat[idx, ], covariates, mediators)$result
    } else legacy$compute_iorw_single(dat[idx, ], covariates, mediators[mi], kind)
    difference <- max(abs(fast - unlist(reference[, names(fast), with = FALSE])))
    validation[[length(validation) + 1L]] <- data.table(check = paste(kind, j, "fast_original_agreement"), max_difference = difference, passed = difference < 1e-8)
    benchmark[[length(benchmark) + 1L]] <- data.table(model = kind, resampled = j == 2L, fast_elapsed_seconds = elapsed)
  }
}
fwrite(rbindlist(validation), file.path(OUT, "reproduction_and_implementation_checks.csv"))
fwrite(rbindlist(benchmark), file.path(OUT, "benchmark.csv"))
if (!all(rbindlist(validation)$passed)) stop("Fast bootstrap does not reproduce original formula coefficients")
paths <- c(old_data_path, new_score_path, new_shared_path,
  file.path(SOURCE_DIR, c("analysis_overlap25_iorw_joint_mediation.R", "analysis_overlap25_iorw_per_protein_mediation.R")))
hashes <- data.table(path = paths, md5 = unname(tools::md5sum(paths)))
fwrite(hashes, file.path(OUT, "input_hashes.csv"))
saveRDS(list(data = dat, design = design, covariates = covariates, proteins = proteins,
             mediators = mediators, input_hashes = hashes, n_boot = 1000L, seed = 20260602L,
             note = "Original cohort, aligned protein values, covariate coding and IORW rules preserved; only pain score and COMP inclusion changed."),
        file.path(OUT, "prepared_inputs.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_preparation.txt"))
print(rbindlist(validation)); print(rbindlist(benchmark))
