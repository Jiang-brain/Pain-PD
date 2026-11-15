source("04_proteomics_mediation/01_independent_pain_score/00_common.R")
PRED_OUT <- sens_output("prediction")
dir.create(PRED_OUT, showWarnings = FALSE)
REF_DIR <- Sys.getenv("PD_PAIN_PREDICTION_SOURCE_DIR", file.path(SENS_DATA_DIR, "prediction_reference"))
SOURCE_DIR <- file.path(CODE_ROOT_SENS, "05_prediction", "source_reference")
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args)) args[1] else "baseline"
if (!phase %in% c("baseline", "new")) stop("Use baseline or new.")
is_function_definition <- function(expr) {
  nzchar(assigned_name(expr)) && is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))
}
load_functions <- function(path, environment) {
  expressions <- parse(path)
  for (expr in expressions) if (is_function_definition(expr)) eval(expr, envir = environment)
  expressions
}
case_dir <- function(name) {
  path <- file.path(PRED_OUT, name)
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  path
}

run_stable <- function(case, score_path, coef_path) {
  r <- new.env(parent = globalenv())
  expressions <- load_functions(file.path(SOURCE_DIR, "analysis_stable_nested_pfs_pd_cindex_base.R"), r)
  r$data_path <- sens_input("PD_data.csv")
  r$protein_path <- file.path(REF_DIR, "protein_new_missing60_filtered.csv")
  r$protein_cov_path <- sens_input("PD_data_matched_protein_covariates.csv")
  r$protein_time_path <- sens_input("PD_data_matched_protein_times.csv")
  r$score_path <- score_path
  r$coef_path <- coef_path
  r$out_path <- file.path(case_dir(case), "cindex.csv")
  r$out_summary_path <- file.path(case_dir(case), "cindex_summary.json")
  lhs <- vapply(expressions, assigned_name, "")
  start <- which(lhs == "data_censoring")[1]
  score_index <- which(lhs == "score_dt")[1]
  cache_path <- file.path(PRED_OUT, "stable_cox_base_data.rds")
  if (file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    r$analysis_data <- cached$data
    r$covariates <- cached$covariates
  } else {
    for (i in seq.int(start, score_index - 1L)) eval(expressions[[i]], envir = r)
    saveRDS(list(data = r$analysis_data, covariates = r$covariates), cache_path)
  }
  for (i in seq.int(score_index, length(expressions))) {
    eval(expressions[[i]], envir = r)
    if (lhs[i] == "coef_dt") r$coef_dt <- r$coef_dt[method == "elasticnet"]
  }
  r$write_json_simple(list(method = "Elastic Net only; original PA5 Cox specification, unchanged.",
                           score_path = score_path, coefficient_path = coef_path,
                           n_stable_proteins = nrow(r$stable_dt)), r$out_summary_path)
  covariate_audit <- data.table(covariate = r$covariates,
    unique_nonmissing = vapply(r$covariates, function(v) length(unique(r$model_data[[v]][!is.na(r$model_data[[v]])])), integer(1)),
    included_in_model = r$covariates %in% r$method_covariates)
  fwrite(covariate_audit, file.path(case_dir(case), "covariates_actually_used.csv"))
  load_functions(file.path(SOURCE_DIR, "analysis_stable_nested_pfs_pd_cindex.R"), r)
  # Reuse the identical fitted models instead of refitting them for NRI.
  nri <- rbindlist(list(
    r$continuous_nri_row("elasticnet", "traditional_covariates_vs_traditional_plus_stable_pfs",
                         r$base_fit, r$pfs_plus, r$model_data, max(r$model_data$time)),
    r$continuous_nri_row("elasticnet", "traditional_covariates_vs_traditional_plus_stable_proteins",
                         r$base_fit, r$protein_plus, r$model_data, max(r$model_data$time))))
  fwrite(nri, file.path(case_dir(case), "nri.csv"))
  list(cindex = r$result_dt, nri = nri)
}

run_shared <- function() {
  r <- new.env(parent = globalenv())
  expressions <- load_functions(file.path(SOURCE_DIR, "analysis_overlap25_pca_prediction_cindex.R"), r)
  r$data_path <- sens_input("PD_data.csv")
  r$intersection_path <- file.path(REF_DIR, "protein_intersection_FDR05.csv")
  r$protein_cov_path <- sens_input("PD_data_matched_protein_covariates.csv")
  r$protein_time_path <- sens_input("PD_data_matched_protein_times.csv")
  r$followup_2y <- 365.25 * 2
  lhs <- vapply(expressions, assigned_name, "")
  start <- which(lhs == "same")[1]
  end <- which(lhs == "covariates")[1]
  for (i in seq.int(start, end)) {
    if (lhs[i] == "pca_scores") {
      # Read ONLY stored individual protein columns and IDs. No PCA values,
      # PCA scores, loadings, or PCA models enter any prediction here.
      r$pca_scores <- fread(file.path(REF_DIR, "overlap25_directional_pca_signature.csv"),
                           select = c("eid", r$protein_cols))
    } else eval(expressions[[i]], envir = r)
  }
  old_proteins <- r$overlap_proteins
  new_proteins <- fread(sens_output("shared_proteins_fdr05_concordant.csv"))$protein
  if (!all(new_proteins %in% old_proteins)) stop("Unexpected new shared proteins; rebuild aligned columns explicitly.")
  task_specs <- list(
    list(name = "incident_PD_within_5y", data = r$make_task(r$analysis_data, 5)),
    list(name = "incident_PD_within_10y", data = r$make_task(r$analysis_data, 10)),
    list(name = "incident_PD_full_followup_after_2y", data = r$make_full_followup_task(r$analysis_data)))
  rows <- list()
  nri_rows <- list()
  for (task in task_specs) {
    # Freeze the original complete-case cohort before dropping COMP, so the
    # comparison isolates the panel change rather than changing participants.
    columns <- c("time_horizon", "status_horizon", r$covariates, r$protein_cols)
    model_data <- as.data.frame(task$data[, columns, with = FALSE])
    model_data <- droplevels(model_data[complete.cases(model_data), ])
    covs <- r$usable_covariates(model_data, r$covariates)
    cov_terms <- paste(covs, collapse = "+")
    base_fit <- coxph(as.formula(paste("Surv(time_horizon,status_horizon)~", cov_terms)), data = model_data)
    for (case in c("shared_original25", "shared_new24")) {
      proteins <- if (case == "shared_original25") old_proteins else new_proteins
      terms <- paste(paste0("aligned_", proteins, "_z"), collapse = "+")
      alone <- coxph(as.formula(paste("Surv(time_horizon,status_horizon)~", terms)), data = model_data)
      combined <- coxph(as.formula(paste("Surv(time_horizon,status_horizon)~", cov_terms, "+", terms)), data = model_data)
      rows[[length(rows) + 1L]] <- cbind(case = case, rbindlist(list(
        r$fit_row(task$name, "traditional_covariates", base_fit, NULL, 0, "traditional_covariates"),
        r$fit_row(task$name, "shared_proteins_alone", alone, NULL, length(proteins), "shared_proteins"),
        r$fit_row(task$name, "traditional_plus_shared_proteins", combined, base_fit, length(proteins), "shared_proteins"))))
      nri_row <- r$continuous_nri_row(task$name, "traditional_vs_plus_shared", base_fit, combined,
                                     model_data, max(model_data$time_horizon))
      nri_row$method <- paste0("continuous NRI for ", length(proteins),
                               " direction-aligned proteins; SE from reclassification proportions")
      nri_rows[[length(nri_rows) + 1L]] <- cbind(case = case, nri_row)
    }
  }
  cindex <- rbindlist(rows)
  nri <- rbindlist(nri_rows)
  fwrite(cindex, file.path(PRED_OUT, "shared_prediction_cindex.csv"))
  fwrite(nri, file.path(PRED_OUT, "shared_prediction_nri.csv"))
  reference <- fread(file.path(REF_DIR, "overlap25_pca_prediction_cindex_5y_10y_full_followup.csv"))
  reference <- reference[model %in% c("traditional_covariates", "overlap25_proteins_alone", "traditional_plus_overlap25_proteins")]
  reference[, model := sub("overlap25", "shared", model)]
  comparison <- merge(cindex[case == "shared_original25"], reference, by = c("task", "model"), suffixes = c("_rerun", "_saved"))
  comparison[, passed := abs(c_index_rerun - c_index_saved) < 1e-10 & n_rerun == n_saved & events_rerun == events_saved]
  fwrite(comparison, file.path(PRED_OUT, "original_shared_cox_reproduction_checks.csv"))
  if (nrow(comparison) != 9L || !all(comparison$passed)) stop("Original shared Cox baseline does not reproduce.")
  print(cindex)
  print(nri)
}

if (phase == "baseline") {
  baseline <- run_stable("stable_original374", file.path(REF_DIR, "pain_related_proteomic_score_stable_nested_pfs.csv"),
                          file.path(REF_DIR, "pain_related_proteomic_score_stable_nested_coefficients.csv"))
  saved <- fread(file.path(REF_DIR, "pd_incidence_stable_nested_pfs_cindex.csv"))[method == "elasticnet"]
  check <- merge(baseline$cindex, saved, by = c("method", "model_type"), suffixes = c("_rerun", "_saved"))
  check[, passed := abs(c_index_rerun - c_index_saved) < 1e-10 & n_rerun == n_saved & events_rerun == events_saved]
  fwrite(check, file.path(PRED_OUT, "original_stable_cox_reproduction_checks.csv"))
  if (nrow(check) != 5L || !all(check$passed)) stop("Original stable Cox baseline does not reproduce.")
  nri_saved <- fread(file.path(REF_DIR, "pd_incidence_stable_nested_pfs_nri.csv"))[method == "elasticnet"]
  nri_check <- merge(baseline$nri, nri_saved, by = c("method", "comparison"), suffixes = c("_rerun", "_saved"))
  nri_check[, passed := abs(nri_rerun - nri_saved) < 1e-10]
  fwrite(nri_check, file.path(PRED_OUT, "original_stable_nri_reproduction_checks.csv"))
  if (!all(nri_check$passed)) stop("Original stable NRI baseline does not reproduce.")
  run_shared()
} else {
  run_stable("stable_new_score", file.path(PRED_OUT, "new_stable_elastic_net_pfs.csv"),
              file.path(PRED_OUT, "new_stable_elastic_net_coefficients.csv"))
}
cat("Cox phase", phase, "completed.\n")
