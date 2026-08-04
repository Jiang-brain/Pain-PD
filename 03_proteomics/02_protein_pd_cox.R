source("00_setup/ukb_helpers.R")
load_packages(c("data.table", "survival"))

out_path <- result_path("protein_PD_cox_results.csv")

outcome_data <- make_outcome_data(UKB_DATA_FILE)
covars <- make_covariates(UKB_DATA_FILE, n_pcs = 10)
analysis_data <- merge(outcome_data, covars, by = "eid")
analysis_data <- apply_main_exclusions(analysis_data)

proteins <- data.table::fread(PROTEIN_MATRIX_FILE)
protein_cov <- data.table::fread(PROTEIN_COVARIATE_FILE)
data.table::setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- data.table::fread(PROTEIN_TIME_FILE)
time_protein_cols <- setdiff(names(protein_times), "eid")
data.table::setnames(protein_times, time_protein_cols, paste0(time_protein_cols, "_protein_time"))
protein_names <- setdiff(names(proteins), c("eid", "ins_index"))

analysis_data <- merge(as.data.table(analysis_data), protein_cov, by = "eid")
analysis_data <- merge(analysis_data, protein_times, by = "eid")
analysis_data <- merge(analysis_data, proteins, by = "eid")

covariates <- c(MAIN_COVARIATES, "protein_cov_zscore", "season_final", "Batch")
covariate_terms <- paste(c(covariates, "protein_time"), collapse = " + ")

results <- vector("list", length(protein_names))
for (j in seq_along(protein_names)) {
  protein <- protein_names[j]
  time_col <- paste0(protein, "_protein_time")
  model_cols <- c("time", "status", covariates, protein, time_col)
  model_data <- as.data.frame(analysis_data[, model_cols, with = FALSE])
  names(model_data)[(ncol(model_data) - 1):ncol(model_data)] <- c("protein_value", "protein_time")
  fit <- try(survival::coxph(
    as.formula(paste("Surv(time, status) ~ protein_value +", covariate_terms)),
    data = model_data
  ), silent = TRUE)
  if (inherits(fit, "try-error") || !("protein_value" %in% rownames(summary(fit)$coefficients))) {
    results[[j]] <- data.frame(protein = protein, n = NA, events = NA, beta = NA, HR = NA, lower = NA, upper = NA, p_value = NA)
  } else {
    fit_summary <- summary(fit)
    coef_row <- fit_summary$coefficients["protein_value", ]
    conf_row <- fit_summary$conf.int["protein_value", ]
    results[[j]] <- data.frame(
      protein = protein,
      n = fit$n,
      events = fit$nevent,
      beta = coef_row["coef"],
      HR = conf_row["exp(coef)"],
      lower = conf_row["lower .95"],
      upper = conf_row["upper .95"],
      p_value = coef_row["Pr(>|z|)"]
    )
  }
  if (j %% 250 == 0) cat("finished", j, "of", length(protein_names), "\n")
}

result <- data.table::rbindlist(results)
result[, fdr := p.adjust(p_value, method = "BH")]
result[, bonferroni := p.adjust(p_value, method = "bonferroni")]
data.table::setorder(result, p_value)
data.table::fwrite(result, out_path)

cat("eligible protein-merged rows:", nrow(analysis_data), "\n")
cat("proteins tested:", nrow(result), "\n")
cat("FDR<0.05:", sum(result$fdr < 0.05, na.rm = TRUE), "\n")
cat("Output:", out_path, "\n")

