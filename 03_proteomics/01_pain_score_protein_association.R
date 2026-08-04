source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

out_path <- result_path("pain_score_protein_association_results.csv")

covars <- make_covariates(UKB_DATA_FILE, n_pcs = 10)
scores <- data.table::fread(COMPOSITE_SCORE_FILE, select = c("eid", "pain_score"))
protein_cov <- data.table::fread(PROTEIN_COVARIATE_FILE)
data.table::setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- data.table::fread(PROTEIN_TIME_FILE)
time_protein_cols <- setdiff(names(protein_times), "eid")
data.table::setnames(protein_times, time_protein_cols, paste0(time_protein_cols, "_protein_time"))
proteins <- data.table::fread(PROTEIN_MATRIX_FILE)
protein_names <- setdiff(names(proteins), c("eid", "ins_index"))

analysis_data <- merge(scores, covars, by = "eid")
analysis_data <- merge(analysis_data, protein_cov, by = "eid")
analysis_data <- merge(analysis_data, protein_times, by = "eid")
analysis_data <- merge(analysis_data, proteins, by = "eid")

covariates <- c(MAIN_COVARIATES, "protein_cov_zscore", "season_final", "Batch")
covariate_terms <- paste(c(covariates, "protein_time"), collapse = " + ")

base_model <- stats::lm(as.formula(paste("pain_score ~", paste(covariates, collapse = " + "))),
                        data = analysis_data, na.action = na.exclude)

results <- vector("list", length(protein_names))
for (j in seq_along(protein_names)) {
  protein <- protein_names[j]
  time_col <- paste0(protein, "_protein_time")
  model_cols <- c("pain_score", covariates, protein, time_col)
  model_data <- as.data.frame(analysis_data[, model_cols, with = FALSE])
  names(model_data)[(ncol(model_data) - 1):ncol(model_data)] <- c("protein_value", "protein_time")
  fit <- try(stats::lm(as.formula(paste("pain_score ~ protein_value +", covariate_terms)), data = model_data),
             silent = TRUE)
  if (inherits(fit, "try-error") || !("protein_value" %in% rownames(summary(fit)$coefficients))) {
    results[[j]] <- data.frame(protein = protein, n = NA, beta = NA, se = NA, t = NA, p_value = NA)
  } else {
    coef_row <- summary(fit)$coefficients["protein_value", ]
    results[[j]] <- data.frame(
      protein = protein,
      n = nobs(fit),
      beta = coef_row["Estimate"],
      se = coef_row["Std. Error"],
      t = coef_row["t value"],
      p_value = coef_row["Pr(>|t|)"]
    )
  }
  if (j %% 250 == 0) cat("finished", j, "of", length(protein_names), "\n")
}

result <- data.table::rbindlist(results)
result[, fdr := p.adjust(p_value, method = "BH")]
result[, bonferroni := p.adjust(p_value, method = "bonferroni")]
data.table::setorder(result, p_value)
data.table::fwrite(result, out_path)

cat("base complete-case n:", nobs(base_model), "\n")
cat("proteins tested:", nrow(result), "\n")
cat("FDR<0.05:", sum(result$fdr < 0.05, na.rm = TRUE), "\n")
cat("Output:", out_path, "\n")

