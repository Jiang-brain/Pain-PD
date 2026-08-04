source("00_setup/config.R")

joint_script <- "04_mediation/01_overlap25_iorw_joint_mediation.R"

load_joint_definitions <- function(path) {
  if (!file.exists(path)) stop("Cannot find joint IORW script: ", path)
  lines <- readLines(path, warn = FALSE)
  main_start <- grep("^covariates\\s*<-", lines)
  if (length(main_start) == 0) stop("Cannot locate main analysis block in: ", path)
  eval(parse(text = paste(lines[seq_len(main_start[1] - 1)], collapse = "\n")), envir = .GlobalEnv)
}

load_joint_definitions(joint_script)

out_data_path <- result_path("overlap25_iorw_per_protein_mediation_analysis_data.csv")
out_result_path <- result_path("overlap25_iorw_per_protein_mediation_results.csv")
out_boot_path <- result_path("overlap25_iorw_per_protein_mediation_bootstrap.csv")
out_summary_path <- result_path("overlap25_iorw_per_protein_mediation_summary.json")

set.seed(20260602)
n_boot <- as.integer(Sys.getenv("IORW_N_BOOT", "1000"))
weight_trim_quantile <- as.numeric(Sys.getenv("IORW_WEIGHT_TRIM_Q", "0.99"))

bt <- function(x) paste0("`", gsub("`", "", x), "`")

rhs_terms <- function(vars) {
  vars <- vars[!is.na(vars) & nzchar(vars)]
  if (length(vars) == 0) return("1")
  paste(bt(vars), collapse = " + ")
}

make_lm_formula <- function(outcome, vars) {
  as.formula(paste(bt(outcome), "~", rhs_terms(vars)))
}

make_cox_formula <- function(vars) {
  as.formula(paste("Surv(time, status) ~", rhs_terms(vars)))
}

extract_coef <- function(fit, term) {
  cm <- summary(fit)$coefficients
  rn <- gsub("`", "", rownames(cm), fixed = TRUE)
  idx <- which(rn == term)
  if (length(idx) == 0) stop("Term not found in model: ", term)
  cm[idx[1], ]
}

coef_value <- function(coef_row, name) {
  if (name %in% names(coef_row)) return(as.numeric(coef_row[name]))
  NA_real_
}

compute_iorw_single <- function(dat, covariates, mediator_col, protein_name,
                                weight_trim_quantile = 0.99) {
  dat <- droplevels(as.data.frame(dat))
  covariates <- usable_covariates(dat, covariates)

  numerator_formula <- make_lm_formula("pain_score", covariates)
  denominator_formula <- make_lm_formula("pain_score", c(covariates, mediator_col))
  total_formula <- make_cox_formula(c("pain_score", covariates))
  naive_direct_formula <- make_cox_formula(c("pain_score", mediator_col, covariates))

  numerator_fit <- lm(numerator_formula, data = dat)
  denominator_fit <- lm(denominator_formula, data = dat)
  numerator_sd <- sqrt(sum(residuals(numerator_fit)^2) / numerator_fit$df.residual)
  denominator_sd <- sqrt(sum(residuals(denominator_fit)^2) / denominator_fit$df.residual)

  numerator_density <- dnorm(dat$pain_score, mean = predict(numerator_fit, newdata = dat), sd = numerator_sd)
  denominator_density <- dnorm(dat$pain_score, mean = predict(denominator_fit, newdata = dat), sd = denominator_sd)
  numerator_density <- pmax(numerator_density, .Machine$double.xmin)
  denominator_density <- pmax(denominator_density, .Machine$double.xmin)

  weights <- numerator_density / denominator_density
  weights[!is.finite(weights)] <- NA_real_
  trim_value <- quantile(weights, probs = weight_trim_quantile, na.rm = TRUE)
  weights_trimmed <- pmin(weights, trim_value)
  weights_trimmed <- weights_trimmed / mean(weights_trimmed, na.rm = TRUE)
  dat$iorw_weight <- weights_trimmed

  total_fit <- coxph(total_formula, data = dat)
  direct_fit <- coxph(total_formula, data = dat, weights = iorw_weight, robust = TRUE)
  naive_direct_fit <- coxph(naive_direct_formula, data = dat)

  total_coef <- extract_coef(total_fit, "pain_score")
  direct_coef <- extract_coef(direct_fit, "pain_score")
  naive_direct_coef <- extract_coef(naive_direct_fit, "pain_score")
  naive_mediator_coef <- extract_coef(naive_direct_fit, mediator_col)
  denominator_coef <- extract_coef(denominator_fit, mediator_col)

  total <- coef_value(total_coef, "coef")
  direct <- coef_value(direct_coef, "coef")
  indirect <- total - direct

  data.table(
    method = "IORW_single_direction_aligned_protein",
    protein = protein_name,
    mediator = mediator_col,
    n = nrow(dat),
    events = sum(dat$status),
    exposure = "pain_score",
    total_logHR = total,
    total_HR = exp(total),
    total_se = coef_value(total_coef, "se(coef)"),
    total_p_value = coef_value(total_coef, "Pr(>|z|)"),
    direct_logHR = direct,
    direct_HR = exp(direct),
    direct_robust_se = coef_value(direct_coef, "robust se"),
    direct_p_value = coef_value(direct_coef, "Pr(>|z|)"),
    indirect_logHR = indirect,
    indirect_HR = exp(indirect),
    proportion_mediated = indirect / total,
    pain_on_mediator_beta = coef_value(denominator_coef, "Estimate"),
    pain_on_mediator_p_value = coef_value(denominator_coef, "Pr(>|t|)"),
    naive_direct_logHR = coef_value(naive_direct_coef, "coef"),
    naive_direct_HR = exp(coef_value(naive_direct_coef, "coef")),
    naive_mediator_logHR = coef_value(naive_mediator_coef, "coef"),
    naive_mediator_HR = exp(coef_value(naive_mediator_coef, "coef")),
    naive_mediator_p_value = coef_value(naive_mediator_coef, "Pr(>|z|)"),
    weight_mean = mean(weights_trimmed, na.rm = TRUE),
    weight_sd = sd(weights_trimmed, na.rm = TRUE),
    weight_min = min(weights_trimmed, na.rm = TRUE),
    weight_p50 = quantile(weights_trimmed, 0.50, na.rm = TRUE),
    weight_p95 = quantile(weights_trimmed, 0.95, na.rm = TRUE),
    weight_p99 = quantile(weights_trimmed, 0.99, na.rm = TRUE),
    weight_max = max(weights_trimmed, na.rm = TRUE),
    weight_effective_n = sum(weights_trimmed, na.rm = TRUE)^2 / sum(weights_trimmed^2, na.rm = TRUE),
    weight_trim_quantile = weight_trim_quantile,
    weight_trim_value_raw = trim_value,
    numerator_R2 = summary(numerator_fit)$r.squared,
    denominator_R2 = summary(denominator_fit)$r.squared,
    c_index_total = unname(summary(total_fit)$concordance["C"]),
    c_index_direct = unname(summary(direct_fit)$concordance["C"])
  )
}

run_observed_per_protein <- function(dat, covariates, mediator_cols, proteins, weight_trim_quantile) {
  rows <- vector("list", length(mediator_cols))
  for (m in seq_along(mediator_cols)) {
    cat("observed IORW", m, "of", length(mediator_cols), proteins[m], "\n")
    flush.console()
    rows[[m]] <- compute_iorw_single(
      dat = dat,
      covariates = covariates,
      mediator_col = mediator_cols[m],
      protein_name = proteins[m],
      weight_trim_quantile = weight_trim_quantile
    )
  }
  rbindlist(rows, fill = TRUE)
}

bootstrap_per_protein <- function(dat, covariates, mediator_cols, proteins,
                                  n_boot, weight_trim_quantile) {
  if (n_boot <= 0) return(data.table())
  n <- nrow(dat)
  rows <- vector("list", length(mediator_cols) * n_boot)
  k <- 1L
  for (m in seq_along(mediator_cols)) {
    cat("bootstrap protein", m, "of", length(mediator_cols), proteins[m], "\n")
    flush.console()
    for (b in seq_len(n_boot)) {
      idx <- sample.int(n, size = n, replace = TRUE)
      fit <- try(
        compute_iorw_single(dat[idx, ], covariates, mediator_cols[m], proteins[m], weight_trim_quantile),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) {
        rows[[k]] <- data.table(protein = proteins[m], mediator = mediator_cols[m], boot = b, error = as.character(fit))
      } else {
        fit[, boot := b]
        rows[[k]] <- fit
      }
      k <- k + 1L
      if (b %% 50 == 0) {
        cat("  bootstrap", b, "of", n_boot, "\n")
        flush.console()
      }
    }
  }
  rbindlist(rows, fill = TRUE)
}

boot_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(list(lower = NA_real_, upper = NA_real_, p_lt0 = NA_real_,
                p_gt0 = NA_real_, p_two_sided = NA_real_, n = 0L))
  }
  p_lt0 <- mean(x < 0)
  p_gt0 <- mean(x > 0)
  list(
    lower = unname(quantile(x, 0.025, na.rm = TRUE)),
    upper = unname(quantile(x, 0.975, na.rm = TRUE)),
    p_lt0 = p_lt0,
    p_gt0 = p_gt0,
    p_two_sided = min(1, 2 * min(p_lt0, p_gt0)),
    n = length(x)
  )
}

add_bootstrap_results <- function(result, boot_dt) {
  result <- copy(result)
  result[, `:=`(
    indirect_lower = NA_real_,
    indirect_upper = NA_real_,
    indirect_HR_lower = NA_real_,
    indirect_HR_upper = NA_real_,
    indirect_p_boot = NA_real_,
    indirect_p_boot_lt0 = NA_real_,
    indirect_p_boot_gt0 = NA_real_,
    indirect_p_boot_two_sided = NA_real_,
    n_boot_success = 0L
  )]

  if (nrow(boot_dt) == 0 || !("indirect_logHR" %in% names(boot_dt))) {
    result[, indirect_p_boot_fdr := NA_real_]
    return(result)
  }

  for (i in seq_len(nrow(result))) {
    valid <- boot_dt[
      protein == result$protein[i] &
        mediator == result$mediator[i] &
        is.finite(indirect_logHR)
    ]
    bs <- boot_summary(valid$indirect_logHR)
    result[i, `:=`(
      indirect_lower = bs$lower,
      indirect_upper = bs$upper,
      indirect_HR_lower = exp(bs$lower),
      indirect_HR_upper = exp(bs$upper),
      indirect_p_boot = bs$p_lt0,
      indirect_p_boot_lt0 = bs$p_lt0,
      indirect_p_boot_gt0 = bs$p_gt0,
      indirect_p_boot_two_sided = bs$p_two_sided,
      n_boot_success = as.integer(bs$n)
    )]
  }
  result[, indirect_p_boot_fdr := p.adjust(indirect_p_boot, method = "BH")]
  result[]
}

covariates <- c(
  "age", "gender", "ethnicity", "deprivation", "income", "university_education",
  "smoking", "alcohol", "sedentary", "obesity", "prs", paste("comp", 1:20, sep = "_"),
  "protein_cov_zscore", "season_final", "Batch", "protein_time"
)

analysis <- build_analysis_data()
model_data <- analysis$data
proteins <- analysis$proteins
mediator_cols <- analysis$mediator_cols

model_cols <- c("eid", "time", "status", "pain_score", covariates, mediator_cols)
iorw_data <- as.data.frame(model_data[, model_cols, with = FALSE])
iorw_data <- iorw_data[complete.cases(iorw_data), ]
iorw_data <- droplevels(iorw_data)
fwrite(as.data.table(iorw_data), out_data_path)

cat("running per-protein IORW mediation for", length(mediator_cols), "proteins...\n")
result <- run_observed_per_protein(iorw_data, covariates, mediator_cols, proteins, weight_trim_quantile)

overlap_info <- as.data.table(analysis$overlap)[, .(
  protein,
  beta_pain_score,
  p_value_pain_score,
  fdr_pain_score,
  beta_PD_cox,
  HR,
  p_value_PD_cox,
  fdr_PD_cox
)]
result <- merge(result, overlap_info, by = "protein", all.x = TRUE)

cat("running bootstrap with", n_boot, "resamples per protein...\n")
boot_dt <- bootstrap_per_protein(iorw_data, covariates, mediator_cols, proteins, n_boot, weight_trim_quantile)
if (nrow(boot_dt) > 0) fwrite(boot_dt, out_boot_path)

result <- add_bootstrap_results(result, boot_dt)
setorder(result, indirect_p_boot, protein)
fwrite(result, out_result_path)

write_json_simple(list(
  method = "Per-protein IORW mediation for continuous pain_score exposure and each of the 25 direction-aligned overlap proteins as a separate mediator. For each protein, stabilized density-ratio weights are f(pain_score|covariates)/f(pain_score|covariates,mediator). The weighted Cox model estimates the natural direct effect; the indirect effect is total minus direct on the log-HR scale. indirect_p_boot is exactly the proportion of bootstrap indirect_logHR values below 0.",
  n = nrow(iorw_data),
  events = sum(iorw_data$status),
  n_proteins = length(proteins),
  proteins = proteins,
  weight_trim_quantile = weight_trim_quantile,
  n_boot_requested_per_protein = n_boot,
  n_boot_total_requested = n_boot * length(proteins),
  n_boot_total_success = ifelse(nrow(boot_dt) > 0, sum(is.finite(boot_dt$indirect_logHR)), 0),
  data_output = out_data_path,
  result_output = out_result_path,
  bootstrap_output = ifelse(nrow(boot_dt) > 0, out_boot_path, NA)
), out_summary_path)

cat("outputs:\n")
cat(out_data_path, "\n")
cat(out_result_path, "\n")
if (nrow(boot_dt) > 0) cat(out_boot_path, "\n")
cat(out_summary_path, "\n")
print(result)
