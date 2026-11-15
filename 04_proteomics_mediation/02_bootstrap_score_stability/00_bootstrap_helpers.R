source("04_proteomics_mediation/01_independent_pain_score/00_common.R")
BOOT_CODE_DIR <- file.path(CODE_ROOT_SENS, "04_proteomics_mediation/02_bootstrap_score_stability")
BOOT_OUT_DIR <- file.path(SENS_OUT_DIR, "bootstrap1000")
dir.create(BOOT_OUT_DIR, showWarnings = FALSE, recursive = TRUE)
boot_output <- function(name) file.path(BOOT_OUT_DIR, name)

fit_bootstrap_matrix <- function(x, y, control) {
  survival::coxph.fit(x = x, y = y, strata = NULL, offset = NULL,
                     init = NULL, control = control, weights = NULL,
                     method = "efron", rownames = NULL, resid = FALSE,
                     nocenter = c(-1, 0, 1))
}
score_covariance_correlation <- function(indicators, normalized_weights) {
  covariance <- stats::cov(indicators)
  gram <- crossprod(normalized_weights, covariance %*% normalized_weights)
  sd <- sqrt(diag(gram))
  if (any(!is.finite(sd)) || any(sd <= 0)) stop("A bootstrap score has zero variance.")
  answer <- gram / outer(sd, sd)
  answer[answer > 1 & answer < 1 + 1e-12] <- 1
  answer[answer < -1 & answer > -1 - 1e-12] <- -1
  diag(answer) <- 1
  answer
}
weighted_midranks <- function(values, frequency) {
  ordering <- order(values)
  sorted_values <- values[ordering]
  groups <- cumsum(c(TRUE, diff(sorted_values) != 0))
  group_n <- as.numeric(rowsum(frequency[ordering], groups, reorder = FALSE))
  group_end <- cumsum(group_n)
  midrank <- group_end - (group_n - 1) / 2
  answer <- numeric(length(values))
  answer[ordering] <- midrank[groups]
  answer
}
weighted_column_correlation <- function(x, frequency) {
  means <- colSums(x * frequency) / sum(frequency)
  centered <- sweep(x, 2L, means, "-") * sqrt(frequency)
  gram <- crossprod(centered)
  sd <- sqrt(diag(gram))
  if (any(sd <= 0) || any(!is.finite(sd))) stop("Zero variance in weighted correlation.")
  answer <- gram / outer(sd, sd)
  answer[answer > 1 & answer < 1 + 1e-12] <- 1
  diag(answer) <- 1
  answer
}
summarize_correlations <- function(values, comparison, population, method) {
  if (any(!is.finite(values))) stop("Non-finite correlation.")
  qs <- quantile(values, c(0, .025, .05, .25, .5, .75, .95, .975, 1), names = FALSE)
  data.table(comparison = comparison, population = population, method = method,
             n_correlations = length(values), mean = mean(values), sd = sd(values),
             minimum = qs[1], p025 = qs[2], p05 = qs[3], p25 = qs[4], median = qs[5],
             p75 = qs[6], p95 = qs[7], p975 = qs[8], maximum = qs[9],
             fraction_ge_09 = mean(values >= .9), fraction_ge_095 = mean(values >= .95))
}
