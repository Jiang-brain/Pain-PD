# Run from code/. Existing analysis inputs and outputs are read-only.
CODE_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
REFERENCE_DIR <- Sys.getenv("PD_PAIN_PREDICTION_SOURCE_DIR", file.path(CODE_ROOT, "data/prediction_reference"))
SOURCE_DIR <- file.path(CODE_ROOT, "04_proteomics_mediation/03_shared_protein_mediation", "source_reference")
.libPaths(c(file.path(REFERENCE_DIR, "R_libs"), file.path(CODE_ROOT, "R_libs"), .libPaths()))
library(data.table)
library(survival)
setDTthreads(1L)
SENS_DIR <- file.path(Sys.getenv("PD_PAIN_RESULTS_DIR", file.path(CODE_ROOT, "results")), "sensitivity_nonproteomic_weights")
OUT <- file.path(SENS_DIR, "mediation24")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load_function_definitions <- function(path, environment) {
  for (expr in parse(path)) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=") &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, envir = environment)
    }
  }
}
legacy <- new.env(parent = globalenv())
load_function_definitions(file.path(SOURCE_DIR, "analysis_overlap25_iorw_joint_mediation.R"), legacy)
load_function_definitions(file.path(SOURCE_DIR, "analysis_overlap25_iorw_per_protein_mediation.R"), legacy)

make_design <- function(dat, covariates, mediators) {
  dat <- droplevels(as.data.frame(dat))
  covariates <- legacy$usable_covariates(dat, covariates)
  covariate_matrix <- model.matrix(reformulate(covariates), dat)
  list(covariates = covariates, covariate_matrix = covariate_matrix,
       mediator_matrix = as.matrix(dat[, mediators, drop = FALSE]),
       pain = dat$pain_score, survival = Surv(dat$time, dat$status),
       n = nrow(dat), mediators = mediators)
}

# The same lm.fit and survival::coxph.fit engines as the original formula fits.
# Bootstrap CIs need coefficients only, not repeated robust variance/concordance
# and naive-adjusted fits. Agreement is tested before full execution.
fast_iorw <- function(design, indices = seq_len(design$n), mediator_indices, joint = FALSE,
                      aggregate_censored = TRUE) {
  x <- design$covariate_matrix[indices, , drop = FALSE]
  y <- design$pain[indices]
  m <- design$mediator_matrix[indices, mediator_indices, drop = FALSE]
  numerator <- lm.fit(x, y)
  denominator <- lm.fit(cbind(x, m), y)
  nsd <- sqrt(sum(numerator$residuals^2) / numerator$df.residual)
  dsd <- sqrt(sum(denominator$residuals^2) / denominator$df.residual)
  nd <- dnorm(y, mean = numerator$fitted.values, sd = nsd)
  dd <- dnorm(y, mean = denominator$fitted.values, sd = dsd)
  if (!joint) {
    nd <- pmax(nd, .Machine$double.xmin)
    dd <- pmax(dd, .Machine$double.xmin)
  }
  weights <- nd / dd
  weights[!is.finite(weights)] <- NA_real_
  cutoff <- quantile(weights, 0.99, na.rm = TRUE)
  weights <- pmin(weights, cutoff)
  weights <- weights / mean(weights, na.rm = TRUE)
  if (any(!is.finite(weights)) || any(weights <= 0)) stop("Invalid bootstrap IORW weights")
  cox_x <- cbind(pain_score = y, x[, -1, drop = FALSE])
  cox_y <- design$survival[indices, , drop = FALSE]
  multiplicity <- rep(1, length(indices))
  if (aggregate_censored) {
    # Only duplicated censoring rows are collapsed into equivalent case weights.
    # Every duplicated event remains a separate event row, preserving Efron ties.
    keep <- !duplicated(indices) | cox_y[, "status"] == 1
    frequency <- tabulate(indices, nbins = design$n)
    multiplicity <- frequency[indices[keep]]
    multiplicity[cox_y[keep, "status"] == 1] <- 1
    cox_x <- cox_x[keep, , drop = FALSE]
    cox_y <- cox_y[keep, , drop = FALSE]
    weights <- weights[keep]
  }
  fit <- function(w) survival::coxph.fit(
    x = cox_x, y = cox_y, strata = NULL, offset = rep(0, nrow(cox_x)),
    init = NULL, control = coxph.control(), weights = w, method = "efron",
    rownames = NULL, resid = FALSE
  )
  total <- fit(multiplicity)$coefficients[1]
  direct <- fit(weights * multiplicity)$coefficients[1]
  if (!is.finite(total) || !is.finite(direct)) stop("Non-finite Cox coefficient")
  c(total_logHR = unname(total), direct_logHR = unname(direct),
    indirect_logHR = unname(total - direct), proportion_mediated = unname((total - direct) / total))
}

capture_fit <- function(expression) {
  warnings <- character()
  value <- tryCatch(withCallingHandlers(expression, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(value, "error")) {
    return(list(value = NULL, error = conditionMessage(value), warnings = warnings))
  }
  list(value = value, error = "", warnings = warnings)
}
