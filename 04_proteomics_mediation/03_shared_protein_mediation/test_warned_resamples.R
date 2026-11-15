source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
inputs <- readRDS(file.path(OUT, "prepared_inputs.rds"))
jobs <- readRDS(file.path(OUT, "bootstrap_jobs.rds"))
files <- list.files(file.path(OUT, "bootstrap_checkpoints"), pattern = "\\.rds$", full.names = TRUE)
checkpoints <- lapply(files, readRDS)
rows <- rbindlist(lapply(checkpoints, function(checkpoint) {
  x <- copy(checkpoint$rows)
  x[, job_index := checkpoint$job_index]
  x
}), fill = TRUE)
examples <- rows[nzchar(warnings)]
examples[, warning_class := ifelse(grepl("38,39,40,41,42,43,44", warnings, fixed = TRUE), "batch_block", "other")]
examples <- examples[, .SD[1], by = warning_class]
checks <- list()
for (i in seq_len(nrow(examples))) {
  row <- examples[i]
  job <- jobs[[row$job_index]]
  assign(".Random.seed", job$rng_state, envir = .GlobalEnv)
  for (b in seq.int(job$first, row$bootstrap)) idx <- sample.int(inputs$design$n, inputs$design$n, replace = TRUE)
  mi <- match(row$model, inputs$proteins)
  if (is.na(mi)) stop("This warning test expects a single-protein example")
  original <- capture_fit(legacy$compute_iorw_single(inputs$data[idx, ], inputs$covariates, inputs$mediators[mi], row$model))
  if (!is.null(original$value)) {
    fields <- c("total_logHR", "direct_logHR", "indirect_logHR", "proportion_mediated")
    difference <- max(abs(unlist(original$value[, fields, with = FALSE]) - unlist(row[, fields, with = FALSE])))
    checks[[length(checks) + 1L]] <- data.table(model = row$model, bootstrap = row$bootstrap,
      warning_class = row$warning_class, max_difference = difference,
      original_formula_warnings = paste(original$warnings, collapse = " | "),
      fast_warnings = row$warnings, coefficients_agree = difference < 1e-8)
  } else stop(original$error)
}
checks <- rbindlist(checks)
fwrite(checks, file.path(OUT, "warned_resample_original_formula_checks.csv"))
fwrite(data.table(variable_number = seq_len(ncol(inputs$design$covariate_matrix)),
  variable = c("pain_score", colnames(inputs$design$covariate_matrix)[-1])),
  file.path(OUT, "Cox_design_variable_map.csv"))
print(checks)
if (!all(checks$coefficients_agree)) stop("Warned resamples do not reproduce the original formula fits")
