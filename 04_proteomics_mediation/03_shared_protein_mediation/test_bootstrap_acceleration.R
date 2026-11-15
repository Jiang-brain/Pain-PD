source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
inputs <- readRDS(file.path(OUT, "prepared_inputs.rds"))
set.seed(inputs$seed)
idx <- sample.int(inputs$design$n, inputs$design$n, replace = TRUE)
checks <- list()
for (mi in c(1L, 5L, 0L)) {
  joint <- mi == 0L
  m <- if (joint) seq_along(inputs$mediators) else mi
  started <- proc.time()[3]
  reduced <- fast_iorw(inputs$design, idx, m, joint, aggregate_censored = TRUE)
  elapsed <- proc.time()[3] - started
  repeated <- fast_iorw(inputs$design, idx, m, joint, aggregate_censored = FALSE)
  difference <- max(abs(reduced - repeated))
  checks[[length(checks) + 1L]] <- data.table(model = if (joint) "joint24" else inputs$proteins[mi],
    max_coefficient_difference = difference, reduced_elapsed_seconds = elapsed, passed = difference < 1e-8)
}
checks <- rbindlist(checks)
fwrite(checks, file.path(OUT, "censored_duplicate_equivalence_checks.csv"))
print(checks)
if (!all(checks$passed)) stop("Censor-only reduction does not preserve Efron results")
