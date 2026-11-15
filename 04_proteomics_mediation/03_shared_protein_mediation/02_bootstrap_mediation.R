source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
args <- commandArgs(trailingOnly = TRUE)
inputs <- readRDS(file.path(OUT, "prepared_inputs.rds"))
implementation_md5 <- unname(tools::md5sum("04_proteomics_mediation/03_shared_protein_mediation/00_common.R"))
if (!all(unname(tools::md5sum(inputs$input_hashes$path)) == inputs$input_hashes$md5)) stop("Input files changed")
jobs_path <- file.path(OUT, "bootstrap_jobs.rds")

if (length(args) && args[1] == "--prepare-jobs") {
  if (file.exists(jobs_path)) stop("Job schedule exists; do not overwrite seeds/checkpoints")
  jobs <- list()
  # Advance the original serial sample.int stream without altering draws.
  # Removing COMP (last in the original order) leaves the other 24 streams intact.
  for (group in c("single", "joint")) {
    set.seed(inputs$seed)
    models <- if (group == "single") inputs$proteins else "joint24"
    for (model in models) {
      for (first in seq(1L, inputs$n_boot, by = 20L)) {
        last <- min(first + 19L, inputs$n_boot)
        jobs[[length(jobs) + 1L]] <- list(model = model, first = first, last = last,
                                         rng_state = .Random.seed)
        for (b in seq.int(first, last)) sample.int(inputs$design$n, inputs$design$n, replace = TRUE)
      }
      cat("Prepared original bootstrap stream:", model, "\n"); flush.console()
    }
  }
  saveRDS(jobs, jobs_path)
  cat("Prepared", length(jobs), "atomic bootstrap jobs.\n")
  quit(save = "no")
}

if (length(args) != 2L) stop("Use --prepare-jobs or WORKER_INDEX WORKER_COUNT")
worker <- as.integer(args[1]); workers <- as.integer(args[2])
if (is.na(worker) || is.na(workers) || worker < 1L || worker > workers) stop("Invalid worker arguments")
jobs <- readRDS(jobs_path)
checkpoint_dir <- file.path(OUT, "bootstrap_checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE)
for (job_index in seq_along(jobs)) {
  if ((job_index - 1L) %% workers != worker - 1L) next
  job <- jobs[[job_index]]
  output <- file.path(checkpoint_dir, sprintf("%s_%04d_%04d.rds", job$model, job$first, job$last))
  if (file.exists(output)) next
  assign(".Random.seed", job$rng_state, envir = .GlobalEnv)
  mediator_indices <- if (job$model == "joint24") seq_along(inputs$mediators) else match(job$model, inputs$proteins)
  rows <- list()
  started <- proc.time()[3]
  for (b in seq.int(job$first, job$last)) {
    idx <- sample.int(inputs$design$n, inputs$design$n, replace = TRUE)
    fitted <- capture_fit(fast_iorw(inputs$design, idx, mediator_indices, joint = job$model == "joint24"))
    row <- data.table(model = job$model, bootstrap = b, valid = !nzchar(fitted$error),
                       error = fitted$error, warnings = paste(fitted$warnings, collapse = " | "))
    if (!is.null(fitted$value)) for (field in names(fitted$value)) row[[field]] <- fitted$value[field]
    rows[[length(rows) + 1L]] <- row
  }
  checkpoint <- list(job_index = job_index, model = job$model, first = job$first,
                     last = job$last, rows = rbindlist(rows, fill = TRUE),
                     elapsed_seconds = proc.time()[3] - started,
                     implementation_md5 = implementation_md5,
                     input_hashes = inputs$input_hashes)
  temporary <- paste0(output, ".tmp")
  saveRDS(checkpoint, temporary)
  if (!file.rename(temporary, output)) stop("Could not atomically save checkpoint")
  cat("worker", worker, "completed", job$model, job$first, "-", job$last,
      "in", round(checkpoint$elapsed_seconds, 1), "seconds\n"); flush.console()
}
