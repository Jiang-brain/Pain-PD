source("04_proteomics_mediation/02_bootstrap_score_stability/00_bootstrap_helpers.R")
input_path <- boot_output("bootstrap_model_input.rds")
if (!file.exists(input_path)) stop("Run bootstrap step 01 first.")
model_input <- readRDS(input_path)
if (model_input$n_bootstrap != 1000L) stop("This analysis requires exactly 1000 replicates.")
args <- commandArgs(trailingOnly = TRUE)
collect_only <- identical(args, "--collect-only")
if (!collect_only && !length(args)) args <- c("1", "1")
if (!collect_only && length(args) != 2L) stop("Use: Rscript step02.R SHARD_INDEX SHARD_COUNT, or --collect-only.")
shard_index <- if (collect_only) 1L else as.integer(args[1])
shard_count <- if (collect_only) 1L else as.integer(args[2])
if (anyNA(c(shard_index, shard_count)) || shard_count < 1L || shard_index < 1L || shard_index > shard_count) {
  stop("Invalid shard index/count.")
}
checkpoint_dir <- boot_output("checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE)
checkpoint_path <- function(b) file.path(checkpoint_dir, sprintf("replicate_%04d.rds", b))
model_hash <- unname(tools::md5sum(input_path))
existing <- file.exists(vapply(seq_len(1000L), checkpoint_path, ""))
for (b in which(existing)) {
  saved <- readRDS(checkpoint_path(b))
  if (!identical(saved$model_hash, model_hash) || saved$seed != model_input$replicate_seeds[b]) {
    stop("Checkpoint input mismatch: ", b)
  }
}
assigned <- which((seq_len(1000L) - 1L) %% shard_count == shard_index - 1L)
pending <- intersect(which(!existing), assigned)
cat("Existing checkpoints:", sum(existing), "; shard:", shard_index, "/", shard_count,
    "; pending in shard:", length(pending), "\n")

bootstrap_worker <- function(b) {
  started <- proc.time()[3]
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(model_input$replicate_seeds[b])
    idx <- sample.int(nrow(model_input$x), nrow(model_input$x), replace = TRUE)
    fit <- fit_bootstrap_matrix(model_input$x[idx, , drop = FALSE], model_input$y[idx, , drop = FALSE],
                                model_input$control)
    beta <- fit$coefficients
    pain_beta <- beta[model_input$pain_names]
    effective <- pmax(pain_beta, 0)
    denominator <- sum(effective)
    valid <- all(is.finite(beta)) && is.finite(denominator) && denominator > 0
    list(raw_beta = beta, effective_beta = effective, denominator = denominator,
         n_sampled = length(idx), n_unique = length(unique(idx)), events = sum(model_input$y[idx, 2]),
         iterations = fit$iter, valid = valid, error = "")
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(valid = FALSE, error = conditionMessage(e)))
  result$warnings <- warnings
  result$valid <- isTRUE(result$valid) && length(warnings) == 0L
  result$replicate <- b
  result$seed <- model_input$replicate_seeds[b]
  result$model_hash <- model_hash
  result$elapsed_seconds <- unname(proc.time()[3] - started)
  destination <- file.path(checkpoint_dir, sprintf("replicate_%04d.rds", b))
  temporary <- paste0(destination, ".tmp")
  saveRDS(result, temporary)
  if (!file.rename(temporary, destination)) stop("Could not commit checkpoint: ", b)
  c(replicate = b, valid = result$valid, elapsed_seconds = result$elapsed_seconds)
}

if (!collect_only) {
  for (i in seq_along(pending)) {
    result <- bootstrap_worker(pending[i])
    if (i %% 10L == 0L || i == length(pending)) {
      cat("Shard", shard_index, ": completed", i, "/", length(pending), "pending fits.\n")
      flush.console()
    }
    if (result[["valid"]] != 1) stop("Fit failed or warned at replicate ", pending[i], "; checkpoint retained.")
  }
  cat("Shard", shard_index, "finished. Run --collect-only after all shards finish.\n")
  quit(save = "no", status = 0L)
}

if (!all(file.exists(vapply(seq_len(1000L), checkpoint_path, "")))) stop("Not all 1000 checkpoints exist yet.")
replicates <- lapply(seq_len(1000L), function(b) readRDS(checkpoint_path(b)))
diagnostics <- rbindlist(lapply(replicates, function(x) data.table(
  replicate = x$replicate, seed = x$seed, valid = x$valid,
  n_sampled = if (is.null(x$n_sampled)) NA_integer_ else x$n_sampled,
  n_unique = if (is.null(x$n_unique)) NA_integer_ else x$n_unique,
  events = if (is.null(x$events)) NA_integer_ else x$events,
  iterations = if (is.null(x$iterations)) NA_integer_ else x$iterations,
  denominator = if (is.null(x$denominator)) NA_real_ else x$denominator,
  elapsed_seconds = x$elapsed_seconds, warnings = paste(x$warnings, collapse = " | "), error = x$error)))
fwrite(diagnostics, boot_output("bootstrap_fit_diagnostics.csv"))
if (!all(diagnostics$valid)) stop("Some bootstrap fits failed or warned; no replicates were silently discarded.")
raw_beta <- vapply(replicates, function(x) x$raw_beta[model_input$all_names], numeric(length(model_input$all_names)))
effective_beta <- vapply(replicates, function(x) x$effective_beta[model_input$pain_names], numeric(14L))
rownames(raw_beta) <- model_input$all_names
rownames(effective_beta) <- model_input$pain_names
colnames(raw_beta) <- colnames(effective_beta) <- sprintf("bootstrap_%04d", seq_len(1000L))
normalized_weights <- sweep(effective_beta, 2L, diagnostics$denominator, "/")
saveRDS(list(raw_beta = raw_beta, effective_beta = effective_beta,
             normalized_weights = normalized_weights, diagnostics = diagnostics),
        boot_output("bootstrap_coefficients.rds"))
fwrite(data.table(coefficient = rownames(raw_beta), as.data.table(raw_beta)), boot_output("all_model_beta_1000_columns.csv"))
pain_raw_beta <- raw_beta[model_input$pain_names, , drop = FALSE]
fwrite(data.table(coefficient = rownames(pain_raw_beta), as.data.table(pain_raw_beta)), boot_output("pain_raw_beta_1000_columns.csv"))
fwrite(data.table(coefficient = rownames(effective_beta), as.data.table(effective_beta)), boot_output("pain_effective_beta_1000_columns.csv"))
fwrite(data.table(coefficient = rownames(normalized_weights), as.data.table(normalized_weights)), boot_output("pain_normalized_weight_1000_columns.csv"))
writeLines(capture.output(sessionInfo()), boot_output("sessionInfo_bootstrap.txt"))
cat("All 1000 bootstrap fits completed without warnings or failed fits.\n")
