source("04_proteomics_mediation/02_bootstrap_score_stability/00_bootstrap_helpers.R")
required <- c("fixed_weights.rds", "derivation_participant_ids.csv", "composite_pain_scores.csv",
              "protein_roster_sources.csv", "validation_checks.csv")
if (!all(file.exists(vapply(required, sens_output, "")))) stop("Run sensitivity branch 07 first.")
if (!all(fread(sens_output("validation_checks.csv"))$passed)) stop("Previous sensitivity checks did not pass.")
frozen <- readRDS(sens_output("fixed_weights.rds"))
training_ids <- fread(sens_output("derivation_participant_ids.csv"))$eid
target_scores <- fread(sens_output("composite_pain_scores.csv"))
assert_unique_eid(target_scores, "Bootstrap scoring group")

expressions <- parse(reference_path("PD_Pain_WQ.R"))
fit_index <- which(vapply(expressions, assigned_name, "") == "hr")
if (length(fit_index) != 1L) stop("Unexpected original script structure.")
original <- new.env(parent = globalenv())
original$data_path <- sens_input("PD_data.csv")
cat("Reconstructing the unchanged PA cohort...\n")
for (i in seq_len(fit_index - 1L)) {
  expr <- expressions[[i]]
  if (assigned_name(expr) == "data_path" || is_libpaths(expr)) next
  eval(expr, envir = original)
}
assert_unique_eid(original$datas, "Original score cohort")
train_rows <- match(training_ids, original$datas$eid)
target_rows <- match(target_scores$eid, original$datas$eid)
if (anyNA(train_rows) || anyNA(target_rows) || anyDuplicated(training_ids)) stop("Cohort IDs do not match.")
roster_paths <- fread(sens_output("protein_roster_sources.csv"))$path
protein_ids <- unique(unlist(lapply(roster_paths, function(path) fread(path, select = "eid")$eid)))
if (any(training_ids %in% protein_ids) || any(training_ids %in% target_scores$eid)) stop("Derivation leakage.")
if (!all(target_scores$eid %in% protein_ids)) stop("Scored group includes participants without proteins.")
model_formula <- as.formula(frozen$formula)
derivation <- original$new_data[train_rows, , drop = FALSE]
cat("Fitting reference model in", nrow(derivation), "participants...\n")
timing <- system.time(reference_fit <- coxph(model_formula, data = derivation, x = TRUE, y = TRUE))
if (reference_fit$n != length(training_ids)) stop("Unexpected exclusions from reference fit.")
pain_names <- frozen$weights$coefficient
reference_beta <- coef(reference_fit)
if (max(abs(reference_beta[pain_names] - frozen$weights$raw_beta)) >= 1e-10) stop("Reference beta does not reproduce.")
target_indicators <- vapply(frozen$weights$term, function(term) {
  as.numeric(as.character(original$new_data[[term]][target_rows]))
}, numeric(length(target_rows)))
colnames(target_indicators) <- pain_names
reference_weights <- setNames(frozen$weights$effective_beta / frozen$denominator, pain_names)
recomputed <- as.numeric(target_indicators %*% reference_weights)
if (max(abs(recomputed - target_scores$pain_score)) >= 1e-12) stop("Reference scores do not reproduce.")

control <- coxph.control()
direct_time <- system.time(direct_fit <- fit_bootstrap_matrix(reference_fit$x, reference_fit$y, control))
direct_difference <- max(abs(direct_fit$coefficients - reference_beta))
if (!is.finite(direct_difference) || direct_difference >= 1e-10) stop("Matrix fit differs from coxph.")

# Validate ordinary row resampling, including duplicated event rows and Efron
# ties, against coxph. Case multiplicity weights are deliberately not used.
seed <- 20261026L
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(seed)
replicate_seeds <- sample.int(.Machine$integer.max, 1000L, replace = FALSE)
set.seed(replicate_seeds[1])
idx <- sample.int(nrow(derivation), nrow(derivation), replace = TRUE)
bootstrap_time <- system.time(test_direct <- fit_bootstrap_matrix(reference_fit$x[idx, , drop = FALSE],
                                                                 reference_fit$y[idx, , drop = FALSE], control))
test_standard <- coxph(model_formula, data = derivation[idx, , drop = FALSE])
bootstrap_difference <- max(abs(test_direct$coefficients - coef(test_standard)))
if (!is.finite(bootstrap_difference) || bootstrap_difference >= 1e-8) stop("Bootstrap matrix fit differs from coxph.")

model_input <- list(x = reference_fit$x, y = reference_fit$y, control = control,
                    pain_names = pain_names, all_names = names(reference_beta),
                    replicate_seeds = replicate_seeds, master_seed = seed, n_bootstrap = 1000L)
target_input <- list(eid = target_scores$eid, indicators = target_indicators,
                     reference_weights = reference_weights, reference_scores = target_scores$pain_score,
                     any_pain = rowSums(target_indicators) > 0,
                     reference_all_beta = reference_beta)
saveRDS(model_input, boot_output("bootstrap_model_input.rds"), compress = FALSE)
saveRDS(target_input, boot_output("bootstrap_scoring_input.rds"), compress = FALSE)
fwrite(data.table(replicate = seq_len(1000L), seed = replicate_seeds), boot_output("bootstrap_seeds.csv"))
fwrite(data.table(n_derivation = length(training_ids), derivation_events = sum(reference_fit$y[, 2]),
                  n_scored = nrow(target_indicators), n_scored_with_pain = sum(target_input$any_pain),
                  derivation_protein_overlap = 0L, n_bootstrap = 1000L, master_seed = seed,
                  baseline_coxph_elapsed_seconds = unname(timing[3]),
                  baseline_matrix_elapsed_seconds = unname(direct_time[3]),
                  one_bootstrap_elapsed_seconds = unname(bootstrap_time[3]),
                  baseline_coxph_matrix_max_beta_difference = direct_difference,
                  bootstrap_coxph_matrix_max_beta_difference = bootstrap_difference),
       boot_output("preparation_checks.csv"))
paths <- c(reference_path("PD_Pain_WQ.R"), sens_output("fixed_weights.rds"),
           sens_output("derivation_participant_ids.csv"), sens_output("composite_pain_scores.csv"))
fwrite(data.table(path = paths, md5 = unname(tools::md5sum(paths))), boot_output("input_source_hashes.csv"))
writeLines(capture.output(sessionInfo()), boot_output("sessionInfo_preparation.txt"))
cat("Prepared 1000 bootstrap seeds; one matrix bootstrap fit took", bootstrap_time[3], "seconds.\n")
