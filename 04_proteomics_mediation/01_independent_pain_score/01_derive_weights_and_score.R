source("04_proteomics_mediation/01_independent_pain_score/00_common.R")

# Reuse the untouched PA cohort construction, not reorganized helper recodings.
expressions <- parse(reference_path("PD_Pain_WQ.R"))
fit_index <- which(vapply(expressions, assigned_name, "") == "hr")
if (length(fit_index) != 1L) stop("Expected exactly one original mutual Cox model.")
original <- new.env(parent = globalenv())
original$data_path <- sens_input("PD_data.csv")
cat("Preparing the original PA score cohort...\n")
for (i in seq_len(fit_index - 1L)) {
  expr <- expressions[[i]]
  if (assigned_name(expr) == "data_path" || is_libpaths(expr)) next
  eval(expr, envir = original)
}
assert_unique_eid(original$datas, "Original score cohort")

# The model uses 14 interleaved acute/chronic terms, whereas the data columns
# have all acute terms first. Match coefficient names instead of column position.
extract_weights <- function(fit) {
  beta <- coef(fit)[seq_len(14L)]
  if (anyNA(beta) || any(!is.finite(beta))) stop("A pain-term coefficient is not estimable.")
  if (any(!grepl("_(acute|chronic)1$", names(beta)))) stop("Unexpected pain coefficient names.")
  effective <- pmax(beta, 0)
  if (!is.finite(sum(effective)) || sum(effective) <= 0) stop("Invalid weight denominator.")
  data.frame(coefficient = names(beta), term = sub("1$", "", names(beta)),
             raw_beta = as.numeric(beta), effective_beta = as.numeric(effective))
}
apply_fixed_weights <- function(pain_data, weights) {
  indicators <- vapply(weights$term, function(term) {
    if (!(term %in% names(pain_data))) stop("Missing pain indicator: ", term)
    as.numeric(as.character(pain_data[[term]]))
  }, numeric(nrow(pain_data)))
  if (anyNA(indicators) || any(!indicators %in% c(0, 1))) stop("Invalid pain indicators.")
  as.numeric(indicators %*% weights$effective_beta) / sum(weights$effective_beta)
}

cat("Reproducing the original all-cohort score as a mandatory check...\n")
eval(expressions[[fit_index]], envir = original)
baseline_weights <- extract_weights(original$hr)
baseline_scores <- data.table(eid = original$datas$eid,
                             reproduced = apply_fixed_weights(original$new_data, baseline_weights))
saved_scores <- fread(sens_input("composite_pain_scores.csv"), select = c("eid", "pain_score"))
assert_unique_eid(saved_scores, "Saved original scores")
if (!setequal(baseline_scores$eid, saved_scores$eid)) stop("Original score cohort IDs do not reproduce.")
comparison <- merge(baseline_scores, saved_scores, by = "eid")
max_difference <- max(abs(comparison$reproduced - comparison$pain_score))
fwrite(data.table(n = nrow(comparison), max_absolute_difference = max_difference,
                  passed = is.finite(max_difference) && max_difference < 1e-12),
       sens_output("baseline_score_reproduction_check.csv"))
if (!is.finite(max_difference) || max_difference >= 1e-12) {
  stop("Original score reproduction failed; no sensitivity weights will be fitted.")
}
cat("Original score reproduced; maximum absolute difference:", max_difference, "\n")

# Exclude every measured participant, including those omitted by protein QC.
# Always include filtered-matrix IDs as a second guard against partial rosters.
membership_files <- sens_input("protein_new_missing60_filtered.csv")
raw_matrix <- file.path(SENS_DATA_DIR, "protein_new.csv")
roster_override <- Sys.getenv("PD_PAIN_PROTEIN_ROSTER", "")
if (nzchar(roster_override)) {
  if (!file.exists(roster_override)) stop("Protein roster does not exist: ", roster_override)
  membership_files <- unique(c(membership_files, roster_override))
} else if (file.exists(raw_matrix)) {
  membership_files <- c(membership_files, raw_matrix)
} else {
  stop("Provide protein_new.csv or PD_PAIN_PROTEIN_ROSTER (a CSV containing eid). ",
       "QC-filtered IDs alone cannot establish absence of measured proteins.")
}
protein_ids <- unique(unlist(lapply(membership_files, function(path) fread(path, select = "eid")$eid)))
if (anyNA(protein_ids)) stop("Missing eid in protein membership roster.")
is_protein <- original$datas$eid %in% protein_ids
if (!any(is_protein) || !any(!is_protein)) stop("Both derivation and protein groups must be nonempty.")
training_ids <- original$datas$eid[!is_protein]
if (length(intersect(training_ids, protein_ids)) != 0L) stop("Derivation overlaps protein roster.")

model_formula <- formula(original$hr)
derivation <- original$new_data[!is_protein, , drop = FALSE]
cat("Fitting fixed weights in", nrow(derivation), "non-proteomic participants...\n")
external_fit <- coxph(model_formula, data = derivation)
if (external_fit$n != nrow(derivation)) stop("Unexpected missingness in derivation model.")
weights <- extract_weights(external_fit)
new_scores <- apply_fixed_weights(original$new_data, weights)
positive_training_scores <- new_scores[!is_protein & new_scores > 0]
tertiles <- quantile(positive_training_scores, c(1 / 3, 2 / 3), na.rm = TRUE)
categories <- rep(0L, length(new_scores))
categories[new_scores > 0 & new_scores < tertiles[1]] <- 1L
categories[new_scores >= tertiles[1] & new_scores <= tertiles[2] & new_scores > 0] <- 2L
categories[new_scores > tertiles[2]] <- 3L
output_scores <- data.table(eid = original$datas$eid[is_protein], pain_score = new_scores[is_protein],
                           weighted_pain_status = categories[is_protein])
fwrite(output_scores, sens_output("composite_pain_scores.csv"))
weight_comparison <- merge(baseline_weights, weights, by = c("coefficient", "term"),
                           suffixes = c("_original", "_nonproteomic"), sort = FALSE)
fwrite(weight_comparison, sens_output("pain_score_fixed_beta_comparison.csv"))
fwrite(data.table(eid = training_ids), sens_output("derivation_participant_ids.csv"))
fwrite(data.table(group = c("eligible_original_cohort", "nonproteomic_derivation", "proteomic_scoring"),
                  n = c(nrow(original$datas), sum(!is_protein), sum(is_protein)),
                  events = c(sum(original$datas$status), sum(original$datas$status[!is_protein]),
                             sum(original$datas$status[is_protein]))), sens_output("cohort_counts.csv"))
fwrite(data.table(n_protein_roster_unique = length(protein_ids), derivation_roster_overlap = 0L,
                  weight_denominator = sum(weights$effective_beta),
                  lower_tertile_derivation = tertiles[1], upper_tertile_derivation = tertiles[2]),
       sens_output("fixed_weight_checks.csv"))
fwrite(comparison[eid %in% output_scores$eid,
                  .(eid, original_pain_score = pain_score)][output_scores, on = "eid"],
       sens_output("proteomic_score_comparison.csv"))
saveRDS(list(weights = weights, denominator = sum(weights$effective_beta),
             formula = paste(deparse(model_formula), collapse = " "),
             derivation_n = nrow(derivation), derivation_events = sum(derivation$status),
             clipped_negative_beta = TRUE, tertiles = tertiles), sens_output("fixed_weights.rds"))
writeLines(paste(deparse(model_formula), collapse = " "), sens_output("weight_model_formula.txt"))
write_provenance(c(reference_path("PD_Pain_WQ.R"), sens_input("composite_pain_scores.csv")),
                 "score_source_hashes.csv")
fwrite(data.table(path = normalizePath(membership_files, winslash = "/", mustWork = TRUE)),
       sens_output("protein_roster_sources.csv"))
writeLines(capture.output(sessionInfo()), sens_output("sessionInfo_score.txt"))
cat("Protein scores:", nrow(output_scores), "; derivation-protein overlap: 0\n")
