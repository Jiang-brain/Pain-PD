source("04_proteomics_mediation/01_independent_pain_score/00_common.R")
new_result_path <- sens_output("pain_score_protein_association_results.csv")
if (!file.exists(new_result_path)) stop("Run sensitivity step 02 first.")
pain <- fread(new_result_path)
original_pain <- fread(sens_input("pain_score_protein_association_results.csv"))
pd <- fread(sens_input("protein_PD_cox_results.csv"))
if (anyDuplicated(pain$protein) || anyDuplicated(pd$protein)) stop("Duplicate protein results.")
combine <- function(pain_result) {
  joined <- merge(pain_result, pd, by = "protein", suffixes = c("_pain", "_PD"))
  joined[, directionally_concordant := is.finite(beta_pain) & is.finite(beta_PD) & beta_pain * beta_PD > 0]
  joined
}
joined <- combine(pain)
baseline <- combine(original_pain)
baseline_shared <- baseline[fdr_pain < 0.05 & fdr_PD < 0.05 & directionally_concordant]
shared <- joined[fdr_pain < 0.05 & fdr_PD < 0.05 & directionally_concordant]
fwrite(joined, sens_output("protein_pain_PD_merged_all.csv"))
for (criterion in c("fdr", "bonferroni")) {
  overlap <- joined[get(paste0(criterion, "_pain")) < 0.05 & get(paste0(criterion, "_PD")) < 0.05]
  fwrite(overlap, sens_output(paste0("protein_intersection_", criterion, "05_all_directions.csv")))
  fwrite(overlap[directionally_concordant == TRUE],
         sens_output(paste0("shared_proteins_", criterion, "05_concordant.csv")))
}
membership <- merge(baseline[, .(protein, original_pain_beta = beta_pain, original_pain_fdr = fdr_pain,
                                 PD_beta = beta_PD, PD_fdr = fdr_PD)],
                    joined[, .(protein, sensitivity_pain_beta = beta_pain, sensitivity_pain_fdr = fdr_pain)],
                    by = "protein", all = TRUE)
membership[, original_shared := protein %in% baseline_shared$protein]
membership[, sensitivity_shared := protein %in% shared$protein]
membership[, shared_status := fifelse(original_shared & sensitivity_shared, "retained",
                                     fifelse(original_shared, "lost", fifelse(sensitivity_shared, "gained", "neither")))]
fwrite(membership, sens_output("shared_protein_membership_comparison.csv"))
fwrite(membership[shared_status != "neither"], sens_output("shared_protein_retained_lost_gained.csv"))
counts <- data.table(
  analysis = c("original", "nonproteomic_fixed_weights"),
  pain_fdr = c(sum(original_pain$fdr < 0.05, na.rm = TRUE), sum(pain$fdr < 0.05, na.rm = TRUE)),
  PD_fdr = sum(pd$fdr < 0.05, na.rm = TRUE),
  overlap_fdr_all_directions = c(nrow(baseline[fdr_pain < 0.05 & fdr_PD < 0.05]),
                                  nrow(joined[fdr_pain < 0.05 & fdr_PD < 0.05])),
  shared_fdr_concordant = c(nrow(baseline_shared), nrow(shared)))
fwrite(counts, sens_output("shared_protein_summary.csv"))
if (nrow(baseline_shared) != 25L) warning("Baseline shared set differs from the manuscript's 25; inspect source results.")
cat("Shared proteins retained:", sum(membership$shared_status == "retained"),
    "; lost:", sum(membership$shared_status == "lost"),
    "; gained:", sum(membership$shared_status == "gained"), "\n")
print(counts)
write_provenance(sens_input("protein_PD_cox_results.csv"), "unchanged_PD_result_hash.csv")
