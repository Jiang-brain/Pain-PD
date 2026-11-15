# Recompute summaries from saved bootstrap effects; do not fit models or write results.
source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
boot <- fread(file.path(OUT, "mediation24_bootstrap.csv"))
observed <- fread(file.path(OUT, "per_protein_observed.csv"))
single_boot <- copy(boot[model != "joint24"])
single_boot[, `:=`(protein = model, mediator = paste0("aligned_", model, "_z"))]
recomputed <- legacy$add_bootstrap_results(observed, single_boot)
saved <- fread(file.path(OUT, "per_protein_mediation_results.csv"))
setorder(recomputed, protein)
setorder(saved, protein)
stopifnot(identical(recomputed$protein, saved$protein))
for (column in c("indirect_logHR", "indirect_lower", "indirect_upper", "proportion_mediated")) {
  stopifnot(max(abs(recomputed[[column]] - saved[[column]])) < 1e-12)
}
flag <- recomputed$indirect_lower > 0 | recomputed$indirect_upper < 0
stopifnot(sum(flag) == 16L)
stopifnot(sum(saved$indirect_p_boot_fdr < 0.05) == 18L)
stopifnot(sum(saved$indirect_p_boot_two_sided_fdr < 0.05) == 14L)
for (protein_name in saved$protein) {
  b <- single_boot[protein == protein_name & valid == TRUE]
  r <- saved[protein == protein_name]
  stopifnot(max(abs(as.numeric(quantile(b$proportion_mediated, c(0.025, 0.975))) -
                    c(r$proportion_lower, r$proportion_upper))) < 1e-12)
}
joint <- fread(file.path(OUT, "joint24_mediation_results.csv"))
stopifnot(joint$indirect_lower > 0, abs(joint$proportion_mediated * 100 - 17.8175100195335) < 1e-10)
cat("PASS: saved effects/CIs unchanged; 16 CI-positive, legacy one-sided FDR 18, two-sided FDR 14. No models refitted.\n")
