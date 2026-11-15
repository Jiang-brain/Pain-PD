source("04_proteomics_mediation/03_shared_protein_mediation/00_common.R")
inputs <- readRDS(file.path(OUT, "prepared_inputs.rds"))
source_path <- file.path(SOURCE_DIR, "iorw_joint_mediation_joint25_and_each_protein_standardized_beta.R")
path_functions <- new.env(parent = globalenv())
load_function_definitions(source_path, path_functions)
dat <- inputs$data
dat$pain_score_z <- path_functions$zscore_numeric(dat$pain_score)
paths <- path_functions$compute_each_protein_standardized_path_abc_beta(
  dat, inputs$covariates, inputs$mediators, proteins = inputs$proteins,
  exposure_z = "pain_score_z"
)
# The immutable helper accepts arbitrary mediator vectors. Update its historical
# output labels, not its fitted models, to match the actual 24-protein input.
names(paths) <- gsub("joint25", "joint24", names(paths), fixed = TRUE)
for (column in names(paths)[vapply(paths, is.character, logical(1))]) {
  paths[[column]] <- gsub("all_25", "all_24", paths[[column]], fixed = TRUE)
}
fwrite(paths, file.path(OUT, "per_protein_standardized_path_coefficients.csv"))
joint <- fread(file.path(OUT, "joint24_observed.csv"))
sd_pain <- sd(dat$pain_score)
joint_path <- data.table(
  n = nrow(dat), events = sum(dat$status), n_mediators = 24L,
  pain_score_sd = sd_pain,
  total_standardized_logHR = joint$total_logHR * sd_pain,
  direct_standardized_logHR = joint$direct_logHR * sd_pain,
  indirect_standardized_logHR = joint$indirect_logHR * sd_pain,
  note = "Cox coefficients are per 1-SD pain score. a/b coefficients are descriptive adjusted associations, not causal paths."
)
if (max(abs(paths$path_c_standardized_beta_logHR - joint_path$total_standardized_logHR)) > 1e-8) stop("Standardized total coefficient mismatch")
fwrite(joint_path, file.path(OUT, "joint24_standardized_path_coefficients.csv"))
fwrite(data.table(path = source_path, md5 = unname(tools::md5sum(source_path))), file.path(OUT, "path_source_hash.csv"))
print(joint_path)
