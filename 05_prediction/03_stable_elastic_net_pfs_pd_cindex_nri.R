source("00_setup/ukb_helpers.R")
source("00_setup/model_helpers.R")
load_packages(c("data.table", "survival"))

score_path <- result_path("pain_related_proteomic_score_stable_elastic_net_pfs.csv")
coef_path <- result_path("pain_related_proteomic_score_stable_elastic_net_coefficients.csv")
cindex_out_path <- result_path("pd_incidence_stable_elastic_net_pfs_cindex.csv")
nri_out_path <- result_path("pd_incidence_stable_elastic_net_pfs_nri.csv")
nri_summary_path <- result_path("pd_incidence_stable_elastic_net_pfs_nri_summary.json")

if (!file.exists(score_path)) stop("Run 05_prediction/02_build_stable_nested_elastic_net_pfs.py first: ", score_path)
if (!file.exists(coef_path)) stop("Run 05_prediction/02_build_stable_nested_elastic_net_pfs.py first: ", coef_path)

outcome_data <- make_outcome_data(UKB_DATA_FILE)
covars <- make_covariates(UKB_DATA_FILE, n_pcs = 10)
protein_cov <- data.table::fread(PROTEIN_COVARIATE_FILE)
data.table::setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- data.table::fread(PROTEIN_TIME_FILE)
time_cols <- setdiff(names(protein_times), "eid")
protein_times[, protein_time := rowMeans(.SD, na.rm = TRUE), .SDcols = time_cols]
protein_times <- protein_times[, .(eid, protein_time)]

score_dt <- data.table::fread(score_path)
stable_dt <- data.table::fread(coef_path)
stable_dt <- stable_dt[method == "elastic_net" & stable_all_10_folds == TRUE]
if (nrow(stable_dt) == 0) stop("No stable Elastic Net proteins found in ", coef_path)

protein_raw <- data.table::fread(PROTEIN_MATRIX_FILE, select = c("eid", "ins_index", stable_dt$protein))
protein_raw <- protein_raw[order(eid, ins_index)]
protein_raw <- protein_raw[!duplicated(eid)]
protein_raw[, ins_index := NULL]
protein_term_names <- make.unique(make.names(paste0("protein_", stable_dt$protein)))
for (i in seq_along(stable_dt$protein)) {
  protein <- stable_dt$protein[i]
  values <- protein_raw[[protein]]
  values[is.na(values)] <- stable_dt$mean_impute_value[i]
  protein_raw[[protein_term_names[i]]] <- as.numeric(scale(values))
}
protein_raw <- protein_raw[, c("eid", protein_term_names), with = FALSE]

analysis_data <- merge(outcome_data, covars, by = "eid")
analysis_data <- apply_main_exclusions(analysis_data)
analysis_data <- merge(as.data.table(analysis_data), protein_cov, by = "eid")
analysis_data <- merge(analysis_data, protein_times, by = "eid")
analysis_data <- merge(analysis_data, score_dt[, .(eid, stable_elastic_net_pfs_z)], by = "eid")
analysis_data <- merge(analysis_data, protein_raw, by = "eid")

covariates <- c(MAIN_COVARIATES, PROTEOMIC_TECHNICAL_COVARIATES)
model_cols <- c("time", "status", covariates, "stable_elastic_net_pfs_z", protein_term_names)
model_data <- as.data.frame(analysis_data[, model_cols, with = FALSE])
model_data <- model_data[complete.cases(model_data), ]
model_data <- droplevels(model_data)
covs <- usable_covariates(model_data, covariates)
cov_terms <- paste(covs, collapse = " + ")
protein_terms <- paste(protein_term_names, collapse = " + ")

base_fit <- survival::coxph(as.formula(paste("Surv(time, status) ~", cov_terms)), data = model_data)
pfs_plus_fit <- survival::coxph(as.formula(paste("Surv(time, status) ~", cov_terms, "+ stable_elastic_net_pfs_z")), data = model_data)
protein_plus_fit <- survival::coxph(
  as.formula(paste("Surv(time, status) ~", cov_terms, "+", protein_terms)),
  data = model_data,
  control = survival::coxph.control(iter.max = 30)
)
protein_alone_fit <- survival::coxph(as.formula(paste("Surv(time, status) ~", protein_terms)), data = model_data)

cindex_result <- data.table::rbindlist(list(
  cox_fit_row("traditional_covariates", base_fit, NULL, "traditional_covariates"),
  cox_fit_row("stable_elastic_net_pfs_alone", survival::coxph(Surv(time, status) ~ stable_elastic_net_pfs_z, data = model_data), NULL, "stable_elastic_net_pfs_z"),
  cox_fit_row("traditional_plus_stable_elastic_net_pfs", pfs_plus_fit, base_fit, "stable_elastic_net_pfs_z"),
  cox_fit_row("stable_elastic_net_proteins_alone", protein_alone_fit, NULL, "stable_elastic_net_protein_panel"),
  cox_fit_row("traditional_plus_stable_elastic_net_proteins", protein_plus_fit, base_fit, "stable_elastic_net_protein_panel")
), fill = TRUE)

eval_time <- max(model_data$time, na.rm = TRUE)
nri_result <- data.table::rbindlist(list(
  continuous_nri_row(
    "stable_elastic_net",
    "traditional_covariates_vs_traditional_plus_stable_elastic_net_pfs",
    base_fit,
    pfs_plus_fit,
    model_data,
    eval_time,
    status_col = "status"
  ),
  continuous_nri_row(
    "stable_elastic_net",
    "traditional_covariates_vs_traditional_plus_stable_elastic_net_proteins",
    base_fit,
    protein_plus_fit,
    model_data,
    eval_time,
    status_col = "status"
  )
), fill = TRUE)

data.table::fwrite(cindex_result, cindex_out_path)
data.table::fwrite(nri_result, nri_out_path)
write_json_simple <- function(x, path) {
  json <- jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null")
  writeLines(json, path)
}
if (requireNamespace("jsonlite", quietly = TRUE)) {
  write_json_simple(list(
    method = "C-index and continuous NRI for stable nested-CV Elastic Net PPS and stable protein panel.",
    output_cindex = cindex_out_path,
    output_nri = nri_out_path,
    n_stable_proteins = nrow(stable_dt),
    inference = "Continuous NRI uses analytic normal-approximation CI from reclassification proportions; not bootstrap."
  ), nri_summary_path)
}

print(cindex_result)
print(nri_result)
cat("C-index output:", cindex_out_path, "\n")
cat("NRI output:", nri_out_path, "\n")

