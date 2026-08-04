source("00_setup/ukb_helpers.R")
load_packages(c("data.table", "survival"))

analysis_path <- result_path("ukb_main_analysis_dataset.csv")
if (!file.exists(analysis_path)) stop("Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first.")

dat <- data.table::fread(analysis_path)
covariate_terms <- paste(MAIN_COVARIATES, collapse = " + ")

fit_status <- function(status_col, label) {
  model_data <- as.data.frame(dat[, c("time", "status", MAIN_COVARIATES, status_col), with = FALSE])
  names(model_data)[ncol(model_data)] <- "pain_status"
  model_data$pain_status <- relevel(factor(model_data$pain_status), ref = "0")
  fit <- survival::coxph(as.formula(paste("Surv(time, status) ~ pain_status +", covariate_terms)), data = model_data)
  s <- summary(fit)
  terms <- intersect(c("pain_status1", "pain_status2"), rownames(s$coefficients))
  data.frame(
    exposure = label,
    contrast = c("single-site vs pain-free", "multisite vs pain-free")[match(terms, c("pain_status1", "pain_status2"))],
    HR = s$conf.int[terms, "exp(coef)"],
    lower = s$conf.int[terms, "lower .95"],
    upper = s$conf.int[terms, "upper .95"],
    p_value = s$coefficients[terms, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

result <- data.table::rbindlist(list(
  fit_status("acute_pain_status", "acute pain"),
  fit_status("chronic_pain_status", "chronic pain")
), fill = TRUE)
result[, fdr := p.adjust(p_value, method = "BH")]
data.table::fwrite(result, result_path("ukb_single_multisite_pain_cox.csv"))
print(result)

