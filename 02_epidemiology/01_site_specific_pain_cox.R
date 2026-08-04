source("00_setup/ukb_helpers.R")
load_packages(c("data.table", "survival"))

analysis_path <- result_path("ukb_main_analysis_dataset.csv")
if (!file.exists(analysis_path)) stop("Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first.")

dat <- data.table::fread(analysis_path)
site_cols <- c("headache", "facial", "neck_shoulder", "back", "stomach", "hip", "knee")
covariates <- MAIN_COVARIATES
covariate_terms <- paste(covariates, collapse = " + ")

extract_pain_term <- function(fit, term, site, duration_label) {
  s <- summary(fit)
  if (!(term %in% rownames(s$coefficients))) {
    return(data.frame(site = site, duration = duration_label, HR = NA, lower = NA, upper = NA, p_value = NA))
  }
  data.frame(
    site = site,
    duration = duration_label,
    HR = s$conf.int[term, "exp(coef)"],
    lower = s$conf.int[term, "lower .95"],
    upper = s$conf.int[term, "upper .95"],
    p_value = s$coefficients[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

results <- list()
for (site in site_cols) {
  model_data <- as.data.frame(dat[, c("time", "status", covariates, site), with = FALSE])
  names(model_data)[ncol(model_data)] <- "pain_site"
  model_data$pain_site <- relevel(factor(model_data$pain_site), ref = "0")
  fit <- survival::coxph(as.formula(paste("Surv(time, status) ~ pain_site +", covariate_terms)), data = model_data)
  results[[length(results) + 1]] <- extract_pain_term(fit, "pain_site1", site, "acute")
  results[[length(results) + 1]] <- extract_pain_term(fit, "pain_site2", site, "chronic")
}

result <- data.table::rbindlist(results, fill = TRUE)
result[, fdr := p.adjust(p_value, method = "BH")]
data.table::fwrite(result, result_path("ukb_site_specific_pain_cox.csv"))
print(result)

