source("00_setup/ukb_helpers.R")
load_packages(c("data.table", "survival", "rms"))

analysis_path <- result_path("ukb_main_analysis_dataset.csv")
if (!file.exists(analysis_path)) stop("Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first.")

dat <- data.table::fread(analysis_path)
site_cols <- c("headache", "facial", "neck_shoulder", "back", "stomach", "hip", "knee")
covariate_terms <- paste(MAIN_COVARIATES, collapse = " + ")

for (site in site_cols) {
  dat[[paste0(site, "_acute")]] <- as.integer(dat[[site]] == 1)
  dat[[paste0(site, "_chronic")]] <- as.integer(dat[[site]] == 2)
}
pain_terms <- as.vector(rbind(paste0(site_cols, "_acute"), paste0(site_cols, "_chronic")))

mutual_formula <- as.formula(paste("Surv(time, status) ~", paste(c(pain_terms, MAIN_COVARIATES), collapse = " + ")))
mutual_fit <- survival::coxph(mutual_formula, data = dat)
beta <- stats::coef(mutual_fit)[pain_terms]
beta[is.na(beta)] <- 0
denominator <- sum(beta)
if (!is.finite(denominator) || denominator == 0) stop("Cannot build composite pain score: sum of pain-term coefficients is zero.")

dat$pain_score <- as.numeric(as.matrix(dat[, pain_terms, with = FALSE]) %*% beta) / denominator
dat$pain_score[dat$any_site_count == 0] <- 0

positive_scores <- dat$pain_score[dat$pain_score > 0]
tertiles <- quantile(positive_scores, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
dat$weighted_pain_status <- 0
dat$weighted_pain_status[dat$pain_score > 0 & dat$pain_score <= tertiles[1]] <- 1
dat$weighted_pain_status[dat$pain_score > tertiles[1] & dat$pain_score <= tertiles[2]] <- 2
dat$weighted_pain_status[dat$pain_score > tertiles[2]] <- 3
dat$weighted_pain_status <- relevel(factor(dat$weighted_pain_status), ref = "0")

score_fit <- survival::coxph(as.formula(paste("Surv(time, status) ~ weighted_pain_status +", covariate_terms)), data = dat)
score_summary <- summary(score_fit)
score_terms <- grep("^weighted_pain_status", rownames(score_summary$coefficients), value = TRUE)
score_result <- data.frame(
  contrast = score_terms,
  HR = score_summary$conf.int[score_terms, "exp(coef)"],
  lower = score_summary$conf.int[score_terms, "lower .95"],
  upper = score_summary$conf.int[score_terms, "upper .95"],
  p_value = score_summary$coefficients[score_terms, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

dd <- rms::datadist(dat)
options(datadist = "dd")
rcs_fit <- rms::cph(as.formula(paste("Surv(time, status) ~ rcs(pain_score, 3) +", covariate_terms)), data = dat, x = TRUE, y = TRUE)
rcs_anova <- as.data.frame(anova(rcs_fit))

data.table::fwrite(dat[, c("eid", "pain_score", "weighted_pain_status"), with = FALSE], COMPOSITE_SCORE_FILE)
data.table::fwrite(as.data.table(score_result), result_path("ukb_composite_pain_score_tertile_cox.csv"))
data.table::fwrite(as.data.table(rcs_anova), result_path("ukb_composite_pain_score_rcs_anova.csv"))

cat("Composite score output:", COMPOSITE_SCORE_FILE, "\n")
cat("Tertile Cox output:", result_path("ukb_composite_pain_score_tertile_cox.csv"), "\n")
cat("RCS output:", result_path("ukb_composite_pain_score_rcs_anova.csv"), "\n")

