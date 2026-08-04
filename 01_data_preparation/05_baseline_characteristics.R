source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

analysis_path <- result_path("ukb_main_analysis_dataset.csv")
if (!file.exists(analysis_path)) stop("Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first.")

dat <- data.table::fread(analysis_path)
dat$any_pain <- factor(ifelse(dat$any_site_count > 0, "pain", "pain_free"))

summarize_numeric <- function(x) {
  sprintf("%.2f (%.2f)", mean(x, na.rm = TRUE), stats::sd(x, na.rm = TRUE))
}

summarize_factor <- function(x) {
  tab <- table(x, useNA = "ifany")
  paste(sprintf("%s: %d (%.2f%%)", names(tab), as.integer(tab), 100 * as.integer(tab) / sum(tab)), collapse = "; ")
}

vars <- c("age", "gender", "ethnicity", "deprivation", "income", "university_education",
          "smoking", "alcohol", "sedentary", "obesity", "prs")
groups <- c("overall", levels(dat$any_pain))

rows <- list()
for (var in vars) {
  for (group in groups) {
    sub <- if (group == "overall") dat else dat[any_pain == group]
    value <- if (is.numeric(sub[[var]])) summarize_numeric(sub[[var]]) else summarize_factor(sub[[var]])
    rows[[length(rows) + 1]] <- data.frame(variable = var, group = group, value = value, stringsAsFactors = FALSE)
  }
}

summary_table <- data.table::rbindlist(rows)
data.table::fwrite(summary_table, result_path("table1_baseline_characteristics_by_pain.csv"))
cat("Output:", result_path("table1_baseline_characteristics_by_pain.csv"), "\n")

