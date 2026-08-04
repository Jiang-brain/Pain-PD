source("00_setup/config.R")

expand_multi_select <- function(x) {
  raw_values <- as.character(x)
  raw_values[is.na(raw_values)] <- ""
  split_values <- strsplit(raw_values, "\\|")
  max_len <- max(lengths(split_values))
  expanded <- matrix(NA_real_, nrow = length(raw_values), ncol = max_len)
  for (i in seq_along(split_values)) {
    values <- suppressWarnings(as.numeric(split_values[[i]]))
    if (length(values) > 0) expanded[i, seq_along(values)] <- values
  }
  as.data.frame(expanded)
}

get_multi_select <- function(dat, field_id) {
  cols <- grep(paste0("^X", field_id, "\\.0\\."), names(dat), value = TRUE)
  if (length(cols) == 1) return(expand_multi_select(dat[[cols]]))
  as.data.frame(dat[, cols, drop = FALSE])
}

make_censoring_data <- function(path = UKB_DATA_FILE) {
  load_packages(c("data.table"))
  first_diagnosis <- data.frame(data.table::fread(
    path,
    select = c("eid", "53-0.0", "42030-0.0", "42031-0.0", "42032-0.0", "42033-0.0")
  ))
  diagnosis_parkinson <- data.frame(
    id = first_diagnosis$eid,
    diagnosis_date = as.Date(NA),
    source = NA,
    time_interval = NA,
    attending_date = first_diagnosis$X53.0.0
  )
  diagnosis_parkinson$diagnosis_date <- pmin(
    first_diagnosis$X42030.0.0,
    first_diagnosis$X42032.0.0,
    na.rm = TRUE
  )
  index_0 <- which(first_diagnosis$X42030.0.0 == diagnosis_parkinson$diagnosis_date)
  diagnosis_parkinson$source[index_0] <- first_diagnosis$X42031.0.0[index_0]
  diagnosis_parkinson$time_interval <- diagnosis_parkinson$diagnosis_date - diagnosis_parkinson$attending_date
  diagnosis_parkinson <- diagnosis_parkinson[!is.na(diagnosis_parkinson$diagnosis_date), ]
  censoring_date <- max(diagnosis_parkinson$diagnosis_date, na.rm = TRUE)

  lost_follow_up <- data.frame(data.table::fread(path, select = c("eid", "190-0.0", "191-0.0", "53-0.0")))
  colnames(lost_follow_up) <- c("id", "lost_reason", "lost_date", "attending_date")
  lost_follow_up <- lost_follow_up[!is.na(lost_follow_up$lost_date), ]
  lost_follow_up$lost_time_interval <- lost_follow_up$lost_date - lost_follow_up$attending_date

  death_information <- data.table::fread(path, select = c("eid", "40000-0.0", "53-0.0"))
  colnames(death_information) <- c("id", "death_date", "attending_date")
  death_information <- death_information[!is.na(death_information$death_date)]
  death_information$time_interval <- death_information$death_date - death_information$attending_date

  datas <- data.frame(data.table::fread(path, select = c("eid", "53-0.0")))
  data_censoring <- data.frame(id = datas$eid, attending_date = datas$X53.0.0)
  data_censoring$censoring <- rep(data.table::as.IDate(censoring_date), nrow(data_censoring)) - data_censoring$attending_date
  data_censoring$lost <- NA
  data_censoring$lost_reason <- NA
  lost_index <- match(lost_follow_up$id, data_censoring$id)
  data_censoring$lost[lost_index] <- lost_follow_up$lost_time_interval
  data_censoring$lost_reason[lost_index] <- lost_follow_up$lost_reason
  data_censoring$death <- NA
  death_index <- match(death_information$id, data_censoring$id)
  data_censoring$death[death_index] <- death_information$time_interval
  data_censoring$diagnosis <- NA
  data_censoring$source <- NA
  diagnosis_index <- match(diagnosis_parkinson$id, data_censoring$id)
  data_censoring$diagnosis[diagnosis_index] <- diagnosis_parkinson$time_interval
  data_censoring$source[diagnosis_index] <- diagnosis_parkinson$source
  data_censoring$label <- 0
  index_diagnosis <- which(!is.na(data_censoring$diagnosis))
  data_censoring$label[index_diagnosis] <- 1
  data_censoring$censoring[index_diagnosis] <- data_censoring$diagnosis[index_diagnosis]
  index_lost_wo_diagnosis <- which(!is.na(data_censoring$lost) & data_censoring$label == 0)
  data_censoring$censoring[index_lost_wo_diagnosis] <- data_censoring$lost[index_lost_wo_diagnosis]
  index_death_wo_diagnosis <- which(!is.na(data_censoring$death) & data_censoring$label == 0)
  data_censoring$censoring[index_death_wo_diagnosis] <- data_censoring$death[index_death_wo_diagnosis]
  data_censoring
}

make_covariates <- function(path = UKB_DATA_FILE, n_pcs = 10) {
  load_packages(c("data.table"))
  covariate_cols <- c(
    "eid", "21003-0.0", "31-0.0", "21000-0.0", "21000-1.0", "21000-2.0",
    "22189-0.0", "738-0.0", "6138-0.0", "20116-0.0", "20117-0.0",
    "1070-0.0", "21001-0.0", "26260-0.0", paste0("22009-0.", 1:max(10, n_pcs))
  )
  datas <- data.frame(data.table::fread(path, select = covariate_cols))
  age <- datas$X21003.0.0
  gender <- datas$X31.0.0

  ethnicity <- cbind(datas$X21000.0.0, datas$X21000.1.0, datas$X21000.2.0)
  index_nan <- which(ethnicity[, 1] == -1 | ethnicity[, 1] == -3)
  ethnicity[index_nan, 1] <- pmax(ethnicity[index_nan, 1], ethnicity[index_nan, 2], ethnicity[index_nan, 3], na.rm = TRUE)
  race <- ethnicity[, 1]
  race[race %in% c(1, 1001, 1002, 1003)] <- 1
  race[!(race %in% c(1, 1001, 1002, 1003))] <- 0
  race[race %in% c(-1, -3)] <- NA

  deprivation <- datas$X22189.0.0
  tertile_point <- quantile(deprivation, probs = seq(0, 1, 1 / 3), na.rm = TRUE)
  deprivation[deprivation < tertile_point[2]] <- 0
  deprivation[deprivation >= tertile_point[2] & deprivation < tertile_point[3]] <- 1
  deprivation[deprivation >= tertile_point[3]] <- 2

  income <- datas$X738.0.0
  income[income == -1 | income == -3] <- 0
  income[income <= 3 & income >= 1] <- 1
  income[income == 4] <- 2
  income[income == 5] <- 3

  education <- get_multi_select(datas, "6138")
  for (i in seq_len(ncol(education))) {
    education[education[, i] == -3, i] <- NA
    education[education[, i] == -7, i] <- 7
    education[education[, i] == 1, i] <- 20
    education[education[, i] == 2, i] <- 13
    education[education[, i] == 3 | education[, i] == 4, i] <- 10
    education[education[, i] == 5, i] <- 19
    education[education[, i] == 6, i] <- 15
  }
  university_education <- do.call(pmax, c(education, list(na.rm = TRUE)))
  university_education[university_education == 7] <- 0
  university_education[university_education < 20 & university_education != 0] <- 1
  university_education[university_education == 20] <- 2

  smoking <- datas$X20116.0.0
  smoking[smoking == -3] <- NA
  smoking[smoking == 2] <- 1
  alcohol <- datas$X20117.0.0
  alcohol[alcohol < 0] <- NA
  alcohol[alcohol == 2] <- 1
  sedentary <- datas$X1070.0.0
  sedentary[sedentary == -1 | sedentary == -3] <- NA
  sedentary[(sedentary < 4 & sedentary >= 0) | sedentary == -10] <- 0
  sedentary[sedentary >= 4] <- 1
  bmi <- datas$X21001.0.0
  obesity <- bmi
  obesity[bmi < 18.5] <- NA
  obesity[bmi >= 18.5 & bmi < 25] <- 0
  obesity[bmi >= 25 & bmi < 30] <- 1
  obesity[bmi >= 30 & bmi < 35] <- 2
  obesity[bmi >= 35] <- 3

  pd_prs <- datas$X26260.0.0
  tertile_point <- quantile(pd_prs, probs = seq(0, 1, 1 / 3), na.rm = TRUE)
  pd_prs[pd_prs < tertile_point[2]] <- 0
  pd_prs[pd_prs >= tertile_point[2] & pd_prs < tertile_point[3]] <- 1
  pd_prs[pd_prs >= tertile_point[3]] <- 2

  covars <- data.frame(
    eid = datas$eid, age = age, gender = factor(gender), ethnicity = factor(race),
    deprivation = factor(deprivation), income = factor(income),
    university_education = factor(university_education), smoking = factor(smoking),
    alcohol = factor(alcohol), sedentary = factor(sedentary), obesity = factor(obesity),
    prs = factor(pd_prs)
  )
  for (k in seq_len(n_pcs)) covars[[paste0("comp_", k)]] <- datas[[paste0("X22009.0.", k)]]
  covars
}

make_pain_duration <- function(path = UKB_DATA_FILE) {
  load_packages(c("data.table"))
  dat <- data.frame(data.table::fread(path))
  any_pain <- dat$X6159.0.0
  any_pain[any_pain == -7] <- 0
  any_pain[any_pain > 0 | any_pain == -3] <- NA
  pain_duration <- data.frame(
    headache = dat$X3799.0.0,
    facial = dat$X4067.0.0,
    neck_shoulder = dat$X3404.0.0,
    back = dat$X3571.0.0,
    stomach = dat$X3741.0.0,
    hip = dat$X3414.0.0,
    knee = dat$X3773.0.0,
    all_over_body = dat$X2956.0.0
  )
  for (i in seq_len(ncol(pain_duration))) {
    chronic_pain <- which(pain_duration[, i] == 1)
    short_pain <- which(pain_duration[, i] == 0)
    not_answer <- which(pain_duration[, i] %in% c(-1, -3))
    no_pain <- which(any_pain == 0)
    pain_duration[chronic_pain, i] <- 2
    pain_duration[short_pain, i] <- 1
    pain_duration[not_answer, i] <- NA
    pain_duration[no_pain, i] <- 0
  }
  data.frame(eid = dat$eid, pain_duration)
}

drop_body_pain_without_sites <- function(dat) {
  if (!("all_over_body" %in% names(dat))) return(dat)
  dat[is.na(dat$all_over_body) | dat$all_over_body == 0, ]
}

make_outcome_data <- function(path = UKB_DATA_FILE) {
  censoring <- make_censoring_data(path)
  data.frame(eid = censoring$id, status = censoring$label, time = as.numeric(censoring$censoring))
}

apply_main_exclusions <- function(dat, followup_exclusion_days = FOLLOWUP_EXCLUSION_DAYS) {
  dat <- dat[complete.cases(dat), ]
  dat <- dat[dat$age >= 40, ]
  dat <- dat[!(dat$time <= followup_exclusion_days & dat$status == 1), ]
  droplevels(dat)
}
