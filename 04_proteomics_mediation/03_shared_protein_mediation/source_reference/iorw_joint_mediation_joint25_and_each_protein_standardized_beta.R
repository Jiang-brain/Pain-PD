find_input <- function(...) {
  candidates <- c(...)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)==0) {
    stop("None of these input files exist: ", paste(candidates, collapse=", "))
  }
  existing[1]
}

data_path <- find_input("PD_data.csv", "../PA/PD_data.csv", "../PA2/PD_data.csv")
protein_path <- "protein_new_missing60_filtered.csv"
pain_score_path <- "composite_pain_scores.csv"
intersection_path <- "protein_intersection_FDR05.csv"
protein_cov_path <- find_input("PD_data_matched_protein_covariates.csv", "../PA/PD_data_matched_protein_covariates.csv", "../PA2/PD_data_matched_protein_covariates.csv")
protein_time_path <- find_input("PD_data_matched_protein_times.csv", "../PA/PD_data_matched_protein_times.csv", "../PA2/PD_data_matched_protein_times.csv")

out_data_path <- "overlap25_iorw_joint_mediation_analysis_data.csv"
out_result_path <- "overlap25_iorw_joint_mediation_results.csv"
out_boot_path <- "overlap25_iorw_joint_mediation_bootstrap.csv"
out_summary_path <- "overlap25_iorw_joint_mediation_summary.json"
out_joint25_path_abc_path <- "overlap25_iorw_joint25_path_abc_beta.csv"

# Standardized beta outputs:
# - Cox beta/logHR is per 1 SD increase in the predictor.
# - path a standardized beta is SD change in joint mediator score per 1 SD pain_score.
out_std_result_path <- "overlap25_iorw_joint_mediation_results_standardized_beta.csv"
out_joint25_std_path_abc_path <- "overlap25_iorw_joint25_path_abc_standardized_beta.csv"
out_each_protein_std_path_abc_path <- "overlap25_iorw_each_protein_path_abc_standardized_beta.csv"

.libPaths(c(normalizePath("R_libs"), .libPaths()))
library(data.table)
library(survival)

set.seed(20260602)

n_boot <- as.integer(Sys.getenv("IORW_N_BOOT", "0"))
weight_trim_quantile <- as.numeric(Sys.getenv("IORW_WEIGHT_TRIM_Q", "0.99"))

expand_multi_select <- function(x) {
  raw_values <- as.character(x)
  raw_values[is.na(raw_values)] <- ""
  split_values <- strsplit(raw_values, "\\|")
  max_len <- max(lengths(split_values))
  expanded <- matrix(NA_real_, nrow=length(raw_values), ncol=max_len)
  for (i in seq_along(split_values)) {
    values <- suppressWarnings(as.numeric(split_values[[i]]))
    if (length(values)>0) expanded[i,seq_along(values)] <- values
  }
  as.data.frame(expanded)
}

get_multi_select <- function(dat, field_id) {
  cols <- grep(paste0("^X",field_id,"\\.0\\."), names(dat), value=TRUE)
  if (length(cols)==1) return(expand_multi_select(dat[[cols]]))
  as.data.frame(dat[,cols,drop=FALSE])
}

make_censoring_data <- function(path) {
  first_diagnosis <- data.frame(fread(path, select=c("eid","53-0.0","42030-0.0","42031-0.0","42032-0.0","42033-0.0")))
  diagnosis_parkinson <- data.frame(
    id=first_diagnosis$eid,
    diagnosis_date=as.Date(NA),
    source=NA,
    time_interval=NA,
    attending_date=first_diagnosis$X53.0.0
  )
  diagnosis_parkinson$diagnosis_date <- pmin(first_diagnosis$X42030.0.0, first_diagnosis$X42032.0.0, na.rm=TRUE)
  index_0 <- which(first_diagnosis$X42030.0.0==diagnosis_parkinson$diagnosis_date)
  diagnosis_parkinson$source[index_0] <- first_diagnosis$X42031.0.0[index_0]
  diagnosis_parkinson$time_interval <- diagnosis_parkinson$diagnosis_date-diagnosis_parkinson$attending_date
  diagnosis_parkinson <- diagnosis_parkinson[which(!is.na(diagnosis_parkinson$diagnosis_date)),]
  censoring_date <- max(diagnosis_parkinson$diagnosis_date, na.rm=TRUE)
  
  lost_follow_up <- data.frame(fread(path, select=c("eid","190-0.0","191-0.0","53-0.0")))
  colnames(lost_follow_up) <- c("id","lost_reason","lost_date","attending_date")
  lost_follow_up <- lost_follow_up[which(!is.na(lost_follow_up$lost_date)),]
  lost_follow_up$lost_time_interval <- lost_follow_up$lost_date-lost_follow_up$attending_date
  
  death_information <- fread(path, select=c("eid","40000-0.0","53-0.0"))
  colnames(death_information) <- c("id","death_date","attending_date")
  death_information <- death_information[which(!is.na(death_information$death_date))]
  death_information$time_interval <- death_information$death_date-death_information$attending_date
  
  datas <- data.frame(fread(path, select=c("eid","53-0.0")))
  data_censoring <- data.frame(id=datas$eid, attending_date=datas$X53.0.0)
  data_censoring$censoring <- rep(as.IDate(censoring_date), nrow(data_censoring))-data_censoring$attending_date
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
  index_lost_wo_diagnosis <- which(!is.na(data_censoring$lost) & data_censoring$label==0)
  data_censoring$censoring[index_lost_wo_diagnosis] <- data_censoring$lost[index_lost_wo_diagnosis]
  index_death_wo_diagnosis <- which(!is.na(data_censoring$death) & data_censoring$label==0)
  data_censoring$censoring[index_death_wo_diagnosis] <- data_censoring$death[index_death_wo_diagnosis]
  data_censoring$label[index_death_wo_diagnosis] <- 0
  data_censoring
}

make_covariates <- function(path) {
  covariate_cols <- c(
    "eid","21003-0.0","31-0.0","21000-0.0","21000-1.0","21000-2.0",
    "22189-0.0","738-0.0","6138-0.0","20116-0.0","20117-0.0",
    "1070-0.0","21001-0.0","26260-0.0",paste0("22009-0.",1:20)
  )
  datas <- data.frame(fread(path, select=covariate_cols))
  age <- datas$X21003.0.0
  gender <- datas$X31.0.0
  ethnicity <- cbind(datas$X21000.0.0,datas$X21000.1.0,datas$X21000.2.0)
  index_nan <- which(ethnicity[,1]==-1|ethnicity[,1]==-3)
  ethnicity[index_nan,1] <- pmax(ethnicity[index_nan,1],ethnicity[index_nan,2],ethnicity[index_nan,3],na.rm=TRUE)
  race <- ethnicity[,1]
  race[which(race %in% c(1,1001,1002,1003))] <- 1
  race[which((race %in% c(1,1001,1002,1003))==FALSE)] <- 0
  race[which(race==-1|race==-3)] <- NA
  
  deprivation <- datas$X22189.0.0
  tertile_point <- quantile(deprivation, probs=seq(0,1,1/3), na.rm=TRUE)
  deprivation[which(deprivation<tertile_point[2])] <- 0
  deprivation[which(deprivation>=tertile_point[2]&deprivation<tertile_point[3])] <- 1
  deprivation[which(deprivation>=tertile_point[3])] <- 2
  
  income <- datas$X738.0.0
  income[which(income==-1)] <- 0
  income[which(income<=3&income>=1)] <- 1
  income[which(income==4)] <- 2
  income[which(income==5)] <- 3
  income[which(income==-3)] <- 0
  
  education <- get_multi_select(datas,"6138")
  for (i in seq_len(ncol(education))) {
    education[which(education[,i]==-3),i] <- NA
    education[which(education[,i]==-7),i] <- 7
    education[which(education[,i]==1),i] <- 20
    education[which(education[,i]==2),i] <- 13
    education[which(education[,i]==3|education[,i]==4),i] <- 10
    education[which(education[,i]==5),i] <- 19
    education[which(education[,i]==6),i] <- 15
  }
  university_education <- do.call(pmax,c(education,list(na.rm=TRUE)))
  university_education[which(university_education==7)] <- 0
  university_education[which(university_education<20&university_education!=0)] <- 1
  university_education[which(university_education==20)] <- 2
  
  smoking <- datas$X20116.0.0
  smoking[which(smoking==-3)] <- NA
  smoking[which(smoking==2)] <- 1
  alcohol <- datas$X20117.0.0
  alcohol[which(alcohol<0)] <- NA
  alcohol[which(alcohol==2)] <- 1
  sedentary <- datas$X1070.0.0
  sedentary[which(sedentary==-1|sedentary==-3)] <- NA
  sedentary[which((sedentary<4&sedentary>=0)|sedentary==-10)] <- 0
  sedentary[which(sedentary>=4)] <- 1
  bmi <- datas$X21001.0.0
  obesity <- bmi
  obesity[which(bmi<18.5)] <- NA
  obesity[which(bmi>=18.5&bmi<25)] <- 0
  obesity[which(bmi>=25&bmi<30)] <- 1
  obesity[which(bmi>=30&bmi<35)] <- 2
  obesity[which(bmi>=35)] <- 3
  pd_prs <- datas$X26260.0.0
  tertile_point <- quantile(pd_prs, probs=seq(0,1,1/3), na.rm=TRUE)
  pd_prs[which(pd_prs<tertile_point[2])] <- 0
  pd_prs[which(pd_prs>=tertile_point[2]&pd_prs<tertile_point[3])] <- 1
  pd_prs[which(pd_prs>=tertile_point[3])] <- 2
  
  covars <- data.frame(
    eid=datas$eid, age=age, gender=factor(gender), ethnicity=factor(race),
    deprivation=factor(deprivation), income=factor(income),
    university_education=factor(university_education), smoking=factor(smoking),
    alcohol=factor(alcohol), sedentary=factor(sedentary), obesity=factor(obesity),
    prs=factor(pd_prs)
  )
  for (k in 1:20) covars[[paste0("comp_",k)]] <- datas[[paste0("X22009.0.",k)]]
  covars
}

usable_covariates <- function(model_data, covariates) {
  keep <- vapply(covariates, function(v) {
    x <- model_data[[v]]
    if (is.factor(x) || is.character(x)) return(nlevels(droplevels(factor(x))) >= 2)
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))
  covariates[keep]
}


zscore_numeric <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm=TRUE)
  if (!is.finite(s) || s == 0) {
    stop("Cannot standardize variable because SD is zero or non-finite.")
  }
  as.numeric((x - mean(x, na.rm=TRUE)) / s)
}

make_model_matrix_terms <- function(dat, vars) {
  vars <- usable_covariates(dat, vars)
  terms <- paste(vars, collapse="+")
  list(vars=vars, terms=terms)
}

write_json_simple <- function(x, path) {
  scalar <- function(v) {
    if (length(v) == 1 && is.na(v)) return("null")
    if (is.character(v)) return(paste0('"', gsub('"', '\\"', gsub('\\\\', '/', v)), '"'))
    if (is.logical(v)) return(ifelse(v, "true", "false"))
    if (is.numeric(v)) return(as.character(v))
    paste0('"', as.character(v), '"')
  }
  lines <- c("{")
  nm <- names(x)
  for (i in seq_along(x)) {
    comma <- ifelse(i < length(x), ",", "")
    value <- x[[i]]
    if (length(value) == 1) lines <- c(lines, paste0('  "', nm[i], '": ', scalar(value), comma))
    else lines <- c(lines, paste0('  "', nm[i], '": [', paste(vapply(value, scalar, character(1)), collapse=", "), ']', comma))
  }
  writeLines(c(lines, "}"), path)
}

build_analysis_data <- function() {
  cat("selecting 25 same-direction FDR overlap proteins...\n")
  intersection <- fread(intersection_path)
  intersection[, same_direction := sign(beta_pain_score) == sign(beta_PD_cox)]
  overlap <- intersection[same_direction==TRUE]
  if (nrow(overlap) != 25) stop("Expected 25 same-direction proteins, got ", nrow(overlap))
  proteins <- overlap$protein
  direction <- sign(overlap$beta_pain_score)
  names(direction) <- proteins
  
  cat("building eligible analysis data...\n")
  data_censoring <- make_censoring_data(data_path)
  covars <- make_covariates(data_path)
  pain_scores <- fread(pain_score_path, select=c("eid","pain_score"))
  protein_cov <- fread(protein_cov_path)
  setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
  protein_cov[, season_final := factor(season_final)]
  protein_cov[, Batch := factor(Batch)]
  protein_times <- fread(protein_time_path)
  time_cols <- setdiff(names(protein_times), "eid")
  protein_times[, protein_time := rowMeans(.SD, na.rm=TRUE), .SDcols=time_cols]
  protein_times <- protein_times[, .(eid, protein_time)]
  
  outcome_data <- data.frame(eid=data_censoring$id, status=data_censoring$label, time=as.numeric(data_censoring$censoring))
  analysis_base <- merge(outcome_data, covars, by="eid")
  analysis_base <- merge(analysis_base, pain_scores, by="eid")
  analysis_base <- analysis_base[complete.cases(analysis_base),]
  analysis_base <- analysis_base[analysis_base$age>=40,]
  analysis_base <- analysis_base[!(analysis_base$time<365*2 & analysis_base$status==1),]
  analysis_base <- merge(as.data.table(analysis_base), protein_cov, by="eid")
  analysis_base <- merge(analysis_base, protein_times, by="eid")
  
  protein_dt <- fread(protein_path, select=c("eid","ins_index",proteins))
  setorder(protein_dt, eid, ins_index)
  protein_dt <- protein_dt[, .SD[1], by=eid]
  protein_dt[, ins_index := NULL]
  model_data <- merge(analysis_base, protein_dt, by="eid")
  
  protein_matrix <- as.matrix(model_data[, ..proteins])
  storage.mode(protein_matrix) <- "double"
  protein_median <- apply(protein_matrix, 2, median, na.rm=TRUE)
  for (j in seq_len(ncol(protein_matrix))) {
    idx <- is.na(protein_matrix[,j])
    if (any(idx)) protein_matrix[idx,j] <- protein_median[j]
  }
  protein_mean <- colMeans(protein_matrix)
  protein_sd <- apply(protein_matrix, 2, sd)
  protein_sd[protein_sd==0 | is.na(protein_sd)] <- 1
  protein_z <- scale(protein_matrix, center=protein_mean, scale=protein_sd)
  aligned_z <- sweep(protein_z, 2, direction[colnames(protein_z)], `*`)
  mediator_cols <- paste0("aligned_", proteins, "_z")
  for (j in seq_along(proteins)) model_data[[mediator_cols[j]]] <- as.numeric(aligned_z[,proteins[j]])
  
  list(
    data=model_data,
    proteins=proteins,
    mediator_cols=mediator_cols,
    overlap=overlap,
    protein_median=protein_median,
    protein_mean=protein_mean,
    protein_sd=protein_sd
  )
}

compute_iorw <- function(dat, covariates, mediator_cols, weight_trim_quantile=0.99) {
  dat <- as.data.frame(dat)
  dat <- droplevels(dat)
  covariates <- usable_covariates(dat, covariates)
  cov_terms <- paste(covariates, collapse="+")
  mediator_terms <- paste(mediator_cols, collapse="+")
  
  numerator_formula <- as.formula(paste("pain_score ~", cov_terms))
  denominator_formula <- as.formula(paste("pain_score ~", paste(c(covariates, mediator_cols), collapse="+")))
  total_formula <- as.formula(paste("Surv(time, status) ~ pain_score +", cov_terms))
  
  numerator_fit <- lm(numerator_formula, data=dat)
  denominator_fit <- lm(denominator_formula, data=dat)
  numerator_sd <- sqrt(sum(residuals(numerator_fit)^2) / numerator_fit$df.residual)
  denominator_sd <- sqrt(sum(residuals(denominator_fit)^2) / denominator_fit$df.residual)
  
  numerator_density <- dnorm(dat$pain_score, mean=predict(numerator_fit, newdata=dat), sd=numerator_sd)
  denominator_density <- dnorm(dat$pain_score, mean=predict(denominator_fit, newdata=dat), sd=denominator_sd)
  weights <- numerator_density / denominator_density
  weights[!is.finite(weights)] <- NA_real_
  trim_value <- quantile(weights, probs=weight_trim_quantile, na.rm=TRUE)
  weights_trimmed <- pmin(weights, trim_value)
  weights_trimmed <- weights_trimmed / mean(weights_trimmed, na.rm=TRUE)
  dat$iorw_weight <- weights_trimmed
  
  total_fit <- coxph(total_formula, data=dat)
  direct_fit <- coxph(total_formula, data=dat, weights=iorw_weight, robust=TRUE)
  naive_direct_fit <- coxph(
    as.formula(paste("Surv(time, status) ~ pain_score +", mediator_terms, "+", cov_terms)),
    data=dat
  )
  
  total_coef <- summary(total_fit)$coefficients["pain_score",]
  direct_coef <- summary(direct_fit)$coefficients["pain_score",]
  naive_direct_coef <- summary(naive_direct_fit)$coefficients["pain_score",]
  total <- total_coef["coef"]
  direct <- direct_coef["coef"]
  indirect <- total - direct
  
  result <- data.table(
    method="IORW_joint_25_direction_aligned_proteins",
    n=nrow(dat),
    events=sum(dat$status),
    n_mediators=length(mediator_cols),
    exposure="pain_score",
    mediator_set="25_direction_aligned_overlap_proteins",
    total_logHR=total,
    total_HR=exp(total),
    total_se=total_coef["se(coef)"],
    total_p_value=total_coef["Pr(>|z|)"],
    direct_logHR=direct,
    direct_HR=exp(direct),
    direct_robust_se=direct_coef["robust se"],
    direct_p_value=direct_coef["Pr(>|z|)"],
    indirect_logHR=indirect,
    indirect_HR=exp(indirect),
    proportion_mediated=indirect / total,
    naive_direct_logHR=naive_direct_coef["coef"],
    naive_direct_HR=exp(naive_direct_coef["coef"]),
    weight_mean=mean(weights_trimmed, na.rm=TRUE),
    weight_sd=sd(weights_trimmed, na.rm=TRUE),
    weight_min=min(weights_trimmed, na.rm=TRUE),
    weight_p50=quantile(weights_trimmed, 0.50, na.rm=TRUE),
    weight_p95=quantile(weights_trimmed, 0.95, na.rm=TRUE),
    weight_p99=quantile(weights_trimmed, 0.99, na.rm=TRUE),
    weight_max=max(weights_trimmed, na.rm=TRUE),
    weight_trim_quantile=weight_trim_quantile,
    weight_trim_value_raw=trim_value,
    numerator_R2=summary(numerator_fit)$r.squared,
    denominator_R2=summary(denominator_fit)$r.squared,
    c_index_total=unname(summary(total_fit)$concordance["C"]),
    c_index_direct=unname(summary(direct_fit)$concordance["C"])
  )
  list(result=result, weights=weights_trimmed, covariates=covariates)
}


compute_iorw_standardized <- function(dat, covariates, mediator_cols, weight_trim_quantile=0.99, exposure="pain_score_z") {
  # Same IORW model as compute_iorw(), but the exposure is standardized.
  # Resulting total/direct/indirect logHR values are beta values per 1 SD increase in exposure.
  dat <- as.data.frame(dat)
  dat <- droplevels(dat)
  covariates <- usable_covariates(dat, covariates)
  cov_terms <- paste(covariates, collapse="+")
  mediator_terms <- paste(mediator_cols, collapse="+")
  
  numerator_formula <- as.formula(paste(exposure, "~", cov_terms))
  denominator_formula <- as.formula(paste(exposure, "~", paste(c(covariates, mediator_cols), collapse="+")))
  total_formula <- as.formula(paste("Surv(time, status) ~", exposure, "+", cov_terms))
  
  numerator_fit <- lm(numerator_formula, data=dat)
  denominator_fit <- lm(denominator_formula, data=dat)
  numerator_sd <- sqrt(sum(residuals(numerator_fit)^2) / numerator_fit$df.residual)
  denominator_sd <- sqrt(sum(residuals(denominator_fit)^2) / denominator_fit$df.residual)
  
  numerator_density <- dnorm(dat[[exposure]], mean=predict(numerator_fit, newdata=dat), sd=numerator_sd)
  denominator_density <- dnorm(dat[[exposure]], mean=predict(denominator_fit, newdata=dat), sd=denominator_sd)
  weights <- numerator_density / denominator_density
  weights[!is.finite(weights)] <- NA_real_
  trim_value <- quantile(weights, probs=weight_trim_quantile, na.rm=TRUE)
  weights_trimmed <- pmin(weights, trim_value)
  weights_trimmed <- weights_trimmed / mean(weights_trimmed, na.rm=TRUE)
  dat$iorw_weight <- weights_trimmed
  
  total_fit <- coxph(total_formula, data=dat)
  direct_fit <- coxph(total_formula, data=dat, weights=iorw_weight, robust=TRUE)
  naive_direct_fit <- coxph(
    as.formula(paste("Surv(time, status) ~", exposure, "+", mediator_terms, "+", cov_terms)),
    data=dat
  )
  
  total_coef <- summary(total_fit)$coefficients[exposure,]
  direct_coef <- summary(direct_fit)$coefficients[exposure,]
  naive_direct_coef <- summary(naive_direct_fit)$coefficients[exposure,]
  total <- total_coef["coef"]
  direct <- direct_coef["coef"]
  indirect <- total - direct
  
  result <- data.table(
    method="IORW_joint_25_direction_aligned_proteins_standardized_exposure",
    beta_scale="standardized_beta_per_1_SD_exposure",
    n=nrow(dat),
    events=sum(dat$status),
    n_mediators=length(mediator_cols),
    exposure=exposure,
    mediator_set="25_direction_aligned_overlap_proteins",
    path_c_standardized_beta_logHR=total,
    path_c_standardized_HR_per_1SD=exp(total),
    path_c_standardized_se=total_coef["se(coef)"],
    path_c_standardized_p_value=total_coef["Pr(>|z|)"],
    path_cprime_standardized_beta_logHR=direct,
    path_cprime_standardized_HR_per_1SD=exp(direct),
    path_cprime_standardized_robust_se=direct_coef["robust se"],
    path_cprime_standardized_p_value=direct_coef["Pr(>|z|)"],
    joint_indirect_standardized_beta_logHR=indirect,
    joint_indirect_standardized_HR_per_1SD=exp(indirect),
    proportion_mediated=indirect / total,
    naive_cprime_standardized_beta_logHR=naive_direct_coef["coef"],
    naive_cprime_standardized_HR_per_1SD=exp(naive_direct_coef["coef"]),
    weight_mean=mean(weights_trimmed, na.rm=TRUE),
    weight_sd=sd(weights_trimmed, na.rm=TRUE),
    weight_min=min(weights_trimmed, na.rm=TRUE),
    weight_p50=quantile(weights_trimmed, 0.50, na.rm=TRUE),
    weight_p95=quantile(weights_trimmed, 0.95, na.rm=TRUE),
    weight_p99=quantile(weights_trimmed, 0.99, na.rm=TRUE),
    weight_max=max(weights_trimmed, na.rm=TRUE),
    weight_trim_quantile=weight_trim_quantile,
    weight_trim_value_raw=trim_value,
    numerator_R2=summary(numerator_fit)$r.squared,
    denominator_R2=summary(denominator_fit)$r.squared,
    c_index_total=unname(summary(total_fit)$concordance["C"]),
    c_index_direct=unname(summary(direct_fit)$concordance["C"])
  )
  list(result=result, weights=weights_trimmed, covariates=covariates)
}

compute_joint25_standardized_path_abc_beta <- function(dat, covariates, mediator_cols, iorw_std_result=NULL,
                                                       exposure_z="pain_score_z",
                                                       joint_mediator_z="joint25_mediator_score_z") {
  # One-row standardized joint 25-mediator path summary.
  # path a: standardized beta from lm(joint_mediator_z ~ exposure_z + covariates)
  # path b: Cox beta/logHR per 1 SD joint mediator score, adjusted for exposure_z + covariates
  # path c: Cox beta/logHR per 1 SD exposure
  # path c': Cox beta/logHR per 1 SD exposure after adding the joint mediator score
  dat <- as.data.frame(dat)
  dat <- droplevels(dat)
  
  needed_cols <- c("time", "status", exposure_z, joint_mediator_z, covariates)
  needed_cols <- unique(needed_cols[needed_cols %in% names(dat)])
  dat <- dat[complete.cases(dat[, needed_cols, drop=FALSE]),]
  dat <- droplevels(dat)
  
  covariates <- usable_covariates(dat, covariates)
  make_rhs <- function(vars) paste(vars, collapse="+")
  
  # path a: exposure_z -> joint_mediator_z
  a_formula <- as.formula(
    paste(joint_mediator_z, "~", make_rhs(c(exposure_z, covariates)))
  )
  a_fit <- lm(a_formula, data=dat)
  a_coef <- summary(a_fit)$coefficients[exposure_z,]
  
  # path c: exposure_z -> outcome, total effect
  c_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure_z, covariates)))
  )
  c_fit <- coxph(c_formula, data=dat)
  c_coef <- summary(c_fit)$coefficients[exposure_z,]
  
  # path b and traditional c': outcome ~ exposure_z + joint_mediator_z + covariates
  bc_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure_z, joint_mediator_z, covariates)))
  )
  bc_fit <- coxph(bc_formula, data=dat)
  bc_coef <- summary(bc_fit)$coefficients
  b_coef <- bc_coef[joint_mediator_z,]
  cprime_traditional_coef <- bc_coef[exposure_z,]
  
  get_iorw_std <- function(name, default=NA_real_) {
    if (is.null(iorw_std_result)) return(default)
    if (name %in% names(iorw_std_result)) return(iorw_std_result[[name]][1])
    default
  }
  
  out <- data.table(
    analysis_level="joint_25_mediators",
    beta_scale="standardized_beta",
    beta_interpretation="path a: SD mediator per 1 SD pain score; Cox paths: logHR per 1 SD predictor",
    mediator_set="25_direction_aligned_overlap_proteins",
    joint_mediator="joint25_mediator_score_z",
    joint_mediator_definition="z-scored rowMeans of the 25 direction-aligned z-scored protein mediators",
    exposure=exposure_z,
    n=nrow(dat),
    events=sum(dat$status),
    n_mediators=length(mediator_cols),
    
    path_a_standardized_beta=a_coef["Estimate"],
    path_a_standardized_se=a_coef["Std. Error"],
    path_a_standardized_t=a_coef["t value"],
    path_a_standardized_p_value=a_coef["Pr(>|t|)"],
    path_a_model="lm(joint25_mediator_score_z ~ pain_score_z + covariates)",
    
    path_b_standardized_beta_logHR=b_coef["coef"],
    path_b_standardized_HR_per_1SD=exp(b_coef["coef"]),
    path_b_standardized_se=b_coef["se(coef)"],
    path_b_standardized_z=b_coef["z"],
    path_b_standardized_p_value=b_coef["Pr(>|z|)"],
    path_b_model="coxph(Surv(time,status) ~ pain_score_z + joint25_mediator_score_z + covariates)",
    
    path_c_standardized_beta_logHR=c_coef["coef"],
    path_c_standardized_HR_per_1SD=exp(c_coef["coef"]),
    path_c_standardized_se=c_coef["se(coef)"],
    path_c_standardized_z=c_coef["z"],
    path_c_standardized_p_value=c_coef["Pr(>|z|)"],
    path_c_model="coxph(Surv(time,status) ~ pain_score_z + covariates)",
    
    path_cprime_standardized_beta_logHR_traditional_joint_score=cprime_traditional_coef["coef"],
    path_cprime_standardized_HR_per_1SD_traditional_joint_score=exp(cprime_traditional_coef["coef"]),
    path_cprime_standardized_se_traditional_joint_score=cprime_traditional_coef["se(coef)"],
    path_cprime_standardized_z_traditional_joint_score=cprime_traditional_coef["z"],
    path_cprime_standardized_p_value_traditional_joint_score=cprime_traditional_coef["Pr(>|z|)"],
    
    product_ab_standardized_beta_logHR=a_coef["Estimate"] * b_coef["coef"],
    difference_c_cprime_standardized_beta_logHR_traditional_joint_score=c_coef["coef"] - cprime_traditional_coef["coef"],
    
    path_c_standardized_beta_logHR_iorw=get_iorw_std("path_c_standardized_beta_logHR"),
    path_c_standardized_HR_per_1SD_iorw=get_iorw_std("path_c_standardized_HR_per_1SD"),
    
    path_cprime_standardized_beta_logHR_iorw=get_iorw_std("path_cprime_standardized_beta_logHR"),
    path_cprime_standardized_HR_per_1SD_iorw=get_iorw_std("path_cprime_standardized_HR_per_1SD"),
    
    joint_indirect_standardized_beta_logHR_iorw=get_iorw_std("joint_indirect_standardized_beta_logHR"),
    joint_indirect_standardized_HR_per_1SD_iorw=get_iorw_std("joint_indirect_standardized_HR_per_1SD"),
    proportion_mediated_iorw=get_iorw_std("proportion_mediated")
  )
  
  out
}



compute_each_protein_standardized_path_abc_beta <- function(dat, covariates, mediator_cols, proteins=NULL,
                                                            exposure_z="pain_score_z") {
  # Per-protein standardized path summary for the 25 direction-aligned proteins.
  # Each mediator is standardized in the final complete-case analysis sample.
  # path a: mediator_z ~ pain_score_z + covariates
  # path b single: Cox beta/logHR for mediator_z in Surv ~ pain_score_z + mediator_z + covariates
  # path b joint25-adjusted: Cox beta/logHR for mediator_z in Surv ~ pain_score_z + all 25 mediator_z + covariates
  # path c: Cox beta/logHR for pain_score_z in Surv ~ pain_score_z + covariates
  dat <- as.data.frame(dat)
  dat <- droplevels(dat)
  
  if (is.null(proteins)) {
    proteins <- gsub("^aligned_", "", mediator_cols)
    proteins <- gsub("_z$", "", proteins)
  }
  if (length(proteins) != length(mediator_cols)) {
    proteins <- mediator_cols
  }
  
  needed_cols <- c("time", "status", exposure_z, covariates, mediator_cols)
  needed_cols <- unique(needed_cols[needed_cols %in% names(dat)])
  dat <- dat[complete.cases(dat[, needed_cols, drop=FALSE]),]
  dat <- droplevels(dat)
  
  covariates <- usable_covariates(dat, covariates)
  make_rhs <- function(vars) paste(vars, collapse="+")
  get_val <- function(x, nm) {
    if (is.null(x) || length(x)==0 || !(nm %in% names(x))) return(NA_real_)
    as.numeric(x[nm])
  }
  get_row <- function(coef_mat, rn) {
    if (is.null(coef_mat) || !(rn %in% rownames(coef_mat))) return(NULL)
    coef_mat[rn,]
  }
  
  # Standardize each protein mediator in the final analytic sample.
  mediator_std_cols <- paste0("mediator_std_", seq_along(mediator_cols))
  for (i in seq_along(mediator_cols)) {
    dat[[mediator_std_cols[i]]] <- zscore_numeric(dat[[mediator_cols[i]]])
  }
  
  # path c: standardized exposure -> outcome, total effect. Same for all proteins.
  c_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure_z, covariates)))
  )
  c_fit <- coxph(c_formula, data=dat)
  c_coef <- summary(c_fit)$coefficients[exposure_z,]
  
  # Joint 25 adjusted model: gives one path b coefficient for each protein,
  # adjusted for pain_score_z, all other 24 proteins, and covariates.
  joint_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure_z, mediator_std_cols, covariates)))
  )
  joint_fit <- coxph(joint_formula, data=dat)
  joint_coef <- summary(joint_fit)$coefficients
  cprime_joint25_coef <- get_row(joint_coef, exposure_z)
  
  rows <- vector("list", length(mediator_cols))
  
  for (i in seq_along(mediator_cols)) {
    mediator <- mediator_cols[i]
    mediator_std <- mediator_std_cols[i]
    protein <- proteins[i]
    
    # path a: pain_score_z -> mediator_z
    a_formula <- as.formula(
      paste(mediator_std, "~", make_rhs(c(exposure_z, covariates)))
    )
    a_fit <- lm(a_formula, data=dat)
    a_coef <- summary(a_fit)$coefficients[exposure_z,]
    
    # Single-protein path b and c': outcome ~ pain_score_z + one mediator_z + covariates
    single_formula <- as.formula(
      paste("Surv(time, status) ~", make_rhs(c(exposure_z, mediator_std, covariates)))
    )
    single_fit <- coxph(single_formula, data=dat)
    single_coef <- summary(single_fit)$coefficients
    b_single_coef <- get_row(single_coef, mediator_std)
    cprime_single_coef <- get_row(single_coef, exposure_z)
    
    # Joint 25-adjusted path b for this protein.
    b_joint25_coef <- get_row(joint_coef, mediator_std)
    
    a_beta <- get_val(a_coef, "Estimate")
    b_single_beta <- get_val(b_single_coef, "coef")
    b_joint25_beta <- get_val(b_joint25_coef, "coef")
    c_beta <- get_val(c_coef, "coef")
    cprime_single_beta <- get_val(cprime_single_coef, "coef")
    cprime_joint25_beta <- get_val(cprime_joint25_coef, "coef")
    
    rows[[i]] <- data.table(
      analysis_level="each_protein_mediator",
      beta_scale="standardized_beta",
      beta_interpretation="path a: SD protein mediator per 1 SD pain score; Cox paths: logHR per 1 SD predictor",
      protein=protein,
      mediator=mediator,
      mediator_standardized_var=mediator_std,
      exposure=exposure_z,
      n=nrow(dat),
      events=sum(dat$status),
      n_mediators_total=length(mediator_cols),
      
      path_a_standardized_beta=a_beta,
      path_a_standardized_se=get_val(a_coef, "Std. Error"),
      path_a_standardized_t=get_val(a_coef, "t value"),
      path_a_standardized_p_value=get_val(a_coef, "Pr(>|t|)"),
      path_a_model="lm(protein_mediator_z ~ pain_score_z + covariates)",
      
      path_b_standardized_beta_logHR_single_protein=b_single_beta,
      path_b_standardized_HR_per_1SD_single_protein=exp(b_single_beta),
      path_b_standardized_se_single_protein=get_val(b_single_coef, "se(coef)"),
      path_b_standardized_z_single_protein=get_val(b_single_coef, "z"),
      path_b_standardized_p_value_single_protein=get_val(b_single_coef, "Pr(>|z|)"),
      path_b_single_model="coxph(Surv(time,status) ~ pain_score_z + protein_mediator_z + covariates)",
      
      path_c_standardized_beta_logHR=c_beta,
      path_c_standardized_HR_per_1SD=exp(c_beta),
      path_c_standardized_se=get_val(c_coef, "se(coef)"),
      path_c_standardized_z=get_val(c_coef, "z"),
      path_c_standardized_p_value=get_val(c_coef, "Pr(>|z|)"),
      path_c_model="coxph(Surv(time,status) ~ pain_score_z + covariates)",
      
      path_cprime_standardized_beta_logHR_single_protein=cprime_single_beta,
      path_cprime_standardized_HR_per_1SD_single_protein=exp(cprime_single_beta),
      path_cprime_standardized_se_single_protein=get_val(cprime_single_coef, "se(coef)"),
      path_cprime_standardized_z_single_protein=get_val(cprime_single_coef, "z"),
      path_cprime_standardized_p_value_single_protein=get_val(cprime_single_coef, "Pr(>|z|)"),
      
      product_ab_standardized_beta_logHR_single_protein=a_beta * b_single_beta,
      difference_c_cprime_standardized_beta_logHR_single_protein=c_beta - cprime_single_beta,
      
      path_b_standardized_beta_logHR_joint25_adjusted=b_joint25_beta,
      path_b_standardized_HR_per_1SD_joint25_adjusted=exp(b_joint25_beta),
      path_b_standardized_se_joint25_adjusted=get_val(b_joint25_coef, "se(coef)"),
      path_b_standardized_z_joint25_adjusted=get_val(b_joint25_coef, "z"),
      path_b_standardized_p_value_joint25_adjusted=get_val(b_joint25_coef, "Pr(>|z|)"),
      path_b_joint25_adjusted_model="coxph(Surv(time,status) ~ pain_score_z + all_25_protein_mediator_z + covariates)",
      
      path_cprime_standardized_beta_logHR_joint25_adjusted=cprime_joint25_beta,
      path_cprime_standardized_HR_per_1SD_joint25_adjusted=exp(cprime_joint25_beta),
      path_cprime_standardized_se_joint25_adjusted=get_val(cprime_joint25_coef, "se(coef)"),
      path_cprime_standardized_z_joint25_adjusted=get_val(cprime_joint25_coef, "z"),
      path_cprime_standardized_p_value_joint25_adjusted=get_val(cprime_joint25_coef, "Pr(>|z|)"),
      
      product_ab_standardized_beta_logHR_joint25_adjusted=a_beta * b_joint25_beta,
      difference_c_cprime_standardized_beta_logHR_joint25_adjusted=c_beta - cprime_joint25_beta
    )
  }
  
  rbindlist(rows, fill=TRUE)
}

compute_joint25_path_abc_beta <- function(dat, covariates, mediator_cols, iorw_result=NULL, exposure="pain_score") {
  # One-row joint 25-mediator path summary.
  # The 25 proteins have already been direction-aligned and z-scored in build_analysis_data().
  # We define a scalar joint mediator score as the mean of these 25 aligned z-scores.
  # This gives one interpretable joint path a and path b beta:
  #   path a: pain_score -> joint25_mediator_score
  #   path b: joint25_mediator_score -> PD, adjusting pain_score + covariates
  # IORW joint mediation still provides the preferred joint direct and indirect effects:
  #   path c      = total logHR from Cox model
  #   path c'     = IORW natural direct logHR from weighted Cox model
  #   indirect    = c - c'
  
  dat <- as.data.frame(dat)
  dat <- droplevels(dat)
  
  mediator_matrix <- as.matrix(dat[, mediator_cols, drop=FALSE])
  storage.mode(mediator_matrix) <- "double"
  dat$joint25_mediator_score <- rowMeans(mediator_matrix, na.rm=FALSE)
  
  needed_cols <- c("time", "status", exposure, "joint25_mediator_score", covariates)
  needed_cols <- unique(needed_cols[needed_cols %in% names(dat)])
  dat <- dat[complete.cases(dat[, needed_cols, drop=FALSE]),]
  dat <- droplevels(dat)
  
  covariates <- usable_covariates(dat, covariates)
  make_rhs <- function(vars) paste(vars, collapse="+")
  
  # path a: exposure -> joint mediator score
  a_formula <- as.formula(
    paste("joint25_mediator_score ~", make_rhs(c(exposure, covariates)))
  )
  a_fit <- lm(a_formula, data=dat)
  a_coef <- summary(a_fit)$coefficients[exposure,]
  
  # path c: exposure -> outcome, total effect
  c_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure, covariates)))
  )
  c_fit <- coxph(c_formula, data=dat)
  c_coef <- summary(c_fit)$coefficients[exposure,]
  
  # Traditional joint-score path b and c': outcome ~ exposure + joint mediator score + covariates
  bc_formula <- as.formula(
    paste("Surv(time, status) ~", make_rhs(c(exposure, "joint25_mediator_score", covariates)))
  )
  bc_fit <- coxph(bc_formula, data=dat)
  bc_coef <- summary(bc_fit)$coefficients
  b_coef <- bc_coef["joint25_mediator_score",]
  cprime_traditional_coef <- bc_coef[exposure,]
  
  get_iorw <- function(name, default=NA_real_) {
    if (is.null(iorw_result)) return(default)
    if (name %in% names(iorw_result)) return(iorw_result[[name]][1])
    default
  }
  
  out <- data.table(
    analysis_level="joint_25_mediators",
    mediator_set="25_direction_aligned_overlap_proteins",
    joint_mediator="joint25_mediator_score",
    joint_mediator_definition="rowMeans of the 25 direction-aligned z-scored protein mediators",
    exposure=exposure,
    n=nrow(dat),
    events=sum(dat$status),
    n_mediators=length(mediator_cols),
    
    path_a_beta=a_coef["Estimate"],
    path_a_se=a_coef["Std. Error"],
    path_a_t=a_coef["t value"],
    path_a_p_value=a_coef["Pr(>|t|)"],
    path_a_model="lm(joint25_mediator_score ~ pain_score + covariates)",
    
    path_b_beta_logHR=b_coef["coef"],
    path_b_HR=exp(b_coef["coef"]),
    path_b_se=b_coef["se(coef)"],
    path_b_z=b_coef["z"],
    path_b_p_value=b_coef["Pr(>|z|)"],
    path_b_model="coxph(Surv(time,status) ~ pain_score + joint25_mediator_score + covariates)",
    
    path_c_beta_logHR=c_coef["coef"],
    path_c_HR=exp(c_coef["coef"]),
    path_c_se=c_coef["se(coef)"],
    path_c_z=c_coef["z"],
    path_c_p_value=c_coef["Pr(>|z|)"],
    path_c_model="coxph(Surv(time,status) ~ pain_score + covariates)",
    
    path_cprime_beta_logHR_traditional_joint_score=cprime_traditional_coef["coef"],
    path_cprime_HR_traditional_joint_score=exp(cprime_traditional_coef["coef"]),
    path_cprime_se_traditional_joint_score=cprime_traditional_coef["se(coef)"],
    path_cprime_z_traditional_joint_score=cprime_traditional_coef["z"],
    path_cprime_p_value_traditional_joint_score=cprime_traditional_coef["Pr(>|z|)"],
    
    product_ab_beta_logHR=a_coef["Estimate"] * b_coef["coef"],
    difference_c_cprime_beta_logHR_traditional_joint_score=c_coef["coef"] - cprime_traditional_coef["coef"],
    
    path_c_beta_logHR_iorw=get_iorw("total_logHR"),
    path_c_HR_iorw=get_iorw("total_HR"),
    path_c_lower_iorw=get_iorw("total_lower"),
    path_c_upper_iorw=get_iorw("total_upper"),
    path_c_p_boot_iorw=get_iorw("total_p_boot"),
    
    path_cprime_beta_logHR_iorw=get_iorw("direct_logHR"),
    path_cprime_HR_iorw=get_iorw("direct_HR"),
    path_cprime_lower_iorw=get_iorw("direct_lower"),
    path_cprime_upper_iorw=get_iorw("direct_upper"),
    path_cprime_p_boot_iorw=get_iorw("direct_p_boot"),
    
    joint_indirect_beta_logHR_iorw=get_iorw("indirect_logHR"),
    joint_indirect_HR_iorw=get_iorw("indirect_HR"),
    joint_indirect_lower_iorw=get_iorw("indirect_lower"),
    joint_indirect_upper_iorw=get_iorw("indirect_upper"),
    joint_indirect_p_boot_iorw=get_iorw("indirect_p_boot"),
    proportion_mediated_iorw=get_iorw("proportion_mediated"),
    n_boot_success=get_iorw("n_boot_success")
  )
  
  out
}

bootstrap_iorw <- function(dat, covariates, mediator_cols, n_boot, weight_trim_quantile) {
  if (n_boot <= 0) return(data.table())
  boot_rows <- vector("list", n_boot)
  n <- nrow(dat)
  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size=n, replace=TRUE)
    boot_dat <- dat[idx,]
    row <- try(compute_iorw(boot_dat, covariates, mediator_cols, weight_trim_quantile)$result, silent=TRUE)
    if (inherits(row, "try-error")) {
      boot_rows[[b]] <- data.table(boot=b, error=as.character(row))
    } else {
      row[, boot := b]
      boot_rows[[b]] <- row
    }
    if (b %% 50 == 0) {
      cat("bootstrap", b, "of", n_boot, "\n")
      flush.console()
    }
  }
  rbindlist(boot_rows, fill=TRUE)
}

add_boot_ci <- function(result, boot_dt) {
  if (nrow(boot_dt)==0 || !("indirect_logHR" %in% names(boot_dt))) {
    result[, `:=`(
      indirect_lower=NA_real_, indirect_upper=NA_real_, indirect_p_boot=NA_real_,
      direct_lower=NA_real_, direct_upper=NA_real_, direct_p_boot=NA_real_,
      total_lower=NA_real_, total_upper=NA_real_, total_p_boot=NA_real_,
      proportion_lower=NA_real_, proportion_upper=NA_real_, proportion_p_boot=NA_real_,
      n_boot_success=0L
    )]
    return(result)
  }
  valid <- boot_dt[!is.na(indirect_logHR)]
  p_two <- function(x, obs) {
    if (length(x)==0 || is.na(obs)) return(NA_real_)
    2 * min(mean(x <= 0, na.rm=TRUE), mean(x >= 0, na.rm=TRUE))
  }
  result[, `:=`(
    indirect_lower=quantile(valid$indirect_logHR, 0.025, na.rm=TRUE),
    indirect_upper=quantile(valid$indirect_logHR, 0.975, na.rm=TRUE),
    indirect_p_boot=p_two(valid$indirect_logHR, indirect_logHR),
    direct_lower=quantile(valid$direct_logHR, 0.025, na.rm=TRUE),
    direct_upper=quantile(valid$direct_logHR, 0.975, na.rm=TRUE),
    direct_p_boot=p_two(valid$direct_logHR, direct_logHR),
    total_lower=quantile(valid$total_logHR, 0.025, na.rm=TRUE),
    total_upper=quantile(valid$total_logHR, 0.975, na.rm=TRUE),
    total_p_boot=p_two(valid$total_logHR, total_logHR),
    proportion_lower=quantile(valid$proportion_mediated, 0.025, na.rm=TRUE),
    proportion_upper=quantile(valid$proportion_mediated, 0.975, na.rm=TRUE),
    proportion_p_boot=p_two(valid$proportion_mediated, proportion_mediated),
    indirect_HR_lower=exp(quantile(valid$indirect_logHR, 0.025, na.rm=TRUE)),
    indirect_HR_upper=exp(quantile(valid$indirect_logHR, 0.975, na.rm=TRUE)),
    direct_HR_lower=exp(quantile(valid$direct_logHR, 0.025, na.rm=TRUE)),
    direct_HR_upper=exp(quantile(valid$direct_logHR, 0.975, na.rm=TRUE)),
    total_HR_lower=exp(quantile(valid$total_logHR, 0.025, na.rm=TRUE)),
    total_HR_upper=exp(quantile(valid$total_logHR, 0.975, na.rm=TRUE)),
    n_boot_success=nrow(valid)
  )]
  result
}

covariates <- c("age","gender","ethnicity","deprivation","income","university_education",
                "smoking","alcohol","sedentary","obesity","prs",paste("comp",1:20,sep="_"),
                "protein_cov_zscore","season_final","Batch","protein_time")

analysis <- build_analysis_data()
model_data <- analysis$data
mediator_cols <- analysis$mediator_cols

model_cols <- c("eid","time","status","pain_score",covariates,mediator_cols)
iorw_data <- as.data.frame(model_data[, model_cols, with=FALSE])
iorw_data <- iorw_data[complete.cases(iorw_data),]
iorw_data <- droplevels(iorw_data)

# Add a single joint mediator score for the 25 aligned proteins.
# This is used only for the one-row path a / path b beta summary.
iorw_data$joint25_mediator_score <- rowMeans(as.matrix(iorw_data[, mediator_cols, drop=FALSE]), na.rm=FALSE)

# Standardized variables for beta values comparable to paper-style path figures.
# pain_score_z: 1 unit = 1 SD higher pain score.
# joint25_mediator_score_z: 1 unit = 1 SD higher joint 25-protein mediator score.
iorw_data$pain_score_z <- zscore_numeric(iorw_data$pain_score)
iorw_data$joint25_mediator_score_z <- zscore_numeric(iorw_data$joint25_mediator_score)

fwrite(as.data.table(iorw_data), out_data_path)

cat("running joint IORW mediation...\n")
iorw <- compute_iorw(iorw_data, covariates, mediator_cols, weight_trim_quantile)
result <- iorw$result

cat("running standardized joint IORW mediation beta models...\n")
iorw_std <- compute_iorw_standardized(
  dat=iorw_data,
  covariates=covariates,
  mediator_cols=mediator_cols,
  weight_trim_quantile=weight_trim_quantile,
  exposure="pain_score_z"
)
std_result <- iorw_std$result
fwrite(std_result, out_std_result_path)

cat("running bootstrap with", n_boot, "resamples...\n")
boot_dt <- bootstrap_iorw(iorw_data, covariates, mediator_cols, n_boot, weight_trim_quantile)
if (nrow(boot_dt)>0) fwrite(boot_dt, out_boot_path)
result <- add_boot_ci(result, boot_dt)
fwrite(result, out_result_path)

cat("computing joint 25 path a, b, c beta results...\n")
joint25_path_abc_result <- compute_joint25_path_abc_beta(
  dat=iorw_data,
  covariates=covariates,
  mediator_cols=mediator_cols,
  iorw_result=result,
  exposure="pain_score"
)
fwrite(joint25_path_abc_result, out_joint25_path_abc_path)

cat("computing standardized joint 25 path a, b, c beta results...\n")
joint25_std_path_abc_result <- compute_joint25_standardized_path_abc_beta(
  dat=iorw_data,
  covariates=covariates,
  mediator_cols=mediator_cols,
  iorw_std_result=std_result,
  exposure_z="pain_score_z",
  joint_mediator_z="joint25_mediator_score_z"
)
fwrite(joint25_std_path_abc_result, out_joint25_std_path_abc_path)

cat("computing each-protein standardized path a, b, c beta results...\n")
each_protein_std_path_abc_result <- compute_each_protein_standardized_path_abc_beta(
  dat=iorw_data,
  covariates=covariates,
  mediator_cols=mediator_cols,
  proteins=analysis$proteins,
  exposure_z="pain_score_z"
)
fwrite(each_protein_std_path_abc_result, out_each_protein_std_path_abc_path)

write_json_simple(list(
  method="Joint IORW mediation for continuous pain_score exposure and 25 direction-aligned overlap protein mediators. Stabilized density-ratio weights are estimated as f(pain_score|covariates)/f(pain_score|covariates,mediators) using Gaussian linear models; weighted Cox model estimates natural direct effect; indirect effect is total minus direct on log-HR scale.",
  n=nrow(iorw_data),
  events=sum(iorw_data$status),
  n_mediators=length(mediator_cols),
  mediators=mediator_cols,
  weight_trim_quantile=weight_trim_quantile,
  n_boot_requested=n_boot,
  n_boot_success=ifelse(nrow(boot_dt)>0, sum(!is.na(boot_dt$indirect_logHR)), 0),
  data_output=out_data_path,
  result_output=out_result_path,
  standardized_result_output=out_std_result_path,
  joint25_path_abc_output=out_joint25_path_abc_path,
  joint25_standardized_path_abc_output=out_joint25_std_path_abc_path,
  each_protein_standardized_path_abc_output=out_each_protein_std_path_abc_path,
  bootstrap_output=ifelse(nrow(boot_dt)>0, out_boot_path, NA)
), out_summary_path)

cat("outputs:\n")
cat(out_data_path, "\n")
cat(out_result_path, "\n")
cat(out_std_result_path, "\n")
cat(out_joint25_path_abc_path, "\n")
cat(out_joint25_std_path_abc_path, "\n")
cat(out_each_protein_std_path_abc_path, "\n")
if (nrow(boot_dt)>0) cat(out_boot_path, "\n")
cat(out_summary_path, "\n")
print(result)
print(std_result)
print(joint25_path_abc_result)
print(joint25_std_path_abc_result)
print(each_protein_std_path_abc_result)
