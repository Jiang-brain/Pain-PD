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

.libPaths(c(normalizePath("R_libs"), .libPaths()))
library(data.table)
library(survival)

set.seed(20260602)

n_boot <- as.integer(Sys.getenv("IORW_N_BOOT", "1000"))
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
fwrite(as.data.table(iorw_data), out_data_path)

cat("running joint IORW mediation...\n")
iorw <- compute_iorw(iorw_data, covariates, mediator_cols, weight_trim_quantile)
result <- iorw$result

cat("running bootstrap with", n_boot, "resamples...\n")
boot_dt <- bootstrap_iorw(iorw_data, covariates, mediator_cols, n_boot, weight_trim_quantile)
if (nrow(boot_dt)>0) fwrite(boot_dt, out_boot_path)
result <- add_boot_ci(result, boot_dt)
fwrite(result, out_result_path)

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
  bootstrap_output=ifelse(nrow(boot_dt)>0, out_boot_path, NA)
), out_summary_path)

cat("outputs:\n")
cat(out_data_path, "\n")
cat(out_result_path, "\n")
if (nrow(boot_dt)>0) cat(out_boot_path, "\n")
cat(out_summary_path, "\n")
print(result)
