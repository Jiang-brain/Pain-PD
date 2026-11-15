find_input <- function(...) {
  candidates <- c(...)
  existing <- candidates[file.exists(candidates)]
  if (length(existing)==0) stop("None of these input files exist: ", paste(candidates, collapse=", "))
  existing[1]
}

data_path <- find_input("PD_data.csv", "../PA/PD_data.csv", "../PA2/PD_data.csv")
pca_score_path <- "overlap25_directional_pca_signature.csv"
intersection_path <- "protein_intersection_FDR05.csv"
protein_cov_path <- find_input("PD_data_matched_protein_covariates.csv", "../PA/PD_data_matched_protein_covariates.csv", "../PA2/PD_data_matched_protein_covariates.csv")
protein_time_path <- find_input("PD_data_matched_protein_times.csv", "../PA/PD_data_matched_protein_times.csv", "../PA2/PD_data_matched_protein_times.csv")
out_path <- "overlap25_pca_prediction_cindex_5y_10y_full_followup.csv"
nri_out_path <- "overlap25_pca_prediction_nri_5y_10y_full_followup.csv"
followup_2y <- 365.25 * 2

.libPaths(c(normalizePath("R_libs"), .libPaths()))
library(data.table)
library(survival)

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

fit_row <- function(task, model_name, fit, ref_fit=NULL, n_extra_terms=0L, extra_terms="") {
  cidx <- unname(summary(fit)$concordance[1])
  se <- unname(summary(fit)$concordance["se(C)"])
  ref <- if (is.null(ref_fit)) NA_real_ else unname(summary(ref_fit)$concordance[1])
  ref_se <- if (is.null(ref_fit)) NA_real_ else unname(summary(ref_fit)$concordance["se(C)"])
  data.frame(
    task=task,
    model=model_name,
    n=fit$n,
    events=fit$nevent,
    n_extra_terms=n_extra_terms,
    extra_terms=extra_terms,
    c_index=cidx,
    c_index_se=se,
    reference_c_index=ref,
    reference_c_index_se=ref_se,
    c_index_delta=ifelse(is.na(ref), NA_real_, cidx-ref),
    stringsAsFactors=FALSE
  )
}

predict_cox_absolute_risk <- function(fit, model_data, eval_time) {
  bh <- basehaz(fit, centered=FALSE)
  bh <- bh[is.finite(bh$time) & is.finite(bh$hazard),]
  bh <- bh[bh$time <= eval_time,]
  if (nrow(bh) == 0) {
    h0 <- 0
  } else {
    h0 <- bh$hazard[nrow(bh)]
  }
  lp <- predict(fit, newdata=model_data, type="lp", reference="zero")
  risk <- 1 - exp(-h0 * exp(lp))
  pmin(pmax(as.numeric(risk), 0), 1)
}

continuous_nri_row <- function(task, comparison, old_fit, new_fit, model_data, eval_time) {
  old_risk <- predict_cox_absolute_risk(old_fit, model_data, eval_time)
  new_risk <- predict_cox_absolute_risk(new_fit, model_data, eval_time)
  status <- model_data$status_horizon
  ok <- is.finite(old_risk) & is.finite(new_risk) & !is.na(status)
  old_risk <- old_risk[ok]
  new_risk <- new_risk[ok]
  status <- status[ok]
  delta <- new_risk - old_risk
  event_delta <- delta[status == 1]
  nonevent_delta <- delta[status == 0]

  prop_up <- function(x) mean(x > 0)
  prop_down <- function(x) mean(x < 0)
  event_up <- prop_up(event_delta)
  event_down <- prop_down(event_delta)
  nonevent_up <- prop_up(nonevent_delta)
  nonevent_down <- prop_down(nonevent_delta)
  nri_event <- event_up - event_down
  nri_nonevent <- nonevent_down - nonevent_up
  nri <- nri_event + nri_nonevent

  n_event <- length(event_delta)
  n_nonevent <- length(nonevent_delta)
  se_event <- sqrt((event_up + event_down - nri_event^2) / n_event)
  se_nonevent <- sqrt((nonevent_up + nonevent_down - nri_nonevent^2) / n_nonevent)
  nri_se <- sqrt(se_event^2 + se_nonevent^2)
  normal_p <- function(estimate, se) {
    if (!is.finite(estimate) || !is.finite(se) || se <= 0) return(NA_real_)
    2 * pnorm(abs(estimate / se), lower.tail=FALSE)
  }
  nri_p <- normal_p(nri, nri_se)
  event_p <- normal_p(nri_event, se_event)
  nonevent_p <- normal_p(nri_nonevent, se_nonevent)

  data.frame(
    task=task,
    comparison=comparison,
    n=length(status),
    events=n_event,
    nonevents=n_nonevent,
    eval_time_days=eval_time,
    old_mean_risk=mean(old_risk),
    new_mean_risk=mean(new_risk),
    nri=nri,
    nri_se=nri_se,
    nri_lower=nri - 1.96 * nri_se,
    nri_upper=nri + 1.96 * nri_se,
    nri_p=nri_p,
    nri_event=nri_event,
    nri_event_se=se_event,
    nri_event_lower=nri_event - 1.96 * se_event,
    nri_event_upper=nri_event + 1.96 * se_event,
    nri_event_p=event_p,
    event_up=event_up,
    event_down=event_down,
    nri_nonevent=nri_nonevent,
    nri_nonevent_se=se_nonevent,
    nri_nonevent_lower=nri_nonevent - 1.96 * se_nonevent,
    nri_nonevent_upper=nri_nonevent + 1.96 * se_nonevent,
    nri_nonevent_p=nonevent_p,
    nonevent_down=nonevent_down,
    nonevent_up=nonevent_up,
    method="continuous NRI for 25 direction-aligned proteins; SE from reclassification proportions",
    stringsAsFactors=FALSE
  )
}

cat("building data...\n")
same <- fread(intersection_path)
same[, same_direction := sign(beta_pain_score) == sign(beta_PD_cox)]
overlap_proteins <- same[same_direction==TRUE, protein]
if (length(overlap_proteins) != 25) stop("Expected 25 same-direction proteins.")
protein_cols <- paste0("aligned_", overlap_proteins, "_z")

data_censoring <- make_censoring_data(data_path)
covars <- make_covariates(data_path)
protein_cov <- fread(protein_cov_path)
setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- fread(protein_time_path)
time_cols <- setdiff(names(protein_times), "eid")
protein_times[, protein_time := rowMeans(.SD, na.rm=TRUE), .SDcols=time_cols]
protein_times <- protein_times[, .(eid, protein_time)]

outcome_data <- data.table(eid=data_censoring$id, status=data_censoring$label, time=as.numeric(data_censoring$censoring))
pca_scores <- fread(pca_score_path, select=c("eid","overlap25_pca_signature",protein_cols))
analysis_data <- merge(outcome_data, covars, by="eid")
analysis_data <- analysis_data[complete.cases(analysis_data),]
analysis_data <- analysis_data[analysis_data$age>=40,]
analysis_data <- analysis_data[!(analysis_data$time<=followup_2y & analysis_data$status==1),]
analysis_data <- merge(as.data.table(analysis_data), protein_cov, by="eid")
analysis_data <- merge(analysis_data, protein_times, by="eid")
analysis_data <- merge(analysis_data, pca_scores, by="eid")

covariates <- c("age","gender","ethnicity","deprivation","income","university_education",
                "smoking","alcohol","sedentary","obesity","prs",paste("comp",1:20,sep="_"),
                "protein_cov_zscore","season_final","Batch","protein_time")

make_task <- function(data, horizon_years) {
  horizon <- 365.25 * horizon_years
  task <- copy(data)
  task[, time_horizon := pmin(time, horizon)]
  task[, status_horizon := as.integer(status==1 & time <= horizon)]
  task[time >= horizon | status_horizon==1]
}

make_full_followup_task <- function(data) {
  task <- copy(data[time > followup_2y])
  task[, time_horizon := time - followup_2y]
  task[, status_horizon := as.integer(status==1 & time > followup_2y)]
  task
}

results <- list()
nri_results <- list()
task_specs <- list(
  list(name="incident_PD_within_5y", data=make_task(analysis_data, 5)),
  list(name="incident_PD_within_10y", data=make_task(analysis_data, 10)),
  list(name="incident_PD_full_followup_after_2y", data=make_full_followup_task(analysis_data))
)

for (task_spec in task_specs) {
  task_name <- task_spec$name
  task_data <- task_spec$data
  model_cols <- c("time_horizon","status_horizon",covariates,"overlap25_pca_signature",protein_cols)
  model_data <- as.data.frame(task_data[, model_cols, with=FALSE])
  model_data <- model_data[complete.cases(model_data),]
  model_data <- droplevels(model_data)
  covs <- usable_covariates(model_data, covariates)
  cov_terms <- paste(covs, collapse="+")
  protein_terms <- paste(protein_cols, collapse="+")
  base_fit <- coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", cov_terms)), data=model_data)
  pca_fit <- coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", cov_terms, "+ overlap25_pca_signature")), data=model_data)
  protein_alone_fit <- coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", protein_terms)), data=model_data)
  protein_plus_fit <- coxph(as.formula(paste("Surv(time_horizon, status_horizon) ~", cov_terms, "+", protein_terms)), data=model_data)
  pca_alone_fit <- coxph(Surv(time_horizon, status_horizon) ~ overlap25_pca_signature, data=model_data)
  results[[length(results)+1]] <- fit_row(task_name, "traditional_covariates", base_fit, NULL, 0, "traditional_covariates")
  results[[length(results)+1]] <- fit_row(task_name, "pca_signature_alone", pca_alone_fit, NULL, 1, "overlap25_pca_signature")
  results[[length(results)+1]] <- fit_row(task_name, "traditional_plus_pca_signature", pca_fit, base_fit, 1, "overlap25_pca_signature")
  results[[length(results)+1]] <- fit_row(task_name, "overlap25_proteins_alone", protein_alone_fit, NULL, length(protein_cols), "25_direction_aligned_overlap_proteins")
  results[[length(results)+1]] <- fit_row(task_name, "traditional_plus_overlap25_proteins", protein_plus_fit, base_fit, length(protein_cols), "25_direction_aligned_overlap_proteins")

  eval_time <- max(model_data$time_horizon, na.rm=TRUE)
  # NRI is restricted to the 25-protein model; PCA is retained only for c-index output.
  nri_results[[length(nri_results)+1]] <- continuous_nri_row(
    task_name,
    "traditional_covariates_vs_traditional_plus_overlap25_proteins",
    base_fit,
    protein_plus_fit,
    model_data,
    eval_time
  )
}

result <- rbindlist(results)
fwrite(result, out_path)
print(result)
cat("output:", out_path, "\n")

nri_result <- rbindlist(nri_results)
fwrite(nri_result, nri_out_path)
print(nri_result)
cat("NRI output:", nri_out_path, "\n")
