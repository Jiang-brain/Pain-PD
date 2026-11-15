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
score_path <- "pain_related_proteomic_score_stable_nested_pfs.csv"
coef_path <- "pain_related_proteomic_score_stable_nested_coefficients.csv"
protein_cov_path <- find_input("PD_data_matched_protein_covariates.csv", "../PA/PD_data_matched_protein_covariates.csv", "../PA2/PD_data_matched_protein_covariates.csv")
protein_time_path <- find_input("PD_data_matched_protein_times.csv", "../PA/PD_data_matched_protein_times.csv", "../PA2/PD_data_matched_protein_times.csv")

out_path <- "pd_incidence_stable_nested_pfs_cindex.csv"
out_summary_path <- "pd_incidence_stable_nested_pfs_summary.json"

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

model_row <- function(method, model_type, fit, reference_fit=NULL, n_extra_terms=NA_integer_, extra_terms="") {
  c_aug <- unname(summary(fit)$concordance["C"])
  se_aug <- unname(summary(fit)$concordance["se(C)"])
  c_base <- if (is.null(reference_fit)) NA_real_ else unname(summary(reference_fit)$concordance["C"])
  se_base <- if (is.null(reference_fit)) NA_real_ else unname(summary(reference_fit)$concordance["se(C)"])
  data.frame(
    method=method,
    model_type=model_type,
    n=fit$n,
    events=fit$nevent,
    n_extra_terms=n_extra_terms,
    extra_terms=extra_terms,
    c_index=c_aug,
    c_index_se=se_aug,
    reference_c_index=c_base,
    reference_c_index_se=se_base,
    c_index_delta=ifelse(is.na(c_base), NA_real_, c_aug-c_base),
    stringsAsFactors=FALSE
  )
}

write_json_simple <- function(x, path) {
  scalar <- function(v) {
    if (length(v) == 1 && is.na(v)) return("null")
    if (is.character(v)) return(paste0('"', gsub('\\\\', '/', v), '"'))
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

cat("building PD Cox base data...\n")
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

outcome_data <- data.frame(eid=data_censoring$id, status=data_censoring$label, time=as.numeric(data_censoring$censoring))
analysis_data <- merge(outcome_data, covars, by="eid")
analysis_data <- analysis_data[complete.cases(analysis_data),]
analysis_data <- analysis_data[analysis_data$age>=40,]
analysis_data <- analysis_data[!(analysis_data$time<365*2 & analysis_data$status==1),]
analysis_data <- merge(as.data.table(analysis_data), protein_cov, by="eid")
analysis_data <- merge(analysis_data, protein_times, by="eid")

covariates <- c("age","gender","ethnicity","deprivation","income","university_education",
                "smoking","alcohol","sedentary","obesity","prs",paste("comp",1:20,sep="_"),
                "protein_cov_zscore","season_final","Batch","protein_time")

score_dt <- fread(score_path)
coef_dt <- fread(coef_path)
stable_dt <- coef_dt[stable_all_10_folds==TRUE]
methods <- sort(unique(stable_dt$method))
stable_proteins <- sort(unique(stable_dt$protein))

cat("reading stable proteins:", length(stable_proteins), "\n")
protein_raw <- fread(protein_path, select=c("eid","ins_index",stable_proteins))
setorder(protein_raw, eid, ins_index)
protein_raw <- protein_raw[, .SD[1], by=eid]
protein_raw[, ins_index := NULL]

results <- list()
for (method_name in methods) {
  cat("Cox models for", method_name, "\n")
  method_score <- paste0("stable_", method_name, "_pfs_z")
  method_coef <- stable_dt[method == method_name]
  method_proteins <- method_coef$protein
  protein_term_names <- make.unique(make.names(paste0("prot_", method_name, "_", method_proteins)))
  protein_dt <- protein_raw[, c("eid", method_proteins), with=FALSE]
  for (i in seq_along(method_proteins)) {
    protein <- method_proteins[i]
    values <- protein_dt[[protein]]
    impute_value <- method_coef[method_coef$protein == protein, mean_impute_value][1]
    values[is.na(values)] <- impute_value
    protein_dt[[protein_term_names[i]]] <- as.numeric(scale(values))
  }
  protein_dt <- protein_dt[, c("eid", protein_term_names), with=FALSE]

  dat <- merge(copy(analysis_data), score_dt[, c("eid", method_score), with=FALSE], by="eid")
  dat <- merge(dat, protein_dt, by="eid")
  model_cols <- c("time","status",covariates,method_score,protein_term_names)
  model_data <- as.data.frame(dat[, model_cols, with=FALSE])
  model_data <- model_data[complete.cases(model_data),]
  model_data <- droplevels(model_data)
  method_covariates <- usable_covariates(model_data, covariates)
  covariate_terms <- paste(method_covariates, collapse="+")
  protein_terms <- paste(protein_term_names, collapse="+")

  base_fit <- coxph(as.formula(paste("Surv(time, status) ~", covariate_terms)), data=model_data)
  pfs_alone <- coxph(as.formula(paste("Surv(time, status) ~", method_score)), data=model_data)
  pfs_plus <- coxph(as.formula(paste("Surv(time, status) ~", covariate_terms, "+", method_score)), data=model_data)
  protein_alone <- coxph(as.formula(paste("Surv(time, status) ~", protein_terms)), data=model_data, control=coxph.control(iter.max=30))
  protein_plus <- coxph(as.formula(paste("Surv(time, status) ~", covariate_terms, "+", protein_terms)), data=model_data, control=coxph.control(iter.max=30))

  pfs_coef <- summary(pfs_plus)$coefficients[method_score,]
  pfs_conf <- summary(pfs_plus)$conf.int[method_score,]
  rows <- rbind(
    model_row(method_name, "traditional_covariates", base_fit, NULL, 0, "traditional_covariates"),
    model_row(method_name, "stable_pfs_alone", pfs_alone, NULL, 1, method_score),
    model_row(method_name, "traditional_plus_stable_pfs", pfs_plus, base_fit, 1, method_score),
    model_row(method_name, "stable_proteins_alone", protein_alone, NULL, length(protein_term_names), paste0(length(protein_term_names), "_stable_proteins")),
    model_row(method_name, "traditional_plus_stable_proteins", protein_plus, base_fit, length(protein_term_names), paste0(length(protein_term_names), "_stable_proteins"))
  )
  rows$stable_pfs_beta <- NA_real_
  rows$stable_pfs_HR <- NA_real_
  rows$stable_pfs_lower <- NA_real_
  rows$stable_pfs_upper <- NA_real_
  rows$stable_pfs_p_value <- NA_real_
  rows[rows$model_type=="traditional_plus_stable_pfs", c("stable_pfs_beta","stable_pfs_HR","stable_pfs_lower","stable_pfs_upper","stable_pfs_p_value")] <-
    list(pfs_coef["coef"], pfs_conf["exp(coef)"], pfs_conf["lower .95"], pfs_conf["upper .95"], pfs_coef["Pr(>|z|)"])
  results[[length(results)+1]] <- rows
}

result_dt <- rbindlist(results, fill=TRUE)
fwrite(result_dt, out_path)
write_json_simple(list(
  method="Cox c-index for stable nested-CV LASSO/ElasticNet PFS and stable proteins, excluding PD events within 2 years.",
  score_path=score_path,
  coefficient_path=coef_path,
  output=out_path,
  methods=methods,
  n_unique_stable_proteins=length(stable_proteins)
), out_summary_path)

print(result_dt)
cat("output:", out_path, "\n")
