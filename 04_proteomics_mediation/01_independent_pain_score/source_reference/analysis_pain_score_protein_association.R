data_path <- "PD_data.csv"
score_path <- "composite_pain_scores.csv"
protein_path <- "protein_new_missing60_filtered.csv"
protein_cov_path <- "PD_data_matched_protein_covariates.csv"
protein_time_path <- "PD_data_matched_protein_times.csv"
out_path <- "pain_score_protein_association_results.csv"

.libPaths(c(normalizePath("R_libs"), .libPaths()))
library(data.table)

expand_multi_select <- function(x) {
  raw_values <- as.character(x)
  raw_values[is.na(raw_values)] <- ""
  split_values <- strsplit(raw_values, "\\|")
  max_len <- max(lengths(split_values))
  expanded <- matrix(NA_real_, nrow=length(raw_values), ncol=max_len)
  for (i in seq_along(split_values)) {
    values <- suppressWarnings(as.numeric(split_values[[i]]))
    if (length(values)>0) {
      expanded[i,seq_along(values)] <- values
    }
  }
  as.data.frame(expanded)
}

get_multi_select <- function(dat, field_id) {
  cols <- grep(paste0("^X",field_id,"\\.0\\."), names(dat), value=TRUE)
  if (length(cols)==1) {
    return(expand_multi_select(dat[[cols]]))
  }
  as.data.frame(dat[,cols,drop=FALSE])
}

covariate_cols <- c(
  "eid","21003-0.0","31-0.0","21000-0.0","21000-1.0","21000-2.0",
  "22189-0.0","738-0.0","6138-0.0","20116-0.0","20117-0.0",
  "1070-0.0","21001-0.0","26260-0.0",paste0("22009-0.",1:20)
)

datas <- data.frame(fread(data_path, select=covariate_cols))

age <- datas$X21003.0.0
gender <- datas$X31.0.0

ethnicity <- cbind(datas$X21000.0.0,datas$X21000.1.0,datas$X21000.2.0)
index_nan <- which(ethnicity[,1]==-1|ethnicity[,1]==-3)
ethnicity[index_nan,1] <- pmax(ethnicity[index_nan,1],ethnicity[index_nan,2],ethnicity[index_nan,3],na.rm=TRUE)
race <- ethnicity[,1]
index_white <- which((race %in% c(1,1001,1002,1003))==TRUE)
index_nonwhite <- which((race %in% c(1,1001,1002,1003))==FALSE)
index_not_answer <- which(race==-1|race==-3)
race[index_white] <- 1
race[index_nonwhite] <- 0
race[index_not_answer] <- NA
ethnicity <- race

deprivation <- datas$X22189.0.0
tertile_point <- quantile(deprivation, probs=seq(0,1,1/3), na.rm=TRUE)
index_low <- which(deprivation<tertile_point[2])
index_middle <- which(deprivation>=tertile_point[2]&deprivation<tertile_point[3])
index_high <- which(deprivation>=tertile_point[3])
deprivation[index_low] <- 0
deprivation[index_middle] <- 1
deprivation[index_high] <- 2

income <- datas$X738.0.0
index_not_answer <- which(income==-3)
index_non <- which(income==-1)
index_low <- which(income<=3&income>=1)
index_medium <- which(income==4)
index_high <- which(income==5)
income[index_non] <- 0
income[index_low] <- 1
income[index_medium] <- 2
income[index_high] <- 3
income[index_not_answer] <- 0

education <- get_multi_select(datas,"6138")
for (i in seq_len(ncol(education))) {
  index_no_answer <- which(education[,i]==-3)
  education[index_no_answer,i] <- NA
  index_none <- which(education[,i]==-7)
  education[index_none,i] <- 7

  index_20 <- which(education[,i]==1)
  index_13 <- which(education[,i]==2)
  index_10 <- which(education[,i]==3|education[,i]==4)
  index_19 <- which(education[,i]==5)
  index_15 <- which(education[,i]==6)
  index_7 <- which(education[,i]==7)
  education[index_20,i] <- 20
  education[index_13,i] <- 13
  education[index_10,i] <- 10
  education[index_19,i] <- 19
  education[index_15,i] <- 15
  education[index_7,i] <- 7
}
university_education <- do.call(pmax,c(education,list(na.rm=TRUE)))
index_no_education <- which(university_education==7)
index_no_university <- which(university_education<20&university_education!=7)
index_university <- which(university_education==20)
university_education[index_no_education] <- 0
university_education[index_no_university] <- 1
university_education[index_university] <- 2

smoking <- datas$X20116.0.0
index_no_answer <- which(smoking==-3)
smoking[index_no_answer] <- NA
index_current <- which(smoking==2)
smoking[index_current] <- 1

alcohol <- datas$X20117.0.0
index_no_answer <- which(alcohol<0)
alcohol[index_no_answer] <- NA
index_current <- which(alcohol==2)
alcohol[index_current] <- 1

sedentary <- datas$X1070.0.0
index_na <- which(sedentary==-1|sedentary==-3)
index_low <- which((sedentary<4&sedentary>=0)|sedentary==-10)
index_high <- which(sedentary>=4)
sedentary[index_na] <- NA
sedentary[index_low] <- 0
sedentary[index_high] <- 1

bmi <- datas$X21001.0.0
obesity <- bmi
index_underweight <- which(bmi<18.5)
index_normal <- which(bmi>=18.5&bmi<25)
index_overweight <- which(bmi>=25&bmi<30)
index_obese1 <- which(bmi>=30&bmi<35)
index_obese2 <- which(bmi>=35)
obesity[index_underweight] <- NA
obesity[index_normal] <- 0
obesity[index_overweight] <- 1
obesity[index_obese1] <- 2
obesity[index_obese2] <- 3

pd_prs <- datas$X26260.0.0
tertile_point <- quantile(pd_prs, probs=seq(0,1,1/3), na.rm=TRUE)
index_low <- which(pd_prs<tertile_point[2])
index_middle <- which(pd_prs>=tertile_point[2]&pd_prs<tertile_point[3])
index_high <- which(pd_prs>=tertile_point[3])
pd_prs[index_low] <- 0
pd_prs[index_middle] <- 1
pd_prs[index_high] <- 2

covars <- data.frame(
  eid=datas$eid, age=age, gender=factor(gender), ethnicity=factor(ethnicity),
  deprivation=factor(deprivation), income=factor(income),
  university_education=factor(university_education), smoking=factor(smoking),
  alcohol=factor(alcohol), sedentary=factor(sedentary), obesity=factor(obesity),
  prs=factor(pd_prs)
)
for (k in 1:20) {
  covars[[paste0("comp_",k)]] <- datas[[paste0("X22009.0.",k)]]
}
rm(datas)

scores <- fread(score_path, select=c("eid","pain_score"))
protein_cov <- fread(protein_cov_path)
setnames(protein_cov, "X74.0.0_zscore", "protein_cov_zscore")
protein_cov[, season_final := factor(season_final)]
protein_cov[, Batch := factor(Batch)]
protein_times <- fread(protein_time_path)
time_protein_cols <- setdiff(names(protein_times), "eid")
setnames(protein_times, time_protein_cols, paste0(time_protein_cols, "_protein_time"))
proteins <- fread(protein_path)
protein_names <- setdiff(names(proteins), c("eid","ins_index"))

analysis_data <- merge(scores, covars, by="eid")
analysis_data <- merge(analysis_data, protein_cov, by="eid")
analysis_data <- merge(analysis_data, protein_times, by="eid")
analysis_data <- merge(analysis_data, proteins, by="eid")
rm(scores,covars,protein_cov,protein_times,proteins)

covariates <- c("age","gender","ethnicity","deprivation","income","university_education",
                "smoking","alcohol","sedentary","obesity","prs",paste("comp",1:20,sep="_"),
                "protein_cov_zscore","season_final","Batch")
covariate_terms <- paste(c(covariates, "protein_time"), collapse="+")

base_formula <- as.formula(paste("pain_score ~", paste(covariates, collapse="+")))
base_model <- lm(base_formula, data=analysis_data, na.action=na.exclude)
base_n <- nobs(base_model)

results <- vector("list", length(protein_names))
for (j in seq_along(protein_names)) {
  protein <- protein_names[j]
  time_col <- paste0(protein, "_protein_time")
  model_cols <- c("pain_score", covariates, protein, time_col)
  model_data <- as.data.frame(analysis_data[, model_cols, with=FALSE])
  names(model_data)[(ncol(model_data)-1):ncol(model_data)] <- c("protein_value","protein_time")
  fit <- try(lm(
    as.formula(paste("pain_score ~ protein_value +", covariate_terms)),
    data=model_data
  ), silent=TRUE)

  if (inherits(fit,"try-error") || !("protein_value" %in% rownames(summary(fit)$coefficients))) {
    results[[j]] <- data.frame(protein=protein,n=NA,beta=NA,se=NA,t=NA,p_value=NA)
  } else {
    coef_row <- summary(fit)$coefficients["protein_value",]
    results[[j]] <- data.frame(
      protein=protein,
      n=nobs(fit),
      beta=coef_row["Estimate"],
      se=coef_row["Std. Error"],
      t=coef_row["t value"],
      p_value=coef_row["Pr(>|t|)"]
    )
  }

  if (j %% 250 == 0) {
    cat("finished", j, "of", length(protein_names), "\n")
  }
}

result <- rbindlist(results)
result[, fdr := p.adjust(p_value, method="BH")]
result[, bonferroni := p.adjust(p_value, method="bonferroni")]
setorder(result, p_value)
fwrite(result, out_path)

cat("base complete-case n:", base_n, "\n")
cat("proteins tested:", nrow(result), "\n")
cat("FDR<0.05:", sum(result$fdr<0.05, na.rm=TRUE), "\n")
cat("Bonferroni<0.05:", sum(result$bonferroni<0.05, na.rm=TRUE), "\n")
cat("output:", out_path, "\n")
