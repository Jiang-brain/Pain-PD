base_script <- "analysis_stable_nested_pfs_pd_cindex_base.R"
if (!file.exists(base_script)) stop("Missing base Cox script: ", base_script)

# Run the original stable nested PFS c-index analysis first.
source(base_script, local=globalenv())

nri_out_path <- "pd_incidence_stable_nested_pfs_nri.csv"
nri_summary_path <- "pd_incidence_stable_nested_pfs_nri_summary.json"

predict_cox_absolute_risk <- function(fit, model_data, eval_time) {
  bh <- basehaz(fit, centered=FALSE)
  bh <- bh[is.finite(bh$time) & is.finite(bh$hazard),]
  bh <- bh[bh$time <= eval_time,]
  h0 <- if (nrow(bh)==0) 0 else bh$hazard[nrow(bh)]
  lp <- predict(fit, newdata=model_data, type="lp", reference="zero")
  risk <- 1 - exp(-h0 * exp(lp))
  pmin(pmax(as.numeric(risk), 0), 1)
}

continuous_nri_row <- function(method, comparison, old_fit, new_fit, model_data, eval_time) {
  old_risk <- predict_cox_absolute_risk(old_fit, model_data, eval_time)
  new_risk <- predict_cox_absolute_risk(new_fit, model_data, eval_time)
  status <- model_data$status
  ok <- is.finite(old_risk) & is.finite(new_risk) & !is.na(status)
  old_risk <- old_risk[ok]
  new_risk <- new_risk[ok]
  status <- status[ok]
  delta <- new_risk - old_risk
  event_delta <- delta[status==1]
  nonevent_delta <- delta[status==0]

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

  data.frame(
    method=method,
    comparison=comparison,
    n=length(status),
    events=n_event,
    nonevents=n_nonevent,
    eval_time_days=eval_time,
    old_mean_risk=mean(old_risk),
    new_mean_risk=mean(new_risk),
    nri=nri,
    nri_se=nri_se,
    nri_lower=nri - 1.96*nri_se,
    nri_upper=nri + 1.96*nri_se,
    nri_p=normal_p(nri, nri_se),
    nri_event=nri_event,
    nri_event_se=se_event,
    nri_event_lower=nri_event - 1.96*se_event,
    nri_event_upper=nri_event + 1.96*se_event,
    nri_event_p=normal_p(nri_event, se_event),
    event_up=event_up,
    event_down=event_down,
    nri_nonevent=nri_nonevent,
    nri_nonevent_se=se_nonevent,
    nri_nonevent_lower=nri_nonevent - 1.96*se_nonevent,
    nri_nonevent_upper=nri_nonevent + 1.96*se_nonevent,
    nri_nonevent_p=normal_p(nri_nonevent, se_nonevent),
    nonevent_down=nonevent_down,
    nonevent_up=nonevent_up,
    method_note="continuous NRI from Cox absolute risk; analytic SE from reclassification proportions",
    stringsAsFactors=FALSE
  )
}

nri_results <- list()
for (method_name in methods) {
  cat("NRI models for", method_name, "\n")
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
  model_cols <- c("time", "status", covariates, method_score, protein_term_names)
  model_data <- as.data.frame(dat[, model_cols, with=FALSE])
  model_data <- model_data[complete.cases(model_data),]
  model_data <- droplevels(model_data)
  method_covariates <- usable_covariates(model_data, covariates)
  covariate_terms <- paste(method_covariates, collapse="+")
  protein_terms <- paste(protein_term_names, collapse="+")

  base_fit <- coxph(as.formula(paste("Surv(time, status) ~", covariate_terms)), data=model_data)
  pfs_plus <- coxph(as.formula(paste("Surv(time, status) ~", covariate_terms, "+", method_score)), data=model_data)
  protein_plus <- coxph(
    as.formula(paste("Surv(time, status) ~", covariate_terms, "+", protein_terms)),
    data=model_data,
    control=coxph.control(iter.max=30)
  )

  eval_time <- max(model_data$time, na.rm=TRUE)
  nri_results[[length(nri_results)+1]] <- continuous_nri_row(
    method_name,
    "traditional_covariates_vs_traditional_plus_stable_pfs",
    base_fit,
    pfs_plus,
    model_data,
    eval_time
  )
  nri_results[[length(nri_results)+1]] <- continuous_nri_row(
    method_name,
    "traditional_covariates_vs_traditional_plus_stable_proteins",
    base_fit,
    protein_plus,
    model_data,
    eval_time
  )
}

nri_result_dt <- rbindlist(nri_results, fill=TRUE)
fwrite(nri_result_dt, nri_out_path)
write_json_simple(list(
  method="Continuous NRI for stable nested-CV LASSO/ElasticNet PFS and stable proteins; full follow-up after exclusion of PD events within 2 years.",
  base_cindex_script=base_script,
  output=nri_out_path,
  methods=methods,
  comparisons=c(
    "traditional_covariates_vs_traditional_plus_stable_pfs",
    "traditional_covariates_vs_traditional_plus_stable_proteins"
  ),
  inference="analytic normal-approximation CI from reclassification proportions; not bootstrap"
), nri_summary_path)

print(nri_result_dt)
cat("NRI output:", nri_out_path, "\n")
