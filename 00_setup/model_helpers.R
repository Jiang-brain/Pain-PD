source("00_setup/config.R")

usable_covariates <- function(model_data, covariates) {
  keep <- vapply(covariates, function(v) {
    x <- model_data[[v]]
    if (is.factor(x) || is.character(x)) return(nlevels(droplevels(factor(x))) >= 2)
    length(unique(x[!is.na(x)])) >= 2
  }, logical(1))
  covariates[keep]
}

cox_fit_row <- function(label, fit, ref_fit = NULL, extra_terms = "") {
  fit_summary <- summary(fit)
  cidx <- unname(fit_summary$concordance[1])
  se <- unname(fit_summary$concordance["se(C)"])
  ref <- if (is.null(ref_fit)) NA_real_ else unname(summary(ref_fit)$concordance[1])
  data.frame(
    model = label,
    n = fit$n,
    events = fit$nevent,
    extra_terms = extra_terms,
    c_index = cidx,
    c_index_se = se,
    reference_c_index = ref,
    c_index_delta = ifelse(is.na(ref), NA_real_, cidx - ref),
    stringsAsFactors = FALSE
  )
}

predict_cox_absolute_risk <- function(fit, model_data, eval_time) {
  bh <- basehaz(fit, centered = FALSE)
  bh <- bh[is.finite(bh$time) & is.finite(bh$hazard), ]
  bh <- bh[bh$time <= eval_time, ]
  h0 <- if (nrow(bh) == 0) 0 else bh$hazard[nrow(bh)]
  lp <- predict(fit, newdata = model_data, type = "lp", reference = "zero")
  risk <- 1 - exp(-h0 * exp(lp))
  pmin(pmax(as.numeric(risk), 0), 1)
}

continuous_nri_row <- function(task, comparison, old_fit, new_fit, model_data, eval_time, status_col = "status_horizon") {
  old_risk <- predict_cox_absolute_risk(old_fit, model_data, eval_time)
  new_risk <- predict_cox_absolute_risk(new_fit, model_data, eval_time)
  status <- model_data[[status_col]]
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
    2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  }
  data.frame(
    task = task,
    comparison = comparison,
    n = length(status),
    events = n_event,
    nonevents = n_nonevent,
    eval_time_days = eval_time,
    old_mean_risk = mean(old_risk),
    new_mean_risk = mean(new_risk),
    nri = nri,
    nri_se = nri_se,
    nri_lower = nri - 1.96 * nri_se,
    nri_upper = nri + 1.96 * nri_se,
    nri_p = normal_p(nri, nri_se),
    nri_event = nri_event,
    nri_event_se = se_event,
    nri_event_lower = nri_event - 1.96 * se_event,
    nri_event_upper = nri_event + 1.96 * se_event,
    nri_event_p = normal_p(nri_event, se_event),
    event_up = event_up,
    event_down = event_down,
    nri_nonevent = nri_nonevent,
    nri_nonevent_se = se_nonevent,
    nri_nonevent_lower = nri_nonevent - 1.96 * se_nonevent,
    nri_nonevent_upper = nri_nonevent + 1.96 * se_nonevent,
    nri_nonevent_p = normal_p(nri_nonevent, se_nonevent),
    nonevent_down = nonevent_down,
    nonevent_up = nonevent_up,
    stringsAsFactors = FALSE
  )
}

make_horizon_task <- function(data, horizon_years) {
  load_packages(c("data.table"))
  horizon <- 365.25 * horizon_years
  task <- data.table::copy(data)
  task[, time_horizon := pmin(time, horizon)]
  task[, status_horizon := as.integer(status == 1 & time <= horizon)]
  task[time >= horizon | status_horizon == 1]
}

make_full_followup_task <- function(data) {
  load_packages(c("data.table"))
  task <- data.table::copy(data[time > FOLLOWUP_EXCLUSION_DAYS])
  task[, time_horizon := time - FOLLOWUP_EXCLUSION_DAYS]
  task[, status_horizon := as.integer(status == 1 & time > FOLLOWUP_EXCLUSION_DAYS)]
  task
}

