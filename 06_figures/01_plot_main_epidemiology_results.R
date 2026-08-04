source("00_setup/config.R")
load_packages(c("data.table", "ggplot2"))

plot_forest <- function(input_file, output_file, label_col = "site", group_col = NULL) {
  if (!file.exists(input_file)) stop("Missing result file: ", input_file)
  dat <- data.table::fread(input_file)
  if (!all(c("HR", "lower", "upper", "p_value") %in% names(dat))) {
    stop("Input must contain HR, lower, upper, and p_value columns: ", input_file)
  }
  dat$label <- if (label_col %in% names(dat)) dat[[label_col]] else seq_len(nrow(dat))
  if (!is.null(group_col) && group_col %in% names(dat)) {
    dat$label <- paste(dat[[group_col]], dat$label, sep = ": ")
  }
  dat$label <- factor(dat$label, levels = rev(unique(dat$label)))
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = HR, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.2) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Hazard ratio", y = NULL) +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(output_file, p, width = 7, height = max(4, 0.28 * nrow(dat) + 1.5), dpi = 300)
  output_file
}

outputs <- c()
site_file <- result_path("ukb_site_specific_pain_cox.csv")
if (file.exists(site_file)) {
  outputs <- c(outputs, plot_forest(site_file, figure_path("fig_site_specific_pain_cox.png"), label_col = "site", group_col = "duration"))
}

status_file <- result_path("ukb_single_multisite_pain_cox.csv")
if (file.exists(status_file)) {
  outputs <- c(outputs, plot_forest(status_file, figure_path("fig_single_multisite_pain_cox.png"), label_col = "contrast", group_col = "exposure"))
}

tertile_file <- result_path("ukb_composite_pain_score_tertile_cox.csv")
if (file.exists(tertile_file)) {
  outputs <- c(outputs, plot_forest(tertile_file, figure_path("fig_composite_pain_score_tertile_cox.png"), label_col = "contrast"))
}

if (length(outputs) == 0) {
  cat("No eligible epidemiology result files found in ", RESULTS_DIR, "\n", sep = "")
} else {
  cat("Figure outputs:\n")
  cat(paste(outputs, collapse = "\n"), "\n")
}

