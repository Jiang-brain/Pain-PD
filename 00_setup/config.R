PROJECT_ROOT <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = FALSE)
CODE_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

DATA_DIR <- normalizePath(Sys.getenv("PD_PAIN_DATA_DIR", file.path(CODE_ROOT, "data")),
                           winslash = "/", mustWork = FALSE)
RESULTS_DIR <- normalizePath(Sys.getenv("PD_PAIN_RESULTS_DIR", file.path(CODE_ROOT, "results")),
                              winslash = "/", mustWork = FALSE)
FIGURE_DIR <- normalizePath(file.path(RESULTS_DIR, "figures"), winslash = "/", mustWork = FALSE)

dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)

input_path <- function(filename) file.path(DATA_DIR, filename)
result_path <- function(filename) file.path(RESULTS_DIR, filename)
figure_path <- function(filename) file.path(FIGURE_DIR, filename)

find_input <- function(...) {
  candidates <- c(...)
  candidates <- c(candidates, file.path(DATA_DIR, basename(candidates)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("None of these input files exist: ", paste(unique(candidates), collapse = ", "))
  }
  existing[1]
}

load_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install them before running this script.")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

UKB_DATA_FILE <- input_path("PD_data.csv")
PROTEIN_MATRIX_FILE <- input_path("protein_new_missing60_filtered.csv")
PROTEIN_COVARIATE_FILE <- input_path("PD_data_matched_protein_covariates.csv")
PROTEIN_TIME_FILE <- input_path("PD_data_matched_protein_times.csv")
COMPOSITE_SCORE_FILE <- result_path("composite_pain_scores.csv")
INTERSECTION_FILE <- input_path("protein_intersection_FDR05.csv")
PPMI_RESULTS_FILE <- input_path("results_PPMI.csv")

FOLLOWUP_EXCLUSION_DAYS <- 365.25 * 2
FOLLOWUP_5Y_DAYS <- 365.25 * 5
FOLLOWUP_10Y_DAYS <- 365.25 * 10

MAIN_COVARIATES <- c(
  "age", "gender", "ethnicity", "deprivation", "income", "university_education",
  "smoking", "alcohol", "sedentary", "obesity", "prs", paste0("comp_", 1:10)
)

GENETIC_PC_COVARIATES <- paste0("comp_", 1:10)
PROTEOMIC_TECHNICAL_COVARIATES <- c("protein_cov_zscore", "season_final", "Batch", "protein_time")
