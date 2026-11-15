# Run from code/. This branch never sources or writes main-analysis outputs.
CODE_ROOT_SENS <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
BRANCH_DIR <- file.path(CODE_ROOT_SENS, "04_proteomics_mediation/01_independent_pain_score")
if (!dir.exists(BRANCH_DIR)) stop("Run this script from the code/ directory.")
.libPaths(c(file.path(CODE_ROOT_SENS, "R_libs"), .libPaths()))
for (pkg in c("data.table", "survival")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
}
library(data.table)
library(survival)

SENS_DATA_DIR <- normalizePath(Sys.getenv("PD_PAIN_DATA_DIR", "data"), winslash = "/", mustWork = TRUE)
SENS_OUT_DIR <- file.path(Sys.getenv("PD_PAIN_RESULTS_DIR", "results"), "sensitivity_nonproteomic_weights")
dir.create(SENS_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SENS_OUT_DIR <- normalizePath(SENS_OUT_DIR, winslash = "/", mustWork = TRUE)
if (identical(SENS_OUT_DIR, SENS_DATA_DIR)) stop("Output directory must differ from input directory.")
sens_input <- function(name) {
  path <- file.path(SENS_DATA_DIR, name)
  if (!file.exists(path)) stop("Missing input: ", path)
  path
}
sens_output <- function(name) file.path(SENS_OUT_DIR, name)
reference_path <- function(name) file.path(BRANCH_DIR, "source_reference", name)
assigned_name <- function(expr) {
  if (is.call(expr) && as.character(expr[[1]]) %in% c("<-", "=") && is.symbol(expr[[2]])) {
    return(as.character(expr[[2]]))
  }
  ""
}
is_libpaths <- function(expr) is.call(expr) && identical(expr[[1]], as.name(".libPaths"))
assert_unique_eid <- function(dat, label) {
  if (anyNA(dat$eid) || anyDuplicated(dat$eid)) stop(label, ": missing or duplicate eid.")
}
write_provenance <- function(paths, filename) {
  fwrite(data.table(path = normalizePath(paths, winslash = "/", mustWork = TRUE),
                    md5 = unname(tools::md5sum(paths))), sens_output(filename))
}
