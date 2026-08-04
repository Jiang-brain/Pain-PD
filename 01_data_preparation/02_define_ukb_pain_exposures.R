source("00_setup/ukb_helpers.R")
load_packages(c("data.table"))

pain <- make_pain_duration(UKB_DATA_FILE)
site_cols <- c("headache", "facial", "neck_shoulder", "back", "stomach", "hip", "knee")

pain[, "body_pain_reported"] <- pain$all_over_body == 2 | pain$all_over_body == 1
pain_no_body <- drop_body_pain_without_sites(pain)

pain_no_body$chronic_site_count <- rowSums(pain_no_body[, site_cols] == 2, na.rm = TRUE)
pain_no_body$acute_site_count <- rowSums(pain_no_body[, site_cols] == 1, na.rm = TRUE)
pain_no_body$any_site_count <- rowSums(as.data.frame(lapply(pain_no_body[, site_cols], function(x) x %in% c(1, 2))), na.rm = TRUE)
pain_no_body$chronic_pain_status <- factor(ifelse(
  pain_no_body$chronic_site_count == 0, 0,
  ifelse(pain_no_body$chronic_site_count == 1, 1, 2)
))
pain_no_body$acute_pain_status <- factor(ifelse(
  pain_no_body$acute_site_count == 0, 0,
  ifelse(pain_no_body$acute_site_count == 1, 1, 2)
))

data.table::fwrite(as.data.table(pain_no_body), result_path("ukb_pain_exposures.csv"))

cat("Rows after excluding all-over-body pain without site information:", nrow(pain_no_body), "\n")
cat("Output:", result_path("ukb_pain_exposures.csv"), "\n")
