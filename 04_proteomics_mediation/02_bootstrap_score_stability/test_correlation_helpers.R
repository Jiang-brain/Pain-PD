source("04_proteomics_mediation/02_bootstrap_score_stability/00_bootstrap_helpers.R")
set.seed(19)
patterns <- matrix(c(0,0,0, 1,0,0, 0,1,0, 0,0,1, 1,1,0, 1,0,1), ncol = 3, byrow = TRUE)
frequency <- c(50, 8, 15, 2, 11, 3)
weights <- matrix(c(0,.4,.6, .1,.2,.7, .3,.3,.4, 0,0,1), nrow = 3)
expanded <- patterns[rep(seq_len(nrow(patterns)), frequency), , drop = FALSE]
scores <- expanded %*% weights
group_scores <- patterns %*% weights
ranks <- vapply(seq_len(ncol(weights)), function(j) weighted_midranks(group_scores[,j], frequency),
                numeric(nrow(patterns)))
pearson_difference <- max(abs(score_covariance_correlation(expanded, weights) - cor(scores)))
spearman_difference <- max(abs(weighted_column_correlation(ranks, frequency) - cor(scores, method = "spearman")))
stopifnot(pearson_difference < 1e-12, spearman_difference < 1e-12)
stopifnot(nrow(summarize_correlations(c(.5,.8,.9), "test", "test", "test")) == 1L)
cat("Synthetic mathematical unit tests passed; Pearson difference:", pearson_difference,
    "; Spearman difference:", spearman_difference, "\n")
