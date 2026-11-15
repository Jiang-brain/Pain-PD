import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from common import BRANCH, OUT, REFERENCE_DIR, SENSITIVITY_DIR, md5


def main():
    checks = []

    def check(name, passed, detail=""):
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

    prepared = json.loads((OUT / "preparation_manifest.json").read_text())
    new_summary = json.loads((OUT / "new_elastic_net_summary.json").read_text())
    old_performance = json.loads((OUT / "original_elastic_net_reconstruction_performance.json").read_text())
    for filename, expected in [
        ("original_elastic_net_fold_reconstruction_checks.csv", 10),
        ("original_stable_cox_reproduction_checks.csv", 5),
        ("original_stable_nri_reproduction_checks.csv", 2),
        ("original_shared_cox_reproduction_checks.csv", 9),
        ("lightgbm_baseline/reproduction_checks.csv", 9),
    ]:
        table = pd.read_csv(OUT / filename)
        check(f"reproduce_{filename}", len(table) == expected and table.passed.all())
    hashes = pd.read_csv(OUT / "input_hashes.csv")
    check("input_files_unchanged", all(md5(row.path) == row.md5 for row in hashes.itertuples()))
    snapshots = list((BRANCH / "source_reference").iterdir())
    check("original_source_snapshots_unchanged", all(md5(path) == md5(REFERENCE_DIR / path.name) for path in snapshots if path.is_file()))
    coefficients = pd.read_csv(OUT / "new_elastic_net_fold_coefficients.csv")
    fold_metrics = pd.read_csv(OUT / "new_elastic_net_fold_metrics.csv")
    check("all_10_complete_outer_folds", set(coefficients.fold) == set(range(1, 11)) and
          len(coefficients) == 10 * prepared["n_proteins"] and len(fold_metrics) == 10)
    check("Elastic_Net_L1_and_L2", fold_metrics.l1_ratio.between(0, 1, inclusive="neither").all())
    check("no_ElasticNet_fit_warnings", fold_metrics.warning_count.eq(0).all(), int(fold_metrics.warning_count.sum()))
    oof = pd.read_csv(OUT / "new_elastic_net_oof_predictions.csv")
    folds = pd.read_csv(OUT / "outer_fold_assignments.csv")
    check("OOF_participants_and_folds_match_original", len(oof) == prepared["n"] and not oof.eid.duplicated().any() and
          np.array_equal(oof.eid, folds.eid) and np.array_equal(oof.outer_fold, folds.outer_fold))
    target = np.load(OUT / "target_new_log.npy")
    check("OOF_target_is_new_pain_score", np.allclose(target, oof.target_log, rtol=0, atol=1e-12))
    check("OOF_predictions_finite", np.isfinite(oof.oof_prediction_log).all())
    new_coef = pd.read_csv(OUT / "new_stable_elastic_net_coefficients.csv")
    selected = new_coef[new_coef.stable_all_10_folds]
    check("stable_selection_all_10_folds", selected.selection_count.eq(10).all() and len(selected) == new_summary["n_stable_new"])
    check("stable_signature_not_forced_to_374", len(selected) > 0, f"Observed new stable count: {len(selected)}")

    comparisons = []

    def compare_table(old, new, keys, metric, analysis):
        merged = old.merge(new, on=keys, suffixes=("_original", "_new"), validate="one_to_one")
        for _, row in merged.iterrows():
            comparisons.append({"analysis": analysis, **{key: row[key] for key in keys}, "metric": metric,
                                "original": row[f"{metric}_original"], "new": row[f"{metric}_new"],
                                "difference_new_minus_original": row[f"{metric}_new"] - row[f"{metric}_original"]})
        return merged

    stable_old = pd.read_csv(OUT / "stable_original374/cindex.csv")
    stable_new = pd.read_csv(OUT / "stable_new_score/cindex.csv")
    stable_comparison = compare_table(stable_old, stable_new, ["model_type"], "c_index", "stable_panel_Cox")
    check("stable_Cox_same_participants_events", stable_comparison.n_original.eq(stable_comparison.n_new).all() and
          stable_comparison.events_original.eq(stable_comparison.events_new).all())
    shared = pd.read_csv(OUT / "shared_prediction_cindex.csv")
    shared_comparison = compare_table(shared[shared.case.eq("shared_original25")], shared[shared.case.eq("shared_new24")],
                                      ["task", "model"], "c_index", "shared_panel_Cox")
    check("shared_Cox_same_participants_events", shared_comparison.n_original.eq(shared_comparison.n_new).all() and
          shared_comparison.events_original.eq(shared_comparison.events_new).all())
    stable_nri_old = pd.read_csv(OUT / "stable_original374/nri.csv")
    stable_nri_new = pd.read_csv(OUT / "stable_new_score/nri.csv")
    compare_table(stable_nri_old, stable_nri_new, ["comparison"], "nri", "stable_panel_NRI")
    shared_nri = pd.read_csv(OUT / "shared_prediction_nri.csv")
    compare_table(shared_nri[shared_nri.case.eq("shared_original25")], shared_nri[shared_nri.case.eq("shared_new24")],
                  ["task"], "nri", "shared_panel_NRI")

    lgbm_old = pd.read_csv(OUT / "lightgbm_baseline/cv_metrics.csv")
    lgbm_new = pd.read_csv(OUT / "lightgbm_new/cv_metrics.csv")
    stable_mapping = {"stable_elasticnet_374": "new_stable_panel",
                       "covariates_plus_stable_elasticnet_374": "covariates_plus_new_stable_panel"}
    lgbm_old.feature_set = lgbm_old.feature_set.replace(stable_mapping)
    stable_lgbm = lgbm_new[lgbm_new.feature_set.isin(lgbm_old.feature_set)].copy()
    for metric in ["auc", "auprc", "brier"]:
        merged = compare_table(lgbm_old, stable_lgbm, ["task", "feature_set"], metric, "stable_panel_LightGBM")
        check(f"stable_LightGBM_same_n_events_{metric}", merged.n_total_original.eq(merged.n_total_new).all() and
              merged.events_total_original.eq(merged.events_total_new).all())
    old_shared = pd.read_csv(REFERENCE_DIR / "pd_prediction_intersection_lgbm_cv_metrics.csv")
    old_shared.feature_set = old_shared.feature_set.replace({"fdr_same_direction_25": "new_shared24",
        "covariates_plus_fdr_same_direction_25": "covariates_plus_new_shared24"})
    new_shared = lgbm_new[lgbm_new.feature_set.isin(old_shared.feature_set)]
    for metric in ["auc", "auprc", "brier"]:
        compare_table(old_shared, new_shared, ["task", "feature_set"], metric, "shared_panel_LightGBM")
    cached_predictions = pd.read_csv(REFERENCE_DIR / "pd_prediction_intersection_lgbm_oof_predictions.csv",
                                     usecols=["task", "feature_set", "outcome", "predicted_risk"])
    saved_metrics = pd.read_csv(REFERENCE_DIR / "pd_prediction_intersection_lgbm_cv_metrics.csv")
    cached_checks = []
    for (task, feature_set), group in cached_predictions.groupby(["task", "feature_set"]):
        metrics = {"auc": roc_auc_score(group.outcome, group.predicted_risk),
                   "auprc": average_precision_score(group.outcome, group.predicted_risk),
                   "brier": brier_score_loss(group.outcome, group.predicted_risk)}
        ref = saved_metrics[saved_metrics.task.eq(task) & saved_metrics.feature_set.eq(feature_set)].iloc[0]
        difference = max(abs(metrics[name] - ref[name]) for name in metrics)
        cached_checks.append({"task": task, "feature_set": feature_set,
                              "max_metric_difference": difference, "passed": difference < 1e-10})
    pd.DataFrame(cached_checks).to_csv(OUT / "original_shared_lightgbm_cached_metric_checks.csv", index=False)
    check("original_shared_LightGBM_cached_metrics_consistent", len(cached_checks) == 9 and all(row["passed"] for row in cached_checks))
    predictions = pd.read_csv(OUT / "lightgbm_new/oof_predictions.csv")
    check("new_PD_predictions_finite_bounded", np.isfinite(predictions.predicted_risk).all() and predictions.predicted_risk.between(0, 1).all())
    check("new_PD_OOF_unique_per_person_task_model", not predictions.duplicated(["eid", "task", "feature_set"]).any())
    check("new_PD_OOF_all_10_folds", predictions.groupby(["task", "feature_set"]).fold.nunique().eq(10).all())
    expected_rows = lgbm_new.set_index(["task", "feature_set"]).n_total
    observed_rows = predictions.groupby(["task", "feature_set"]).size()
    check("new_PD_OOF_complete", observed_rows.sort_index().equals(expected_rows.sort_index().rename(None)))
    comparison = pd.DataFrame(comparisons)
    comparison.to_csv(OUT / "prediction_comparison.csv", index=False)
    covariate_audit = pd.read_csv(OUT / "legacy_Cox_covariate_variation_audit.csv")
    constant_covariates = covariate_audit.loc[covariate_audit.unique_nonmissing.lt(2), "covariate"].tolist()
    pd.DataFrame([
        {"issue": "Legacy Cox recoding collapses covariates", "severity": "high",
         "detail": "Constant and automatically dropped: " + ", ".join(constant_covariates)},
        {"issue": "Inner preprocessing is not fully nested", "severity": "high", "detail": "Outer train imputer/scaler precede ElasticNetCV inner splitting."},
        {"issue": "Protein selection precedes PD CV", "severity": "high", "detail": "PD feature selection is not independently repeated inside each PD training fold."},
        {"issue": "Cox apparent performance", "severity": "high", "detail": "C-index and NRI are fitted-cohort estimates, not external validation."},
        {"issue": "Near-term Cox warning", "severity": "high", "detail": "Some legacy Cox models warn of potentially infinite covariate coefficients."},
        {"issue": "Long-window wording", "severity": "medium", "detail": "Original LightGBM after5y is not bounded above at ten years."},
    ]).to_csv(OUT / "methodological_warnings.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "validation_checks.csv", index=False)
    if not all(row["passed"] for row in checks):
        raise ValueError("Some prediction validation checks failed; inspect validation_checks.csv")

    def rows(analysis, metric, selected_models=None):
        table = comparison[comparison.analysis.eq(analysis) & comparison.metric.eq(metric)]
        if selected_models is not None:
            column, values = selected_models
            table = table[table[column].isin(values)]
        result = []
        for _, row in table.iterrows():
            label = " / ".join(str(row[key]) for key in ["task", "model_type", "model", "feature_set", "comparison"] if key in row and pd.notna(row[key]))
            result.append(f"| {label} | {row.original:.4f} | {row.new:.4f} | {row.difference_new_minus_original:+.4f} |")
        return result

    header = ["| Model / window | Original | New score | Difference |", "| --- | ---: | ---: | ---: |"]
    cox_combined = stable_comparison.loc[stable_comparison.model_type.eq("traditional_plus_stable_proteins")].iloc[0]
    lgbm_combined = comparison[comparison.analysis.eq("stable_panel_LightGBM") & comparison.metric.eq("auc") &
                               comparison.feature_set.eq("covariates_plus_new_stable_panel")]
    headline = [
        "## Main Findings",
        f"- New stable panel: {new_summary['n_stable_new']} proteins, not an exact reproduction of the original 374-protein list.",
        f"- Apparent Cox C-index for panel plus covariates: {cox_combined.c_index_original:.4f} to {cox_combined.c_index_new:.4f}; numerical change is small.",
    ]
    for row in lgbm_combined.itertuples(index=False):
        headline.append(f"- LightGBM combined OOF AUC, {row.task}: {row.original:.4f} to {row.new:.4f} ({row.difference_new_minus_original:+.4f}).")
    headline += [
        "- Cox results are close; combined LightGBM AUC is lower in all three windows. This is not exact numerical replication or a statistical equivalence result.",
        "- The single PPS has near-chance standalone PD concordance; the higher values belong to multivariable protein panels.",
        "- Verified original Cox recoding and validation-design issues remain; numerical checks do not certify the original method as correct.", "",
    ]
    report = [
        "# Prediction Sensitivity With New Pain Scores", "", "## Material Passport",
        "- Status: COMPLETED; original numerical baselines reproduced and output checks passed.",
        "- Methodological warnings remain, especially the verified original Cox deprivation/PRS recoding bug.",
        "- Sources: unchanged PA5 prediction scripts/results and branch-07 non-proteomic-derived pain scores.",
        "- Original data, manuscript, main scripts, and prior sensitivity outputs were not overwritten.", "",
        *headline,
        "## Target and Model Re-estimation",
        f"The same {prepared['n']:,} participants and {prepared['n_proteins']:,} proteins were used for Elastic Net.",
        "Only the pain-score target changed. Outer and inner folds, seeds, log(score+0.1) transform, alpha grid,",
        "and L1/L2 penalty candidates were held fixed. Compute-only job/thread limits were used for memory control.",
        f"The new stable panel contained {new_summary['n_stable_new']} proteins (original: 374):",
        f"{new_summary['retained']} retained, {new_summary['lost']} lost, and {new_summary['gained']} gained.",
        f"Original vs new pain-score correlation: {prepared['pain_score_correlation_old_new']:.6f}.", "",
        "### Genuine OOF Pain-Score Prediction",
        f"Original outer-fold OOF Pearson r: {old_performance['original_oof_target_log']['pearson_r']:.4f};",
        f"new outer-fold OOF Pearson r: {new_summary['new_elastic_net_oof_target_log']['pearson_r']:.4f}.",
        "These refer to predicting log-transformed pain, not predicting PD.",
        "The original and new pain targets differ; these r values are descriptive comparisons, not a test that one target is easier or clinically superior.",
        f"Stable mean-coefficient apparent Pearson r: original {old_performance['original_stable_score_target_log']['pearson_r']:.4f},",
        f"new {new_summary['new_stable_score_target_log_apparent']['pearson_r']:.4f}. This is not the OOF result.", "",
        "## Stable-Panel Cox Apparent C-Index", *header, *rows("stable_panel_Cox", "c_index"), "",
        "The protein-panel models have one predictor per protein. The single PPS models use just one proteomic",
        "pain-score predictor. The original 0.831/0.888 values describe the panel, not the single PPS.", "",
        "## Shared-Panel Cox Apparent C-Index", *header,
        *rows("shared_panel_Cox", "c_index", ("model", ["shared_proteins_alone", "traditional_plus_shared_proteins"])), "",
        "The shared panel changed from 25 to 24 by dropping COMP. Participant cohorts were held fixed.", "",
        "## Continuous NRI", *header, *rows("stable_panel_NRI", "nri"), *rows("shared_panel_NRI", "nri"), "",
        "## Stable-Panel LightGBM OOF AUC", *header, *rows("stable_panel_LightGBM", "auc"), "",
        "## Shared-Panel LightGBM OOF AUC", *header, *rows("shared_panel_LightGBM", "auc"), "",
        "AUPRC and Brier comparisons are in prediction_comparison.csv. Full OOF predictions, fold results,",
        "participant IDs, and cohort/event counts are retained under lightgbm_baseline/ and lightgbm_new/.", "",
        "## Replication and Limitations",
        "- IMPORTANT VERIFIED SOURCE BUG: sequential recoding in the original R Cox prediction code collapses",
        "  deprivation and PRS to one observed level. usable_covariates then drops them. The original R models",
        "  therefore do NOT actually adjust for these two variables, despite listing them as candidate covariates.",
        "  This controlled sensitivity preserves that behavior to isolate the pain-target change. A corrected",
        "  covariate analysis must be a separately labelled rerun, not silently mixed with this comparison.",
        "- Original stable and shared Cox C-indices, stable NRI, original stable PPS, and original fold predictions were reproduced.",
        "- The original reported PPS-pain r near 0.303 comes from the stable mean-coefficient population fit, not the genuine",
        "  outer-fold OOF prediction (r near 0.256). Manuscript wording should distinguish these quantities.",
        "- Original 374-protein LightGBM models were refitted and their AUC/AUPRC/Brier checked numerically.",
        "- Original shared-25 LightGBM metrics were verified against saved OOF predictions (not a second refit of the 25-protein models).",
        "- No new count was forced to equal the original 374. No PCA model, signature, or loading was used.",
        "- Exact legacy specification was preserved for comparability: imputation/scaling precede INNER ElasticNetCV splitting",
        "  within each outer training fold. This is not a newly corrected fully nested preprocessing analysis.",
        "- Shared/stable panel selection precedes PD CV. Selection was not repeated independently within PD folds.",
        "- Cox C-indices/NRI are apparent fitted-cohort estimates, not held-out or externally validated performance.",
        "- The original near-term shared Cox models warn of potentially infinite coefficients in some covariates; numerical",
        "  reproduction does not remove this concern. The original technical factors and covariates were not silently changed.",
        "- Original LightGBM after5y means all follow-up after 5 years, NOT a 5-10-year window. Cox uses within-5/10y horizons.",
        "- Original analytic NRI ignores censoring adjustments in reclassification proportions; threshold-based classification",
        "  metrics use data-derived Youden cutoffs. Clinical calibration, utility, and external validity remain unestablished.",
        "  LightGBM preserves the original binary outcome construction rather than fitting a censoring-adjusted survival model.",
        "- A sensitivity comparison showing similar numerical performance is not a statistical equivalence test or causal evidence.",
        f"- All {len(checks)} implemented output/reproduction checks passed. Source files and input hashes were verified.", "",
        "## Files and Privacy",
        "See prediction_comparison.csv, validation_checks.csv, new_elastic_net_summary.json, stable_protein_membership_comparison.csv,",
        "new_stable_elastic_net_coefficients.csv, new_stable_elastic_net_pfs.csv, and the model-specific subdirectories.",
        "Participant-level files and arrays must remain private and excluded from GitHub.",
    ]
    (OUT / "PREDICTION_SENSITIVITY_RESULTS.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(json.dumps({"checks_passed": len(checks), "new_stable_proteins": len(selected),
                      "report": str(OUT / "PREDICTION_SENSITIVITY_RESULTS.md")}, indent=2), flush=True)
    print(comparison.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
