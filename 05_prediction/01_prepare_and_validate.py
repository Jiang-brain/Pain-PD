import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.model_selection import KFold

from common import BRANCH, DATA_DIR, NEW_SCORE_PATH, OUT, REFERENCE_DIR, load_reference, md5


def main():
    original = load_reference("build_stable_nested_elastic_lasso_pfs.py", "original_elastic_net")
    original.PROTEIN_PATH = REFERENCE_DIR / "protein_new_missing60_filtered.csv"
    original.SCORE_PATH = REFERENCE_DIR / "composite_pain_scores.csv"
    print("Reading original prediction cohort and all proteins...", flush=True)
    data, protein_names = original.read_inputs()
    new_scores = pd.read_csv(NEW_SCORE_PATH, usecols=["eid", "pain_score", "weighted_pain_status"])
    if data.eid.duplicated().any() or new_scores.eid.duplicated().any():
        raise ValueError("Duplicate participant IDs")
    matched = data[["eid", "pain_score"]].merge(new_scores, on="eid", how="left", validate="one_to_one",
                                                suffixes=("_original", "_new"), sort=False)
    if not np.array_equal(matched.eid, data.eid) or matched.pain_score_new.isna().any():
        raise ValueError("New scores must cover every original Elastic Net participant in the same order")
    x = data[protein_names].to_numpy(dtype=np.float32)
    np.save(OUT / "protein_matrix.npy", x, allow_pickle=False)
    old_log = np.log(matched.pain_score_original.to_numpy(dtype=float) + original.TARGET_OFFSET)
    new_log = np.log(matched.pain_score_new.to_numpy(dtype=float) + original.TARGET_OFFSET)
    np.save(OUT / "target_original_log.npy", old_log, allow_pickle=False)
    np.save(OUT / "target_new_log.npy", new_log, allow_pickle=False)
    matched["weighted_pain_status_new"] = new_scores.set_index("eid").loc[matched.eid, "weighted_pain_status"].to_numpy()
    matched.to_csv(OUT / "matched_pain_scores.csv", index=False)
    (OUT / "protein_names.json").write_text(json.dumps(protein_names, indent=2), encoding="utf-8")
    cv = KFold(n_splits=10, shuffle=True, random_state=original.RANDOM_STATE)
    folds = np.zeros(len(data), dtype=np.int16)
    for fold, (_, test) in enumerate(cv.split(x), 1):
        folds[test] = fold
    pd.DataFrame({"eid": data.eid, "outer_fold": folds}).to_csv(OUT / "outer_fold_assignments.csv", index=False)

    coefficients = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_stable_nested_coefficients.csv")
    stable = coefficients[coefficients.method.eq("elasticnet") & coefficients.stable_all_10_folds].copy()
    if len(stable) != 374:
        raise ValueError("Original reference does not contain the reported 374-protein Elastic Net panel")
    positions = {name: i for i, name in enumerate(protein_names)}
    predicted = np.full(len(data), float(stable.mean_intercept_raw.mean()))
    for row in stable.itertuples(index=False):
        values = x[:, positions[row.protein]].astype(float)
        values[np.isnan(values)] = row.mean_impute_value
        predicted += values * row.mean_coef_raw_imputed
    saved = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_stable_nested_pfs.csv").set_index("eid").loc[data.eid]
    maximum_difference = float(np.max(np.abs(predicted - saved.stable_elasticnet_pfs_log.to_numpy())))
    if maximum_difference >= 1e-10:
        raise ValueError(f"Original stable score fails reconstruction: {maximum_difference}")

    fold_coefficients = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_nested_elasticnet_fold_coefficients.csv")
    reference_fold_metrics = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_nested_elasticnet_fold_metrics.csv")
    oof = np.full(len(data), np.nan)
    validation_rows = []
    for fold in range(1, 11):
        rows = np.flatnonzero(folds == fold)
        coef = fold_coefficients[fold_coefficients.fold.eq(fold)]
        prediction = np.full(len(rows), float(coef.model_intercept_raw_imputed_space.iloc[0]))
        for row in coef[coef.selected].itertuples(index=False):
            values = x[rows, positions[row.protein]].astype(float)
            values[np.isnan(values)] = row.protein_median_impute
            prediction += values * row.coef_raw_imputed_protein_to_log_pain
        oof[rows] = prediction
        metric = original.evaluate(old_log[rows], prediction)
        ref = reference_fold_metrics.loc[reference_fold_metrics.fold.eq(fold)].iloc[0]
        difference = max(abs(metric[key] - float(ref[key])) for key in ["r2", "rmse", "mae", "pearson_r", "spearman_r"])
        validation_rows.append({"fold": fold, "max_metric_difference": difference, "passed": difference < 1e-5})
    if not all(row["passed"] for row in validation_rows):
        raise ValueError("Original OOF prediction does not reproduce from saved fold coefficients")
    pd.DataFrame(validation_rows).to_csv(OUT / "original_elastic_net_fold_reconstruction_checks.csv", index=False)
    pd.DataFrame({"eid": data.eid, "outer_fold": folds, "original_oof_prediction_log": oof}).to_csv(
        OUT / "original_elastic_net_oof_reconstructed.csv", index=False)
    original_performance = {
        "original_oof_target_log": original.evaluate(old_log, oof),
        "original_stable_score_target_log": original.evaluate(old_log, predicted),
        "frozen_original_stable_score_new_target_log": original.evaluate(new_log, predicted),
    }
    (OUT / "original_elastic_net_reconstruction_performance.json").write_text(
        json.dumps(original_performance, indent=2), encoding="utf-8")
    paths = [NEW_SCORE_PATH, original.PROTEIN_PATH, original.SCORE_PATH,
             REFERENCE_DIR / "pain_related_proteomic_score_stable_nested_coefficients.csv"]
    pd.DataFrame({"path": [str(p) for p in paths], "md5": [md5(p) for p in paths]}).to_csv(OUT / "input_hashes.csv", index=False)
    manifest = {
        "n": len(data), "n_proteins": len(protein_names), "original_stable_proteins": len(stable),
        "all_new_scores_available": True, "original_and_new_participants_identical": True,
        "original_stable_score_max_difference": maximum_difference,
        "pain_score_correlation_old_new": float(matched.pain_score_original.corr(matched.pain_score_new)),
        "random_state": original.RANDOM_STATE, "outer_folds": 10, "inner_folds": 10,
        "l1_ratios": [0.5, 0.75, 0.9], "n_alphas": 70, "target_offset": 0.1,
        "python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
        "scipy": scipy.__version__, "sklearn": sklearn.__version__,
        "warning": "Exact legacy replication: preprocessing is inside each OUTER fold, but precedes ElasticNetCV inner splitting.",
    }
    (OUT / "preparation_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
