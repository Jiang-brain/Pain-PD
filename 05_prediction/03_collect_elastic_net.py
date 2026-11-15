import json

import numpy as np
import pandas as pd

from common import OUT, REFERENCE_DIR, load_reference


def main():
    original = load_reference("build_stable_nested_elastic_lasso_pfs.py", "original_elastic_collect")
    directories = [OUT / "elastic_net_folds" / f"fold_{fold:02d}" for fold in range(1, 11)]
    if not all((directory / "completed.json").exists() for directory in directories):
        raise FileNotFoundError("Require all 10 completed outer folds")
    coefficients = pd.concat([pd.read_csv(directory / "coefficients.csv") for directory in directories], ignore_index=True)
    metrics = pd.concat([pd.read_csv(directory / "metrics.csv") for directory in directories], ignore_index=True)
    oof = pd.concat([pd.read_csv(directory / "oof_predictions.csv") for directory in directories], ignore_index=True)
    cohort = pd.read_csv(OUT / "matched_pain_scores.csv")
    oof = cohort[["eid"]].merge(oof, on="eid", how="left", validate="one_to_one", sort=False)
    if oof.oof_prediction_log.isna().any() or len(oof) != len(cohort):
        raise ValueError("Incomplete OOF predictions")
    if not metrics.l1_ratio.between(0, 1, inclusive="neither").all():
        raise ValueError("The model must contain both L1 and L2 penalties")
    summary = original.summarize_fold_coefficients(coefficients, "elasticnet")
    stable = summary[summary.stable_all_10_folds].copy()
    if stable.empty:
        raise ValueError("No proteins selected in all 10 folds; do not invent a stable panel")
    x = np.load(OUT / "protein_matrix.npy", mmap_mode="r")
    names = json.loads((OUT / "protein_names.json").read_text())
    positions = {name: i for i, name in enumerate(names)}
    predicted = np.full(len(cohort), float(stable.mean_intercept_raw.mean()))
    for row in stable.itertuples(index=False):
        values = x[:, positions[row.protein]].astype(float)
        values[np.isnan(values)] = row.mean_impute_value
        predicted += values * row.mean_coef_raw_imputed
    scores = pd.DataFrame({
        "eid": cohort.eid, "pain_score": cohort.pain_score_new,
        "pain_score_log_plus_0p1": np.log(cohort.pain_score_new.to_numpy() + .1),
        "weighted_pain_status": cohort.weighted_pain_status_new,
        "stable_elasticnet_pfs_log": predicted,
        "stable_elasticnet_pfs_pain_scale": np.exp(predicted) - .1,
        "stable_elasticnet_pfs_z": (predicted - predicted.mean()) / predicted.std(ddof=0),
    })
    coefficients.to_csv(OUT / "new_elastic_net_fold_coefficients.csv", index=False)
    metrics.to_csv(OUT / "new_elastic_net_fold_metrics.csv", index=False)
    summary.to_csv(OUT / "new_stable_elastic_net_coefficients.csv", index=False)
    scores.to_csv(OUT / "new_stable_elastic_net_pfs.csv", index=False)
    oof.to_csv(OUT / "new_elastic_net_oof_predictions.csv", index=False)
    baseline = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_stable_nested_coefficients.csv")
    baseline_names = set(baseline.loc[baseline.method.eq("elasticnet") & baseline.stable_all_10_folds, "protein"])
    new_names = set(stable.protein)
    membership = pd.DataFrame({"protein": sorted(baseline_names | new_names)})
    membership["original_stable"] = membership.protein.isin(baseline_names)
    membership["new_stable"] = membership.protein.isin(new_names)
    membership.to_csv(OUT / "stable_protein_membership_comparison.csv", index=False)
    performance = {
        "n": len(cohort), "n_stable_original": len(baseline_names), "n_stable_new": len(new_names),
        "retained": len(baseline_names & new_names), "lost": len(baseline_names - new_names),
        "gained": len(new_names - baseline_names),
        "new_elastic_net_oof_target_log": original.evaluate(oof.target_log.to_numpy(), oof.oof_prediction_log.to_numpy()),
        "new_stable_score_target_log_apparent": original.evaluate(scores.pain_score_log_plus_0p1.to_numpy(), predicted),
        "total_convergence_and_other_warnings": int(metrics.warning_count.sum()),
        "note": "Stable mean-coefficient scores are not outer-fold OOF predictions; their performance is reported separately.",
    }
    (OUT / "new_elastic_net_summary.json").write_text(json.dumps(performance, indent=2), encoding="utf-8")
    print(json.dumps(performance, indent=2), flush=True)


if __name__ == "__main__":
    main()
