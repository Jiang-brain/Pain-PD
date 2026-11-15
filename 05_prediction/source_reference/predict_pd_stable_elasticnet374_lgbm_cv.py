"""10-fold LightGBM prediction using 374 stable ElasticNet pain-related proteins.

This reuses the preprocessing, missing-value handling, outcome windows, and
LightGBM settings from predict_pd_intersection_lgbm_cv.py. The only change is
that the 25-protein feature set is replaced by the 374 proteins selected by
nested ElasticNet in all 10 outer folds.
"""

import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BASE_SCRIPT = ROOT / "predict_pd_intersection_lgbm_cv.py"
STABLE_COEF_PATH = ROOT / "pain_related_proteomic_score_stable_nested_coefficients.csv"
PYDEPS = ROOT / "pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import pandas as pd


def load_base_module():
    spec = importlib.util.spec_from_file_location("predict_pd_intersection_lgbm_cv_base", BASE_SCRIPT)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load base LightGBM script: {BASE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_stable_elasticnet_proteins():
    if not STABLE_COEF_PATH.exists():
        raise FileNotFoundError(f"Missing stable ElasticNet coefficient file: {STABLE_COEF_PATH}")
    coef = pd.read_csv(STABLE_COEF_PATH)
    required = {"protein", "method", "stable_all_10_folds"}
    missing = required.difference(coef.columns)
    if missing:
        raise ValueError(
            f"Missing columns in {STABLE_COEF_PATH.name}: {', '.join(sorted(missing))}"
        )
    stable_flag = coef["stable_all_10_folds"].astype(str).str.lower().eq("true")
    proteins = (
        coef.loc[stable_flag & coef["method"].eq("elasticnet"), "protein"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    if len(proteins) != 374:
        raise ValueError(
            "Expected 374 stable ElasticNet proteins, "
            f"but found {len(proteins)} in {STABLE_COEF_PATH}"
        )
    return proteins


def configure_outputs(base):
    base.OUT_FEATURES = ROOT / "pd_prediction_stable_elasticnet374_lgbm_feature_sets.json"
    base.OUT_DATA_FULL = ROOT / "pd_prediction_stable_elasticnet374_lgbm_dataset_full_followup.csv"
    base.OUT_DATA_5Y = ROOT / "pd_prediction_stable_elasticnet374_lgbm_dataset_2to5yr.csv"
    base.OUT_DATA_AFTER5Y = ROOT / "pd_prediction_stable_elasticnet374_lgbm_dataset_after5yr.csv"
    base.OUT_OVERALL_METRICS = ROOT / "pd_prediction_stable_elasticnet374_lgbm_cv_metrics.csv"
    base.OUT_FOLD_METRICS = ROOT / "pd_prediction_stable_elasticnet374_lgbm_fold_metrics.csv"
    base.OUT_PREDICTIONS = ROOT / "pd_prediction_stable_elasticnet374_lgbm_oof_predictions.csv"
    base.OUT_MISSINGNESS = ROOT / "pd_prediction_stable_elasticnet374_lgbm_missingness.csv"


def write_comparison(overall_df):
    comparison_path = ROOT / "pd_prediction_lgbm_25_vs_stable_elasticnet374_comparison.csv"
    original_path = ROOT / "pd_prediction_intersection_lgbm_cv_metrics.csv"
    if not original_path.exists():
        print(f"Comparison skipped; missing original metrics: {original_path}", flush=True)
        return None

    original = pd.read_csv(original_path)
    comparisons = [
        (
            "covariates_only_reproducibility_check",
            "covariates_only",
            "covariates_only",
        ),
        (
            "protein_alone_25_vs_stable_elasticnet374",
            "fdr_same_direction_25",
            "stable_elasticnet_374",
        ),
        (
            "covariates_plus_proteins_25_vs_stable_elasticnet374",
            "covariates_plus_fdr_same_direction_25",
            "covariates_plus_stable_elasticnet_374",
        ),
    ]
    metric_names = ["auc", "auprc", "brier"]
    rows = []
    for comparison_name, old_set, new_set in comparisons:
        old = original[original["feature_set"].eq(old_set)].copy()
        new = overall_df[overall_df["feature_set"].eq(new_set)].copy()
        merged = old.merge(new, on="task", suffixes=("_25", "_stable374"))
        for _, row in merged.iterrows():
            result = {
                "comparison": comparison_name,
                "task": row["task"],
                "feature_set_25": old_set,
                "feature_set_stable374": new_set,
                "n_25": int(row["n_total_25"]),
                "n_stable374": int(row["n_total_stable374"]),
                "events_25": int(row["events_total_25"]),
                "events_stable374": int(row["events_total_stable374"]),
            }
            for metric in metric_names:
                result[f"{metric}_25"] = row[f"{metric}_25"]
                result[f"{metric}_stable374"] = row[f"{metric}_stable374"]
                result[f"{metric}_delta_stable374_minus_25"] = (
                    row[f"{metric}_stable374"] - row[f"{metric}_25"]
                )
            rows.append(result)
    comparison_df = pd.DataFrame(rows).sort_values(["comparison", "task"])
    comparison_df.to_csv(comparison_path, index=False)
    return comparison_path


def main():
    base = load_base_module()
    configure_outputs(base)
    stable_proteins = read_stable_elasticnet_proteins()

    protein_sets = {"stable_elasticnet_374": stable_proteins}
    feature_sets = {
        "covariates_only": base.COVARIATE_FEATURES,
        "stable_elasticnet_374": stable_proteins,
        "covariates_plus_stable_elasticnet_374": [
            *base.COVARIATE_FEATURES,
            *stable_proteins,
        ],
    }

    data, all_proteins, data_path = base.load_prediction_base(protein_sets)
    missing_proteins = sorted(set(all_proteins).difference(data.columns))
    if missing_proteins:
        raise ValueError(
            "Stable ElasticNet protein columns missing from merged data: "
            + ", ".join(missing_proteins)
        )

    tasks = base.make_tasks(data)
    tasks["incident_pd_full_followup"][
        ["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]
    ].to_csv(base.OUT_DATA_FULL, index=False)
    tasks["incident_pd_2to5y"][
        ["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]
    ].to_csv(base.OUT_DATA_5Y, index=False)
    tasks["incident_pd_after5y"][
        ["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]
    ].to_csv(base.OUT_DATA_AFTER5Y, index=False)

    base.OUT_FEATURES.write_text(
        json.dumps(
            {
                "base_script": str(BASE_SCRIPT),
                "data_path": str(data_path),
                "protein_path": str(base.PROTEIN_PATH),
                "stable_coefficient_path": str(STABLE_COEF_PATH),
                "stable_selection_rule": (
                    "method == elasticnet and stable_all_10_folds == TRUE"
                ),
                "n_stable_elasticnet_proteins": len(stable_proteins),
                "excluded_rule": "diagnosis_time <= 365.25*2 days, including prevalent PD",
                "protein_sets": protein_sets,
                "covariate_features": base.COVARIATE_FEATURES,
                "feature_sets": feature_sets,
                "time_windows": {
                    "full_followup": "time > 2 years; any PD after 2 years is an event",
                    "2to5yr": "PD onset > 2 and <= 5 years; censoring at >= 5 years",
                    "after5yr": "time > 5 years; PD onset after 5 years is an event",
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    base.write_missingness(tasks, feature_sets)
    overall_df, _ = base.run_cross_validation(tasks, feature_sets)
    comparison_path = write_comparison(overall_df)

    print("Stable ElasticNet protein source:", STABLE_COEF_PATH, flush=True)
    print("Stable ElasticNet proteins: 374", flush=True)
    print("Excluded rule: diagnosis_time <= 2 years", flush=True)
    for task_name, task_data in tasks.items():
        print(
            f"{task_name}: n={len(task_data)}, events={int(task_data['outcome'].sum())}",
            flush=True,
        )
    print(overall_df.to_string(index=False), flush=True)
    print("Outputs:", flush=True)
    for path in [
        base.OUT_OVERALL_METRICS,
        base.OUT_FOLD_METRICS,
        base.OUT_PREDICTIONS,
        base.OUT_FEATURES,
        base.OUT_MISSINGNESS,
        base.OUT_DATA_FULL,
        base.OUT_DATA_5Y,
        base.OUT_DATA_AFTER5Y,
    ]:
        print(f"  {path}", flush=True)
    if comparison_path is not None:
        print(f"  {comparison_path}", flush=True)


if __name__ == "__main__":
    main()
