"""Original LightGBM settings and folds, with checkpointed controlled panels."""
import argparse
import json
import os
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from common import DATA_DIR, OUT, REFERENCE_DIR, SENSITIVITY_DIR, load_reference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=["baseline", "new"])
    parser.add_argument("--task", choices=["incident_pd_2to5y", "incident_pd_after5y", "incident_pd_full_followup"])
    parser.add_argument("--collect-only", action="store_true")
    args = parser.parse_args()
    base = load_reference("predict_pd_intersection_lgbm_cv.py", "original_lightgbm")
    base.PROTEIN_PATH = REFERENCE_DIR / "protein_new_missing60_filtered.csv"
    base.DATA_PATH_CANDIDATES = [DATA_DIR / "PD_data.csv"]
    phase_directory = OUT / f"lightgbm_{args.phase}"
    phase_directory.mkdir(exist_ok=True)
    if args.collect_only:
        task_names = ["incident_pd_2to5y", "incident_pd_after5y", "incident_pd_full_followup"]
        for filename in ["cv_metrics.csv", "fold_metrics.csv", "oof_predictions.csv"]:
            parts = [pd.read_csv(phase_directory / task / filename) for task in task_names]
            pd.concat(parts, ignore_index=True).to_csv(phase_directory / filename, index=False)
        print(pd.read_csv(phase_directory / "cv_metrics.csv").to_string(index=False), flush=True)
        return
    destination = phase_directory / args.task if args.task else phase_directory
    destination.mkdir(exist_ok=True)
    if args.phase == "baseline":
        coefficients = pd.read_csv(REFERENCE_DIR / "pain_related_proteomic_score_stable_nested_coefficients.csv")
        panel = sorted(coefficients.loc[coefficients.method.eq("elasticnet") & coefficients.stable_all_10_folds, "protein"].unique())
        if len(panel) != 374:
            raise ValueError("Original panel must contain 374 proteins")
        feature_sets = {
            "covariates_only": base.COVARIATE_FEATURES,
            "stable_elasticnet_374": panel,
            "covariates_plus_stable_elasticnet_374": [*base.COVARIATE_FEATURES, *panel],
        }
        protein_sets = {"original_stable": panel}
    else:
        coefficients = pd.read_csv(OUT / "new_stable_elastic_net_coefficients.csv")
        panel = sorted(coefficients.loc[coefficients.method.eq("elasticnet") & coefficients.stable_all_10_folds, "protein"].unique())
        shared = sorted(pd.read_csv(SENSITIVITY_DIR / "shared_proteins_fdr05_concordant.csv").protein.unique())
        feature_sets = {
            "covariates_only": base.COVARIATE_FEATURES,
            "new_stable_panel": panel,
            "covariates_plus_new_stable_panel": [*base.COVARIATE_FEATURES, *panel],
            "new_shared24": shared,
            "covariates_plus_new_shared24": [*base.COVARIATE_FEATURES, *shared],
        }
        protein_sets = {"new_stable": panel, "new_shared": shared}
    print(f"Loading original LightGBM prediction cohort; phase={args.phase}, stable proteins={len(panel)}", flush=True)
    data, all_proteins, _ = base.load_prediction_base(protein_sets)
    tasks = base.make_tasks(data)
    if args.task:
        tasks = {args.task: tasks[args.task]}
    feature_metadata = {
        "phase": args.phase, "n_stable": len(panel), "protein_sets": protein_sets,
        "feature_sets": feature_sets, "random_state": base.RANDOM_STATE, "n_splits": base.N_SPLITS,
        "time_windows": "Original full follow-up, >2 to <=5 years, and >5 years (not capped at 10 years).",
        "execution_threads": int(os.environ.get("PD_PAIN_LGBM_THREADS", "4")),
        "warning": "Protein panels are selected before PD cross-validation, as in the original analysis; selection is not nested within PD folds.",
    }
    (destination / "feature_sets.json").write_text(json.dumps(feature_metadata, indent=2), encoding="utf-8")
    fold_rows, overall_rows, prediction_frames = [], [], []
    for task_name, task_data in tasks.items():
        task_data = task_data.dropna(subset=["outcome"]).reset_index(drop=True)
        task_data[["eid", "status", "time", "diagnosis_time", "outcome"]].to_csv(destination / f"cohort_{task_name}.csv", index=False)
        y = task_data.outcome.to_numpy()
        for name, features in feature_sets.items():
            oof = np.full(len(task_data), np.nan)
            fold_ids = np.full(len(task_data), -1, dtype=int)
            cv = StratifiedKFold(n_splits=base.N_SPLITS, shuffle=True, random_state=base.RANDOM_STATE)
            started = time.monotonic()
            for fold, (train, test) in enumerate(cv.split(task_data[features], y), 1):
                checkpoint = destination / f"{task_name}_{name}_fold{fold:02d}.csv"
                if args.phase == "new" and name == "covariates_only" and not checkpoint.exists():
                    baseline_checkpoint = OUT / "lightgbm_baseline" / checkpoint.name
                    if baseline_checkpoint.exists():
                        # Same covariates, phenotype, participant order, and CV
                        # split: reuse the validated control instead of refitting.
                        control = pd.read_csv(baseline_checkpoint)
                        if np.array_equal(control.eid.to_numpy(), task_data.eid.iloc[test].to_numpy()) and np.array_equal(control.outcome.to_numpy(), y[test]):
                            control.to_csv(checkpoint, index=False)
                if checkpoint.exists():
                    saved = pd.read_csv(checkpoint)
                    if not np.array_equal(saved.eid.to_numpy(), task_data.eid.iloc[test].to_numpy()) or not np.array_equal(saved.outcome.to_numpy(), y[test]):
                        raise ValueError(f"Checkpoint cohort mismatch: {checkpoint}")
                    probability = saved.predicted_risk.to_numpy()
                else:
                    model = base.make_model(y[train])
                    model.set_params(n_jobs=int(os.environ.get("PD_PAIN_LGBM_THREADS", "4")))
                    model.fit(task_data.iloc[train][features], y[train])
                    probability = model.predict_proba(task_data.iloc[test][features])[:, 1]
                    saved = pd.DataFrame({"eid": task_data.eid.iloc[test].to_numpy(), "outcome": y[test], "predicted_risk": probability})
                    saved.to_csv(checkpoint, index=False)
                oof[test] = probability
                fold_ids[test] = fold
                fold_rows.append({"task": task_name, "feature_set": name, "fold": fold,
                                  "n_train": len(train), "n_valid": len(test),
                                  **base.evaluate_predictions(y[test], probability)})
                print(f"{args.phase}: {task_name}, {name}, fold {fold}/10 complete", flush=True)
            overall_rows.append({"task": task_name, "feature_set": name, "model": "lightgbm",
                                 "n_total": len(y), "events_total": int(y.sum()), "n_splits": base.N_SPLITS,
                                 "elapsed_seconds": time.monotonic() - started, **base.evaluate_predictions(y, oof)})
            prediction_frames.append(pd.DataFrame({"eid": task_data.eid, "task": task_name, "feature_set": name,
                                                   "fold": fold_ids, "outcome": y, "predicted_risk": oof}))
    overall = pd.DataFrame(overall_rows)
    overall.to_csv(destination / "cv_metrics.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(destination / "fold_metrics.csv", index=False)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    predictions.to_csv(destination / "oof_predictions.csv", index=False)
    if args.phase == "baseline":
        reference = pd.read_csv(REFERENCE_DIR / "pd_prediction_stable_elasticnet374_lgbm_cv_metrics.csv")
        check = overall.merge(reference, on=["task", "feature_set"], suffixes=("_rerun", "_saved"), validate="one_to_one")
        check["passed"] = (check.n_total_rerun.eq(check.n_total_saved) & check.events_total_rerun.eq(check.events_total_saved))
        for metric in ["auc", "auprc", "brier"]:
            check[f"{metric}_absolute_difference"] = abs(check[f"{metric}_rerun"] - check[f"{metric}_saved"])
            check["passed"] &= check[f"{metric}_absolute_difference"] < 1e-8
        check.to_csv(destination / "reproduction_checks.csv", index=False)
        if len(check) != 9 or not check.passed.all():
            raise ValueError("Original LightGBM baseline is not numerically reproduced; inspect saved checks")
    print(overall.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
