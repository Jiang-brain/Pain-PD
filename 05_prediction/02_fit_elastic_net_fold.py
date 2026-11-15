"""Refit one complete original 10-inner-fold Elastic Net outer fold."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
import json
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import KFold
from threadpoolctl import threadpool_limits

from common import OUT, load_reference


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fold", type=int, choices=range(1, 11))
    args = parser.parse_args()
    fold = args.fold
    destination = OUT / "elastic_net_folds" / f"fold_{fold:02d}"
    destination.mkdir(parents=True, exist_ok=True)
    if (destination / "completed.json").exists():
        raise FileExistsError(f"Completed fold already exists; inspect provenance before reuse: {destination}")
    original = load_reference("build_stable_nested_elastic_lasso_pfs.py", "original_elastic_net_fold")
    x = np.load(OUT / "protein_matrix.npy", mmap_mode="r")
    target = np.load(OUT / "target_new_log.npy")
    folds = pd.read_csv(OUT / "outer_fold_assignments.csv")
    names = json.loads((OUT / "protein_names.json").read_text())
    train = np.flatnonzero(folds.outer_fold.to_numpy() != fold)
    test = np.flatnonzero(folds.outer_fold.to_numpy() == fold)
    inner = KFold(n_splits=10, shuffle=True, random_state=original.RANDOM_STATE + fold)
    pipe = original.make_elastic_pipeline(inner, fold)
    # Resource setting only: grid, folds, penalties, seeds, and solver remain
    # unchanged. Limit simultaneous inner-fit copies on the 16-GB workstation.
    inner_jobs = int(os.environ.get("PD_PAIN_EN_INNER_JOBS", "2"))
    pipe.named_steps["model"].set_params(n_jobs=inner_jobs)
    print(f"Outer fold {fold}: train={len(train)}, test={len(test)}, full inner 10-fold ElasticNetCV", flush=True)
    started = time.monotonic()
    captured = []
    with warnings.catch_warnings(record=True) as warning_list, threadpool_limits(limits=1):
        warnings.simplefilter("always", ConvergenceWarning)
        pipe.fit(x[train], target[train])
        prediction = pipe.predict(x[test])
        captured = [{"category": w.category.__name__, "message": str(w.message)} for w in warning_list]
    model = pipe.named_steps["model"]
    imputer = pipe.named_steps["imputer"]
    scaler = pipe.named_steps["scaler"]
    intercept, raw = original.raw_scale_intercept(model.intercept_, model.coef_, scaler.mean_, scaler.scale_)
    metrics = {"fold": fold, "n_train": len(train), "n_test": len(test), "alpha": float(model.alpha_),
               "l1_ratio": float(model.l1_ratio_), "n_selected_proteins": int((np.abs(model.coef_) > 1e-12).sum()),
               "elapsed_seconds": time.monotonic() - started, "warning_count": len(captured),
               "inner_execution_jobs": inner_jobs,
               **original.evaluate(target[test], prediction)}
    coefficient = pd.DataFrame({
        "fold": fold, "protein": names,
        "coef_standardized_protein_to_log_pain": model.coef_,
        "coef_raw_imputed_protein_to_log_pain": raw, "selected": np.abs(model.coef_) > 1e-12,
        "protein_median_impute": imputer.statistics_, "protein_mean_after_impute": scaler.mean_,
        "protein_sd_after_impute": scaler.scale_, "model_intercept_standardized_space": float(model.intercept_),
        "model_intercept_raw_imputed_space": intercept, "alpha": float(model.alpha_), "l1_ratio": float(model.l1_ratio_),
    })
    coefficient.to_csv(destination / "coefficients.csv", index=False)
    pd.DataFrame({"eid": folds.eid.iloc[test], "outer_fold": fold, "target_log": target[test],
                  "oof_prediction_log": prediction}).to_csv(destination / "oof_predictions.csv", index=False)
    pd.DataFrame([metrics]).to_csv(destination / "metrics.csv", index=False)
    (destination / "warnings.json").write_text(json.dumps(captured, indent=2), encoding="utf-8")
    (destination / "completed.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()
