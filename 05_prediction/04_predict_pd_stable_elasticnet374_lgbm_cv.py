"""LightGBM prediction of incident PD using the stable Elastic Net protein panel.

This script implements the manuscript's exploratory classification-based
prediction analysis for three windows:

- full follow-up after excluding PD events within 2 years after baseline
- PD diagnosed from >2 to <=5 years
- PD diagnosed from >5 to <=10 years

Feature sets:

- conventional covariates only
- stable Elastic Net protein panel only
- conventional covariates plus the stable Elastic Net protein panel

The final PPS model is Elastic Net, not LASSO, and this script does not use PCA.
"""

import json
import os
import sys
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("PD_PAIN_DATA_DIR", CODE_ROOT / "data"))
RESULTS_DIR = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE_ROOT / "results"))
PYDEPS = CODE_ROOT / "pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

try:
    from lightgbm import LGBMClassifier
except ImportError as exc:
    raise ImportError(
        "Missing Python package 'lightgbm'. Install it before running this script."
    ) from exc


RANDOM_STATE = 20260721
N_SPLITS = 10
TWO_YEARS = 365.25 * 2
FIVE_YEARS = 365.25 * 5
TEN_YEARS = 365.25 * 10

ANALYSIS_DATA_PATH = RESULTS_DIR / "ukb_main_analysis_dataset.csv"
PROTEIN_PATH = DATA_DIR / "protein_new_missing60_filtered.csv"
STABLE_COEF_PATH = RESULTS_DIR / "pain_related_proteomic_score_stable_elastic_net_coefficients.csv"

OUT_FEATURES = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_feature_sets.json"
OUT_DATA_FULL = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_dataset_full_followup.csv"
OUT_DATA_5Y = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_dataset_2to5yr.csv"
OUT_DATA_10Y = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_dataset_5to10yr.csv"
OUT_OVERALL_METRICS = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_cv_metrics.csv"
OUT_FOLD_METRICS = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_fold_metrics.csv"
OUT_PREDICTIONS = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_oof_predictions.csv"
OUT_MISSINGNESS = RESULTS_DIR / "pd_prediction_stable_elasticnet374_lgbm_missingness.csv"

COVARIATE_FEATURES = [
    "age",
    "gender",
    "ethnicity",
    "deprivation",
    "income",
    "university_education",
    "smoking",
    "alcohol",
    "sedentary",
    "obesity",
    "prs",
    *[f"comp_{i}" for i in range(1, 11)],
]


def read_stable_elastic_net_proteins():
    if not STABLE_COEF_PATH.exists():
        raise FileNotFoundError(
            f"Missing stable Elastic Net coefficient file: {STABLE_COEF_PATH}. "
            "Run 05_prediction/02_build_stable_nested_elastic_net_pfs.py first."
        )
    coef = pd.read_csv(STABLE_COEF_PATH)
    required = {"protein", "method", "stable_all_10_folds"}
    missing = required.difference(coef.columns)
    if missing:
        raise ValueError(f"Missing columns in {STABLE_COEF_PATH.name}: {sorted(missing)}")
    stable_flag = coef["stable_all_10_folds"].astype(str).str.lower().eq("true")
    method = coef["method"].astype(str).str.lower().str.replace("-", "_", regex=False)
    proteins = (
        coef.loc[stable_flag & method.eq("elastic_net"), "protein"]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .sort_values()
        .tolist()
    )
    if len(proteins) != 374:
        raise ValueError(
            f"Expected 374 stable Elastic Net proteins, found {len(proteins)} in {STABLE_COEF_PATH}."
        )
    return proteins


def load_prediction_data(stable_proteins):
    if not ANALYSIS_DATA_PATH.exists():
        raise FileNotFoundError(
            f"Missing analysis dataset: {ANALYSIS_DATA_PATH}. "
            "Run 01_data_preparation/04_prepare_ukb_analysis_dataset.R first."
        )
    if not PROTEIN_PATH.exists():
        raise FileNotFoundError(f"Missing protein matrix: {PROTEIN_PATH}")
    data = pd.read_csv(ANALYSIS_DATA_PATH)
    proteins = pd.read_csv(PROTEIN_PATH, usecols=["eid", "ins_index", *stable_proteins])
    proteins = (
        proteins.sort_values(["eid", "ins_index"])
        .drop_duplicates("eid", keep="first")
        .drop(columns=["ins_index"])
        .reset_index(drop=True)
    )
    data = data.merge(proteins, on="eid", how="inner")
    missing_covariates = sorted(set(COVARIATE_FEATURES).difference(data.columns))
    if missing_covariates:
        raise ValueError(f"Missing covariate columns: {missing_covariates}")
    return data


def make_tasks(data):
    full = data[data["time"] > TWO_YEARS].copy()
    full["outcome"] = ((full["status"] == 1) & (full["time"] > TWO_YEARS)).astype(int)
    full["task"] = "incident_pd_full_followup"

    two_to_five = data[((data["status"] == 1) & (data["time"] > TWO_YEARS) & (data["time"] <= FIVE_YEARS)) | (data["time"] >= FIVE_YEARS)].copy()
    two_to_five["outcome"] = ((two_to_five["status"] == 1) & (two_to_five["time"] > TWO_YEARS) & (two_to_five["time"] <= FIVE_YEARS)).astype(int)
    two_to_five["task"] = "incident_pd_2to5y"

    five_to_ten = data[((data["status"] == 1) & (data["time"] > FIVE_YEARS) & (data["time"] <= TEN_YEARS)) | (data["time"] >= TEN_YEARS)].copy()
    five_to_ten["outcome"] = ((five_to_ten["status"] == 1) & (five_to_ten["time"] > FIVE_YEARS) & (five_to_ten["time"] <= TEN_YEARS)).astype(int)
    five_to_ten["task"] = "incident_pd_5to10y"

    return {
        "incident_pd_full_followup": full,
        "incident_pd_2to5y": two_to_five,
        "incident_pd_5to10y": five_to_ten,
    }


def make_pipeline(features):
    numeric_features = [f for f in features if f == "age" or f.startswith("comp_") or f not in COVARIATE_FEATURES]
    categorical_features = [f for f in features if f in COVARIATE_FEATURES and f not in numeric_features]
    transformers = []
    if numeric_features:
        transformers.append(
            (
                "numeric",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                numeric_features,
            )
        )
    if categorical_features:
        transformers.append(
            (
                "categorical",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                categorical_features,
            )
        )
    preprocessor = ColumnTransformer(transformers=transformers)
    model = LGBMClassifier(
        objective="binary",
        n_estimators=300,
        learning_rate=0.03,
        num_leaves=31,
        subsample=0.8,
        colsample_bytree=0.8,
        min_child_samples=20,
        random_state=RANDOM_STATE,
        n_jobs=4,
        verbose=-1,
    )
    return Pipeline(steps=[("preprocess", preprocessor), ("model", model)])


def compute_metrics(y_true, pred):
    return {
        "auc": float(roc_auc_score(y_true, pred)),
        "auprc": float(average_precision_score(y_true, pred)),
        "brier": float(brier_score_loss(y_true, pred)),
    }


def write_missingness(tasks, feature_sets):
    rows = []
    for task_name, task_data in tasks.items():
        for feature_set, features in feature_sets.items():
            for feature in features:
                rows.append(
                    {
                        "task": task_name,
                        "feature_set": feature_set,
                        "feature": feature,
                        "missing_n": int(task_data[feature].isna().sum()),
                        "missing_fraction": float(task_data[feature].isna().mean()),
                    }
                )
    pd.DataFrame(rows).to_csv(OUT_MISSINGNESS, index=False)


def run_cross_validation(tasks, feature_sets):
    overall_rows = []
    fold_rows = []
    prediction_rows = []
    for task_name, task_data in tasks.items():
        y = task_data["outcome"].astype(int).to_numpy()
        if y.sum() < N_SPLITS or (len(y) - y.sum()) < N_SPLITS:
            raise ValueError(f"Task {task_name} has too few events/non-events for {N_SPLITS}-fold CV.")
        cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)
        for feature_set, features in feature_sets.items():
            oof = np.full(len(task_data), np.nan, dtype=float)
            for fold, (train_idx, test_idx) in enumerate(cv.split(task_data[features], y), start=1):
                pipe = make_pipeline(features)
                pipe.fit(task_data.iloc[train_idx][features], y[train_idx])
                pred = pipe.predict_proba(task_data.iloc[test_idx][features])[:, 1]
                oof[test_idx] = pred
                fold_metric = compute_metrics(y[test_idx], pred)
                fold_rows.append(
                    {
                        "task": task_name,
                        "feature_set": feature_set,
                        "fold": fold,
                        "n_test": int(len(test_idx)),
                        "events_test": int(y[test_idx].sum()),
                        **fold_metric,
                    }
                )
                prediction_rows.append(
                    pd.DataFrame(
                        {
                            "task": task_name,
                            "feature_set": feature_set,
                            "fold": fold,
                            "eid": task_data.iloc[test_idx]["eid"].to_numpy(),
                            "outcome": y[test_idx],
                            "predicted_probability": pred,
                        }
                    )
                )
            overall_rows.append(
                {
                    "task": task_name,
                    "feature_set": feature_set,
                    "n_total": int(len(y)),
                    "events_total": int(y.sum()),
                    **compute_metrics(y, oof),
                }
            )
    overall = pd.DataFrame(overall_rows).sort_values(["task", "feature_set"])
    folds = pd.DataFrame(fold_rows).sort_values(["task", "feature_set", "fold"])
    predictions = pd.concat(prediction_rows, ignore_index=True)
    overall.to_csv(OUT_OVERALL_METRICS, index=False)
    folds.to_csv(OUT_FOLD_METRICS, index=False)
    predictions.to_csv(OUT_PREDICTIONS, index=False)
    return overall, folds


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stable_proteins = read_stable_elastic_net_proteins()
    data = load_prediction_data(stable_proteins)
    tasks = make_tasks(data)
    feature_sets = {
        "covariates_only": COVARIATE_FEATURES,
        "stable_elastic_net_374": stable_proteins,
        "covariates_plus_stable_elastic_net_374": [*COVARIATE_FEATURES, *stable_proteins],
    }

    tasks["incident_pd_full_followup"][["eid", "status", "time", "outcome", *stable_proteins]].to_csv(OUT_DATA_FULL, index=False)
    tasks["incident_pd_2to5y"][["eid", "status", "time", "outcome", *stable_proteins]].to_csv(OUT_DATA_5Y, index=False)
    tasks["incident_pd_5to10y"][["eid", "status", "time", "outcome", *stable_proteins]].to_csv(OUT_DATA_10Y, index=False)

    OUT_FEATURES.write_text(
        json.dumps(
            {
                "stable_coefficient_path": str(STABLE_COEF_PATH),
                "stable_selection_rule": "method == elastic_net and stable_all_10_folds == TRUE",
                "n_stable_elastic_net_proteins": len(stable_proteins),
                "random_state": RANDOM_STATE,
                "n_splits": N_SPLITS,
                "uses_pca": False,
                "model": "LightGBM binary classifier with stratified 10-fold cross-validation",
                "time_windows": {
                    "full_followup": "time > 2 years; any PD after 2 years is an event",
                    "2to5yr": "PD onset >2 and <=5 years vs no PD through at least 5 years",
                    "5to10yr": "PD onset >5 and <=10 years vs no PD through at least 10 years",
                },
                "feature_sets": feature_sets,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    write_missingness(tasks, feature_sets)
    overall, _ = run_cross_validation(tasks, feature_sets)
    print(overall.to_string(index=False), flush=True)
    print("Outputs:", flush=True)
    for path in [
        OUT_OVERALL_METRICS,
        OUT_FOLD_METRICS,
        OUT_PREDICTIONS,
        OUT_FEATURES,
        OUT_MISSINGNESS,
        OUT_DATA_FULL,
        OUT_DATA_5Y,
        OUT_DATA_10Y,
    ]:
        print(f"  {path}", flush=True)


if __name__ == "__main__":
    main()
