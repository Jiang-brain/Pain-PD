import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYDEPS = ROOT / "pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold


DATA_PATH_CANDIDATES = [
    ROOT / "PD_data.csv",
    ROOT.parent / "PA" / "PD_data.csv",
    ROOT.parent / "PA2" / "PD_data.csv",
]
PROTEIN_PATH = ROOT / "protein_new_missing60_filtered.csv"
BONFERRONI_PATH = ROOT / "protein_intersection_Bonferroni05.csv"
FDR_PATH = ROOT / "protein_intersection_FDR05.csv"

OUT_FEATURES = ROOT / "pd_prediction_intersection_lgbm_feature_sets.json"
OUT_DATA_FULL = ROOT / "pd_prediction_intersection_lgbm_dataset_full_followup.csv"
OUT_DATA_5Y = ROOT / "pd_prediction_intersection_lgbm_dataset_2to5yr.csv"
OUT_DATA_AFTER5Y = ROOT / "pd_prediction_intersection_lgbm_dataset_after5yr.csv"
OUT_OVERALL_METRICS = ROOT / "pd_prediction_intersection_lgbm_cv_metrics.csv"
OUT_FOLD_METRICS = ROOT / "pd_prediction_intersection_lgbm_fold_metrics.csv"
OUT_PREDICTIONS = ROOT / "pd_prediction_intersection_lgbm_oof_predictions.csv"
OUT_MISSINGNESS = ROOT / "pd_prediction_intersection_lgbm_missingness.csv"

RANDOM_STATE = 20260523
FOLLOWUP_2Y_DAYS = 365.25 * 2
FOLLOWUP_5Y_DAYS = 365.25 * 5
N_SPLITS = 10

COVARIATE_FEATURES = [
    "age",
    "gender",
    "ethnicity",
    "income",
    "university_education",
    "smoking",
    "alcohol",
    "sedentary",
    "obesity",
    "deprivation",
    "prs",
    *[f"comp_{i}" for i in range(1, 21)],
]


def find_input(candidates):
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("None of these files exists: " + ", ".join(map(str, candidates)))


def make_ethnicity(df):
    ethnicity = df[["21000-0.0", "21000-1.0", "21000-2.0"]].copy()
    primary = ethnicity["21000-0.0"]
    invalid = primary.isin([-1, -3]) | primary.isna()
    fallback = ethnicity.max(axis=1, skipna=True)
    race = primary.mask(invalid, fallback)
    white = race.isin([1, 1001, 1002, 1003])
    out = pd.Series(np.where(white, 1, 0), index=df.index, dtype="float")
    out[race.isin([-1, -3]) | race.isna()] = np.nan
    return out.astype("category")


def make_income(x):
    out = pd.Series(np.nan, index=x.index, dtype="float")
    out[x == -1] = 0
    out[(x >= 1) & (x <= 3)] = 1
    out[x == 4] = 2
    out[x == 5] = 3
    out[x == -3] = 0
    return out.astype("category")


def make_education(x):
    def code_one(value):
        if pd.isna(value):
            return np.nan
        vals = []
        for part in str(value).split("|"):
            try:
                vals.append(float(part))
            except ValueError:
                pass
        if not vals:
            return np.nan
        mapped = []
        for val in vals:
            if val == -3:
                continue
            if val == -7:
                mapped.append(7)
            elif val == 1:
                mapped.append(20)
            elif val == 2:
                mapped.append(13)
            elif val in (3, 4):
                mapped.append(10)
            elif val == 5:
                mapped.append(19)
            elif val == 6:
                mapped.append(15)
        if not mapped:
            return np.nan
        highest = max(mapped)
        if highest == 7:
            return 0
        if highest < 20:
            return 1
        return 2

    return x.map(code_one).astype("category")


def make_binary_current(x):
    out = x.copy().astype("float")
    out[out < 0] = np.nan
    out[out == 2] = 1
    return out.astype("category")


def make_sedentary(x):
    out = pd.Series(np.nan, index=x.index, dtype="float")
    out[((x < 4) & (x >= 0)) | (x == -10)] = 0
    out[x >= 4] = 1
    return out.astype("category")


def make_obesity(bmi):
    out = pd.Series(np.nan, index=bmi.index, dtype="float")
    out[(bmi >= 18.5) & (bmi < 25)] = 0
    out[(bmi >= 25) & (bmi < 30)] = 1
    out[(bmi >= 30) & (bmi < 35)] = 2
    out[bmi >= 35] = 3
    return out.astype("category")


def tertile_code(x):
    q = x.quantile([1 / 3, 2 / 3])
    out = pd.Series(np.nan, index=x.index, dtype="float")
    out[x < q.iloc[0]] = 0
    out[(x >= q.iloc[0]) & (x < q.iloc[1])] = 1
    out[x >= q.iloc[1]] = 2
    return out.astype("category")


def load_pd_outcomes_and_covariates(data_path):
    pc_cols = [f"22009-0.{i}" for i in range(1, 21)]
    cols = [
        "eid",
        "53-0.0",
        "42030-0.0",
        "42032-0.0",
        "190-0.0",
        "191-0.0",
        "40000-0.0",
        "21003-0.0",
        "31-0.0",
        "21000-0.0",
        "21000-1.0",
        "21000-2.0",
        "22189-0.0",
        "738-0.0",
        "6138-0.0",
        "20116-0.0",
        "20117-0.0",
        "1070-0.0",
        "21001-0.0",
        "26260-0.0",
        *pc_cols,
    ]
    df = pd.read_csv(data_path, usecols=cols)
    for col in ["53-0.0", "42030-0.0", "42032-0.0", "191-0.0", "40000-0.0"]:
        df[col] = pd.to_datetime(df[col], errors="coerce")

    diagnosis_date = df[["42030-0.0", "42032-0.0"]].min(axis=1)
    diagnosis_time = (diagnosis_date - df["53-0.0"]).dt.days

    last_diagnosis_date = diagnosis_date.max()
    censoring_time = (last_diagnosis_date - df["53-0.0"]).dt.days
    lost_time = (df["191-0.0"] - df["53-0.0"]).dt.days
    death_time = (df["40000-0.0"] - df["53-0.0"]).dt.days

    time = censoring_time.copy()
    status = diagnosis_time.notna() & (diagnosis_time > 0)
    time[status] = diagnosis_time[status]
    no_event = ~status
    time[no_event & lost_time.notna()] = lost_time[no_event & lost_time.notna()]
    time[no_event & death_time.notna()] = death_time[no_event & death_time.notna()]

    out = pd.DataFrame(
        {
            "eid": df["eid"],
            "status": status.astype(int),
            "time": time,
            "diagnosis_time": diagnosis_time,
            "age": df["21003-0.0"],
            "gender": df["31-0.0"].astype("category"),
            "ethnicity": make_ethnicity(df),
            "deprivation_raw": df["22189-0.0"],
            "income": make_income(df["738-0.0"]),
            "university_education": make_education(df["6138-0.0"]),
            "smoking": make_binary_current(df["20116-0.0"]),
            "alcohol": make_binary_current(df["20117-0.0"]),
            "sedentary": make_sedentary(df["1070-0.0"]),
            "obesity": make_obesity(df["21001-0.0"]),
            "prs_raw": df["26260-0.0"],
        }
    )
    for i, col in enumerate(pc_cols, start=1):
        out[f"comp_{i}"] = df[col]

    out["deprivation"] = tertile_code(out["deprivation_raw"])
    out["prs"] = tertile_code(out["prs_raw"])
    out = out.drop(columns=["deprivation_raw", "prs_raw"])
    out = out[(out["age"] >= 40) & out["time"].notna() & (out["time"] >= 0)]
    return out


def same_direction_proteins(path):
    df = pd.read_csv(path)
    selected = df[np.sign(df["beta_pain_score"]) == np.sign(df["beta_PD_cox"])].copy()
    return selected["protein"].dropna().tolist()


def load_prediction_base(protein_sets):
    data_path = find_input(DATA_PATH_CANDIDATES)
    proteins = sorted({protein for values in protein_sets.values() for protein in values})
    protein_df = pd.read_csv(PROTEIN_PATH, usecols=["eid", "ins_index", *proteins])
    protein_df = (
        protein_df.sort_values(["eid", "ins_index"])
        .drop_duplicates("eid", keep="first")
        .drop(columns=["ins_index"])
    )
    cov = load_pd_outcomes_and_covariates(data_path)
    data = cov.merge(protein_df, on="eid", how="inner")
    early_or_prevalent_pd = data["diagnosis_time"].notna() & (data["diagnosis_time"] <= FOLLOWUP_2Y_DAYS)
    data = data.loc[~early_or_prevalent_pd].reset_index(drop=True)
    return data, proteins, data_path


def make_tasks(data):
    full_followup = data[data["time"] > FOLLOWUP_2Y_DAYS].copy()
    full_followup["outcome"] = (
        (full_followup["status"] == 1)
        & (full_followup["time"] > FOLLOWUP_2Y_DAYS)
    ).astype(int)
    full_followup["task"] = "incident_pd_full_followup"

    within_5y = data[
        ((data["status"] == 1) & (data["time"] > FOLLOWUP_2Y_DAYS) & (data["time"] <= FOLLOWUP_5Y_DAYS))
        | (data["time"] >= FOLLOWUP_5Y_DAYS)
    ].copy()
    within_5y["outcome"] = (
        (within_5y["status"] == 1)
        & (within_5y["time"] > FOLLOWUP_2Y_DAYS)
        & (within_5y["time"] <= FOLLOWUP_5Y_DAYS)
    ).astype(int)
    within_5y["task"] = "incident_pd_2to5y"

    after_5y = data[data["time"] > FOLLOWUP_5Y_DAYS].copy()
    after_5y["outcome"] = ((after_5y["status"] == 1) & (after_5y["time"] > FOLLOWUP_5Y_DAYS)).astype(int)
    after_5y["task"] = "incident_pd_after5y"
    return {
        "incident_pd_2to5y": within_5y,
        "incident_pd_after5y": after_5y,
        "incident_pd_full_followup": full_followup,
    }


def make_model(y_train):
    positives = int(y_train.sum())
    negatives = int(len(y_train) - positives)
    pos_weight = negatives / max(positives, 1)
    return LGBMClassifier(
        n_estimators=500,
        learning_rate=0.02,
        num_leaves=15,
        min_child_samples=30,
        subsample=0.85,
        colsample_bytree=0.85,
        reg_alpha=0.1,
        reg_lambda=1.0,
        objective="binary",
        scale_pos_weight=pos_weight,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        verbose=-1,
    )


def youden_threshold(y_true, prob):
    fpr, tpr, thresholds = roc_curve(y_true, prob)
    idx = int(np.argmax(tpr - fpr))
    return float(thresholds[idx])


def evaluate_predictions(y_true, prob):
    threshold = youden_threshold(y_true, prob)
    pred = (prob >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, pred).ravel()
    return {
        "auc": roc_auc_score(y_true, prob),
        "auprc": average_precision_score(y_true, prob),
        "brier": brier_score_loss(y_true, prob),
        "threshold_youden": threshold,
        "sensitivity": tp / (tp + fn) if (tp + fn) else np.nan,
        "specificity": tn / (tn + fp) if (tn + fp) else np.nan,
        "precision": tp / (tp + fp) if (tp + fp) else np.nan,
        "f1": f1_score(y_true, pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, pred),
        "tp": int(tp),
        "fp": int(fp),
        "tn": int(tn),
        "fn": int(fn),
    }


def run_cross_validation(tasks, feature_sets):
    fold_rows = []
    overall_rows = []
    prediction_frames = []

    for task_name, task_data in tasks.items():
        task_data = task_data.dropna(subset=["outcome"]).reset_index(drop=True)
        for feature_set_name, feature_cols in feature_sets.items():
            y = task_data["outcome"].to_numpy()
            X = task_data[feature_cols]
            oof_prob = np.full(len(task_data), np.nan, dtype=float)
            oof_fold = np.full(len(task_data), -1, dtype=int)
            cv = StratifiedKFold(n_splits=N_SPLITS, shuffle=True, random_state=RANDOM_STATE)

            for fold, (train_idx, test_idx) in enumerate(cv.split(X, y), start=1):
                model = make_model(y[train_idx])
                model.fit(X.iloc[train_idx], y[train_idx])
                prob = model.predict_proba(X.iloc[test_idx])[:, 1]
                oof_prob[test_idx] = prob
                oof_fold[test_idx] = fold
                fold_metrics = evaluate_predictions(y[test_idx], prob)
                fold_rows.append(
                    {
                        "task": task_name,
                        "feature_set": feature_set_name,
                        "model": "lightgbm",
                        "fold": fold,
                        "n_train": int(len(train_idx)),
                        "events_train": int(y[train_idx].sum()),
                        "n_valid": int(len(test_idx)),
                        "events_valid": int(y[test_idx].sum()),
                        **fold_metrics,
                    }
                )

            overall_metrics = evaluate_predictions(y, oof_prob)
            overall_rows.append(
                {
                    "task": task_name,
                    "feature_set": feature_set_name,
                    "model": "lightgbm",
                    "n_total": int(len(task_data)),
                    "events_total": int(y.sum()),
                    "n_splits": N_SPLITS,
                    **overall_metrics,
                }
            )
            prediction_frames.append(
                pd.DataFrame(
                    {
                        "eid": task_data["eid"].to_numpy(),
                        "task": task_name,
                        "feature_set": feature_set_name,
                        "model": "lightgbm",
                        "fold": oof_fold,
                        "outcome": y,
                        "predicted_risk": oof_prob,
                    }
                )
            )

    overall_df = pd.DataFrame(overall_rows).sort_values(["task", "feature_set"])
    fold_df = pd.DataFrame(fold_rows).sort_values(["task", "feature_set", "fold"])
    predictions_df = pd.concat(prediction_frames, ignore_index=True)
    overall_df.to_csv(OUT_OVERALL_METRICS, index=False)
    fold_df.to_csv(OUT_FOLD_METRICS, index=False)
    predictions_df.to_csv(OUT_PREDICTIONS, index=False)
    return overall_df, fold_df


def write_missingness(tasks, feature_sets):
    rows = []
    for task_name, task_data in tasks.items():
        for feature_set_name, feature_cols in feature_sets.items():
            for feature in feature_cols:
                rows.append(
                    {
                        "task": task_name,
                        "feature_set": feature_set_name,
                        "feature": feature,
                        "n": int(len(task_data)),
                        "missing": int(task_data[feature].isna().sum()),
                        "missing_fraction": float(task_data[feature].isna().mean()),
                    }
                )
    pd.DataFrame(rows).to_csv(OUT_MISSINGNESS, index=False)


def main():
    protein_sets = {
        "fdr_same_direction_25": same_direction_proteins(FDR_PATH),
    }
    if len(protein_sets["fdr_same_direction_25"]) != 25:
        raise ValueError(f"Expected 25 FDR same-direction proteins, got {len(protein_sets['fdr_same_direction_25'])}")

    fdr25_proteins = protein_sets["fdr_same_direction_25"]
    feature_sets = {
        "covariates_only": COVARIATE_FEATURES,
        "fdr_same_direction_25": fdr25_proteins,
        "covariates_plus_fdr_same_direction_25": [
            *COVARIATE_FEATURES,
            *fdr25_proteins,
        ],
    }

    data, all_proteins, data_path = load_prediction_base(protein_sets)
    missing_proteins = sorted(set(all_proteins) - set(data.columns))
    if missing_proteins:
        raise ValueError("Protein columns missing from merged data: " + ", ".join(missing_proteins))

    tasks = make_tasks(data)
    tasks["incident_pd_full_followup"][["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]].to_csv(
        OUT_DATA_FULL, index=False
    )
    tasks["incident_pd_2to5y"][["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]].to_csv(
        OUT_DATA_5Y, index=False
    )
    tasks["incident_pd_after5y"][["eid", "status", "time", "diagnosis_time", "outcome", *all_proteins]].to_csv(
        OUT_DATA_AFTER5Y, index=False
    )
    OUT_FEATURES.write_text(
        json.dumps(
            {
                "data_path": str(data_path),
                "protein_path": str(PROTEIN_PATH),
                "excluded_rule": "diagnosis_time <= 365.25*2 days, including prevalent PD",
                "protein_sets": protein_sets,
                "covariate_features": COVARIATE_FEATURES,
                "feature_sets": feature_sets,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    write_missingness(tasks, feature_sets)
    overall_df, fold_df = run_cross_validation(tasks, feature_sets)

    print("Feature sets:", flush=True)
    for name, features in feature_sets.items():
        print(f"  {name} ({len(features)} features)", flush=True)
        if name != "covariates_only":
            proteins = [feature for feature in features if feature not in COVARIATE_FEATURES]
            print(f"    proteins ({len(proteins)}): {', '.join(proteins)}", flush=True)
    for task_name, task_data in tasks.items():
        print(f"{task_name}: n={len(task_data)}, events={int(task_data['outcome'].sum())}", flush=True)
    print(overall_df.to_string(index=False), flush=True)
    print("Outputs:", flush=True)
    for path in [OUT_OVERALL_METRICS, OUT_FOLD_METRICS, OUT_PREDICTIONS, OUT_FEATURES, OUT_MISSINGNESS]:
        print(f"  {path}", flush=True)


if __name__ == "__main__":
    main()
