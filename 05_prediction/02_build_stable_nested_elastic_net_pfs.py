import json
import os
import sys
import warnings
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(os.environ.get("PD_PAIN_DATA_DIR", CODE_ROOT / "data"))
RESULTS_DIR = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE_ROOT / "results"))
PYDEPS = CODE_ROOT / "pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


PROTEIN_PATH = DATA_DIR / "protein_new_missing60_filtered.csv"
SCORE_PATH = RESULTS_DIR / "composite_pain_scores.csv"

ELASTIC_NET_FOLD_COEF = RESULTS_DIR / "pain_related_proteomic_score_nested_elastic_net_fold_coefficients.csv"
ELASTIC_NET_COEF_SUMMARY = RESULTS_DIR / "pain_related_proteomic_score_nested_elastic_net_coefficient_summary.csv"
ELASTIC_NET_FOLD_METRICS = RESULTS_DIR / "pain_related_proteomic_score_nested_elastic_net_fold_metrics.csv"
OUT_STABLE_COEF = RESULTS_DIR / "pain_related_proteomic_score_stable_elastic_net_coefficients.csv"
OUT_STABLE_SCORE = RESULTS_DIR / "pain_related_proteomic_score_stable_elastic_net_pfs.csv"
OUT_SUMMARY = RESULTS_DIR / "pain_related_proteomic_score_stable_elastic_net_summary.json"

RANDOM_STATE = 20260525
N_OUTER_FOLDS = 10
N_INNER_FOLDS = 10
TARGET_OFFSET = 0.1
L1_RATIOS = [0.5, 0.75, 0.9]


def read_inputs(protein_subset=None):
    usecols = ["eid", "ins_index"]
    if protein_subset is None:
        header = pd.read_csv(PROTEIN_PATH, nrows=0).columns.tolist()
        protein_cols = [col for col in header if col not in {"eid", "ins_index"}]
    else:
        protein_cols = list(protein_subset)
    dtype = {col: "float32" for col in protein_cols}
    proteins = pd.read_csv(PROTEIN_PATH, usecols=usecols + protein_cols, dtype=dtype)
    proteins = (
        proteins.sort_values(["eid", "ins_index"])
        .drop_duplicates("eid", keep="first")
        .drop(columns=["ins_index"])
        .reset_index(drop=True)
    )
    scores = pd.read_csv(SCORE_PATH, usecols=["eid", "pain_score", "weighted_pain_status"])
    data = scores.merge(proteins, on="eid", how="inner")
    data = data.dropna(subset=["pain_score"]).reset_index(drop=True)
    return data, protein_cols


def evaluate(y_true, y_pred):
    pearson = pearsonr(y_true, y_pred)
    spearman = spearmanr(y_true, y_pred)
    return {
        "r2": float(r2_score(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_r": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def raw_scale_intercept(model_intercept, coef_standardized, scaler_mean, scaler_scale):
    coef_raw = coef_standardized / scaler_scale
    intercept_raw = model_intercept - float(np.sum(coef_standardized * scaler_mean / scaler_scale))
    return float(intercept_raw), coef_raw


def make_elastic_net_pipeline(inner_cv, fold):
    model = ElasticNetCV(
        l1_ratio=L1_RATIOS,
        eps=1e-3,
        alphas=70,
        cv=inner_cv,
        max_iter=8000,
        tol=1e-3,
        n_jobs=4,
        random_state=RANDOM_STATE + fold,
        selection="random",
    )
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def summarize_fold_coefficients(coef_long):
    out = (
        coef_long.groupby("protein", as_index=False)
        .agg(
            selection_count=("selected", "sum"),
            mean_coef_standardized=("coef_standardized_protein_to_log_pain", "mean"),
            sd_coef_standardized=("coef_standardized_protein_to_log_pain", "std"),
            mean_abs_coef_standardized=("coef_standardized_protein_to_log_pain", lambda s: float(np.mean(np.abs(s)))),
            mean_coef_raw_imputed=("coef_raw_imputed_protein_to_log_pain", "mean"),
            sd_coef_raw_imputed=("coef_raw_imputed_protein_to_log_pain", "std"),
            mean_impute_value=("protein_median_impute", "mean"),
            mean_intercept_raw=("model_intercept_raw_imputed_space", "mean"),
            mean_alpha=("alpha", "mean"),
            mean_l1_ratio=("l1_ratio", "mean"),
        )
    )
    out["method"] = "elastic_net"
    out["selection_frequency"] = out["selection_count"] / N_OUTER_FOLDS
    out["stable_all_10_folds"] = out["selection_count"] == N_OUTER_FOLDS
    out = out.sort_values(
        ["stable_all_10_folds", "selection_count", "mean_abs_coef_standardized"],
        ascending=[False, False, False],
    )
    return out


def run_nested_elastic_net():
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print("reading full protein matrix for nested Elastic Net...", flush=True)
    data, protein_cols = read_inputs()
    x = data[protein_cols]
    y_log = np.log(data["pain_score"].to_numpy(dtype=np.float64) + TARGET_OFFSET)
    outer_cv = KFold(n_splits=N_OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    coef_frames = []
    fold_metric_rows = []
    oof_pred = np.full(len(data), np.nan, dtype=np.float64)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(x), start=1):
        print(f"Elastic Net outer fold {fold}/{N_OUTER_FOLDS}: fitting inner {N_INNER_FOLDS}-fold model...", flush=True)
        inner_cv = KFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE + fold)
        pipe = make_elastic_net_pipeline(inner_cv, fold)
        pipe.fit(x.iloc[train_idx], y_log[train_idx])
        model = pipe.named_steps["model"]
        imputer = pipe.named_steps["imputer"]
        scaler = pipe.named_steps["scaler"]
        pred = pipe.predict(x.iloc[test_idx])
        oof_pred[test_idx] = pred
        selected = np.abs(model.coef_) > 1e-12
        fold_metrics = evaluate(y_log[test_idx], pred)
        fold_metric_rows.append(
            {
                "fold": fold,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "alpha": float(model.alpha_),
                "l1_ratio": float(model.l1_ratio_),
                "n_selected_proteins": int(selected.sum()),
                **fold_metrics,
            }
        )
        intercept_raw, coef_raw = raw_scale_intercept(model.intercept_, model.coef_, scaler.mean_, scaler.scale_)
        coef_frames.append(
            pd.DataFrame(
                {
                    "fold": fold,
                    "protein": protein_cols,
                    "coef_standardized_protein_to_log_pain": model.coef_,
                    "coef_raw_imputed_protein_to_log_pain": coef_raw,
                    "selected": selected,
                    "protein_median_impute": imputer.statistics_,
                    "protein_mean_after_impute": scaler.mean_,
                    "protein_sd_after_impute": scaler.scale_,
                    "model_intercept_standardized_space": float(model.intercept_),
                    "model_intercept_raw_imputed_space": intercept_raw,
                    "alpha": float(model.alpha_),
                    "l1_ratio": float(model.l1_ratio_),
                }
            )
        )
        print(
            f"Elastic Net fold {fold}: alpha={model.alpha_:.6g}, l1_ratio={model.l1_ratio_}, "
            f"selected={int(selected.sum())}, test r={fold_metrics['pearson_r']:.4f}",
            flush=True,
        )

    fold_metrics_df = pd.DataFrame(fold_metric_rows)
    fold_metrics_df.to_csv(ELASTIC_NET_FOLD_METRICS, index=False)
    coef_long = pd.concat(coef_frames, ignore_index=True)
    coef_long.to_csv(ELASTIC_NET_FOLD_COEF, index=False)
    summary = summarize_fold_coefficients(coef_long)
    summary.to_csv(ELASTIC_NET_COEF_SUMMARY, index=False)
    print("Elastic Net out-of-fold performance:", json.dumps(evaluate(y_log, oof_pred), indent=2), flush=True)
    return summary


def build_stable_scores(stable_coef):
    print("building stable Elastic Net PPS scores from mean fold coefficients...", flush=True)
    stable_coef.to_csv(OUT_STABLE_COEF, index=False)
    stable <- stable_coef[stable_coef["stable_all_10_folds"]].copy()
    data, _ = read_inputs(stable["protein"].tolist())
    score_out = pd.DataFrame(
        {
            "eid": data["eid"].to_numpy(),
            "pain_score": data["pain_score"].to_numpy(),
            "pain_score_log_plus_0p1": np.log(data["pain_score"].to_numpy(dtype=np.float64) + TARGET_OFFSET),
            "weighted_pain_status": data["weighted_pain_status"].to_numpy(),
        }
    )
    score = np.full(len(data), float(stable["mean_intercept_raw"].mean()), dtype=np.float64)
    for row in stable.itertuples(index=False):
        values = data[row.protein].to_numpy(dtype=np.float64, copy=True)
        values[np.isnan(values)] = float(row.mean_impute_value)
        score += values * float(row.mean_coef_raw_imputed)
    score_out["stable_elastic_net_pfs_log"] = score
    score_out["stable_elastic_net_pfs_pain_scale"] = np.exp(score) - TARGET_OFFSET
    score_out["stable_elastic_net_pfs_z"] = (score - score.mean()) / score.std(ddof=0)
    score_out.to_csv(OUT_STABLE_SCORE, index=False)

    summary = {
        "method": "Stable nested-CV Elastic Net PPS using proteins selected in all 10 outer folds and mean raw-space fold coefficients.",
        "regularization": "Elastic Net with both L1 and L2 penalties; l1_ratio candidates are 0.5, 0.75, and 0.9.",
        "target": f"log(pain_score + {TARGET_OFFSET})",
        "random_state": RANDOM_STATE,
        "n_outer_folds": N_OUTER_FOLDS,
        "n_inner_folds": N_INNER_FOLDS,
        "n": int(len(data)),
        "n_stable_proteins": int(len(stable)),
        "performance_in_score_population": evaluate(score_out["pain_score_log_plus_0p1"].to_numpy(), score),
        "score_output": str(OUT_STABLE_SCORE),
        "coefficient_output": str(OUT_STABLE_COEF),
        "fold_coefficients_output": str(ELASTIC_NET_FOLD_COEF),
        "fold_metrics_output": str(ELASTIC_NET_FOLD_METRICS),
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


def main():
    stable_coef = run_nested_elastic_net()
    build_stable_scores(stable_coef)


if __name__ == "__main__":
    main()
