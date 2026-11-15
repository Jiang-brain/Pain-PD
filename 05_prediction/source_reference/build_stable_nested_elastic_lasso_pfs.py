import json
import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PYDEPS = ROOT / "pydeps"
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


PROTEIN_PATH = ROOT / "protein_new_missing60_filtered.csv"
SCORE_PATH = ROOT / "composite_pain_scores.csv"

LASSO_FOLD_COEF = ROOT / "pain_related_proteomic_score_nested_lasso_fold_coefficients.csv"
ELASTIC_FOLD_COEF = ROOT / "pain_related_proteomic_score_nested_elasticnet_fold_coefficients.csv"
ELASTIC_COEF_SUMMARY = ROOT / "pain_related_proteomic_score_nested_elasticnet_coefficient_summary.csv"
ELASTIC_FOLD_METRICS = ROOT / "pain_related_proteomic_score_nested_elasticnet_fold_metrics.csv"
OUT_STABLE_COEF = ROOT / "pain_related_proteomic_score_stable_nested_coefficients.csv"
OUT_STABLE_SCORE = ROOT / "pain_related_proteomic_score_stable_nested_pfs.csv"
OUT_SUMMARY = ROOT / "pain_related_proteomic_score_stable_nested_summary.json"

RANDOM_STATE = 20260525
N_OUTER_FOLDS = 10
N_INNER_FOLDS = 10
TARGET_OFFSET = 0.1


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


def make_elastic_pipeline(inner_cv, fold):
    model = ElasticNetCV(
        l1_ratio=[0.5, 0.75, 0.9],
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


def run_nested_elasticnet():
    if ELASTIC_FOLD_COEF.exists() and ELASTIC_COEF_SUMMARY.exists():
        print("reusing existing nested ElasticNet outputs", flush=True)
        return
    warnings.filterwarnings("ignore", category=ConvergenceWarning)
    print("reading full protein matrix for nested ElasticNet...", flush=True)
    data, protein_cols = read_inputs()
    x = data[protein_cols]
    y_log = np.log(data["pain_score"].to_numpy(dtype=np.float64) + TARGET_OFFSET)
    outer_cv = KFold(n_splits=N_OUTER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    coef_frames = []
    fold_metric_rows = []
    oof_pred = np.full(len(data), np.nan, dtype=np.float64)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(x), start=1):
        print(f"elastic outer fold {fold}/{N_OUTER_FOLDS}: fitting inner {N_INNER_FOLDS}-fold ElasticNetCV...", flush=True)
        inner_cv = KFold(n_splits=N_INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE + fold)
        pipe = make_elastic_pipeline(inner_cv, fold)
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
            f"elastic fold {fold}: alpha={model.alpha_:.6g}, l1_ratio={model.l1_ratio_}, "
            f"selected={int(selected.sum())}, test r={fold_metrics['pearson_r']:.4f}",
            flush=True,
        )

    fold_metrics_df = pd.DataFrame(fold_metric_rows)
    fold_metrics_df.to_csv(ELASTIC_FOLD_METRICS, index=False)
    coef_long = pd.concat(coef_frames, ignore_index=True)
    coef_long.to_csv(ELASTIC_FOLD_COEF, index=False)
    summary = summarize_fold_coefficients(coef_long, "elasticnet")
    summary.to_csv(ELASTIC_COEF_SUMMARY, index=False)
    print("elastic OOF:", json.dumps(evaluate(y_log, oof_pred), indent=2), flush=True)


def summarize_fold_coefficients(coef_long, method):
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
        )
    )
    if "l1_ratio" in coef_long.columns:
        l1 = coef_long.groupby("protein", as_index=False).agg(mean_l1_ratio=("l1_ratio", "mean"))
        out = out.merge(l1, on="protein", how="left")
    out["method"] = method
    out["selection_frequency"] = out["selection_count"] / N_OUTER_FOLDS
    out["stable_all_10_folds"] = out["selection_count"] == N_OUTER_FOLDS
    out = out.sort_values(["stable_all_10_folds", "selection_count", "mean_abs_coef_standardized"], ascending=[False, False, False])
    return out


def build_stable_scores():
    print("building stable PFS scores from mean fold coefficients...", flush=True)
    lasso_long = pd.read_csv(LASSO_FOLD_COEF)
    elastic_long = pd.read_csv(ELASTIC_FOLD_COEF)
    summaries = [
        summarize_fold_coefficients(lasso_long, "lasso"),
        summarize_fold_coefficients(elastic_long, "elasticnet"),
    ]
    stable_coef = pd.concat(summaries, ignore_index=True)
    stable_coef.to_csv(OUT_STABLE_COEF, index=False)

    stable_by_method = {
        method: df[df["stable_all_10_folds"]].copy()
        for method, df in stable_coef.groupby("method", sort=False)
    }
    union_proteins = sorted({p for df in stable_by_method.values() for p in df["protein"].tolist()})
    data, _ = read_inputs(union_proteins)
    score_out = pd.DataFrame(
        {
            "eid": data["eid"].to_numpy(),
            "pain_score": data["pain_score"].to_numpy(),
            "pain_score_log_plus_0p1": np.log(data["pain_score"].to_numpy(dtype=np.float64) + TARGET_OFFSET),
            "weighted_pain_status": data["weighted_pain_status"].to_numpy(),
        }
    )
    summary_rows = []
    for method, coef_df in stable_by_method.items():
        score = np.full(len(data), float(coef_df["mean_intercept_raw"].mean()), dtype=np.float64)
        for row in coef_df.itertuples(index=False):
            values = data[row.protein].to_numpy(dtype=np.float64, copy=True)
            values[np.isnan(values)] = float(row.mean_impute_value)
            score += values * float(row.mean_coef_raw_imputed)
        score_out[f"stable_{method}_pfs_log"] = score
        score_out[f"stable_{method}_pfs_pain_scale"] = np.exp(score) - TARGET_OFFSET
        score_out[f"stable_{method}_pfs_z"] = (score - score.mean()) / score.std(ddof=0)
        perf = evaluate(score_out["pain_score_log_plus_0p1"].to_numpy(), score)
        summary_rows.append(
            {
                "method": method,
                "n_stable_proteins": int(len(coef_df)),
                "mean_intercept_raw": float(coef_df["mean_intercept_raw"].mean()),
                **perf,
            }
        )
    score_out.to_csv(OUT_STABLE_SCORE, index=False)

    summary = {
        "method": "Stable nested-CV PFS using proteins selected in all 10 outer folds and mean raw-space fold coefficients.",
        "target": f"log(pain_score + {TARGET_OFFSET})",
        "n": int(len(data)),
        "performance": summary_rows,
        "score_output": str(OUT_STABLE_SCORE),
        "coefficient_output": str(OUT_STABLE_COEF),
        "elastic_fold_coefficients_output": str(ELASTIC_FOLD_COEF),
    }
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


def main():
    if not LASSO_FOLD_COEF.exists():
        raise FileNotFoundError(f"Missing existing nested LASSO coefficients: {LASSO_FOLD_COEF}")
    run_nested_elasticnet()
    build_stable_scores()


if __name__ == "__main__":
    main()
