import math
import re
import sys
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = Path(__import__("os").environ.get("PD_PAIN_DATA_DIR", CODE_ROOT / "data"))
RESULTS_DIR = Path(__import__("os").environ.get("PD_PAIN_RESULTS_DIR", CODE_ROOT / "results"))
PYDEPS = CODE_ROOT / "pydeps"
if PYDEPS.exists():
    sys.path.insert(0, str(PYDEPS))

import numpy as np
import pandas as pd


PPMI_PATH = DATA_DIR / "results_PPMI.csv"
UKB_ALL_PATH = RESULTS_DIR / "protein_PD_cox_results.csv"
UKB_SIG_PATH = RESULTS_DIR / "protein_PD_cox_significant_FDR05.csv"
OUT_DIR = RESULTS_DIR / "protein_result_figures"
OUT_DIR.mkdir(exist_ok=True)


def parse_first_float(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN"}:
        return np.nan
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return float(match.group(0)) if match else np.nan


def bh_fdr(p_values):
    p = pd.Series(p_values, dtype="float64")
    out = pd.Series(np.nan, index=p.index, dtype="float64")
    valid = p.dropna()
    if valid.empty:
        return out

    order = valid.sort_values().index
    ranked = valid.loc[order].to_numpy()
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out.loc[order] = adjusted
    return out


def read_ppmi(path):
    raw = pd.read_csv(path, header=None)
    data = raw.iloc[4:].copy()
    data = data.rename(
        columns={
            0: "protein",
            19: "ppmi_or_ci",
            20: "ppmi_p",
            21: "ppmi_p_bonferroni",
        }
    )
    data = data[["protein", "ppmi_or_ci", "ppmi_p", "ppmi_p_bonferroni"]]
    data["protein"] = data["protein"].astype(str).str.strip()
    data = data[data["protein"].ne("") & data["protein"].ne("nan")].copy()
    data["ppmi_or"] = data["ppmi_or_ci"].map(parse_first_float)
    data["ppmi_p"] = pd.to_numeric(data["ppmi_p"], errors="coerce")
    data["ppmi_p_bonferroni"] = pd.to_numeric(data["ppmi_p_bonferroni"], errors="coerce")
    data["ppmi_fdr_all_ppmi"] = bh_fdr(data["ppmi_p"])
    data["ppmi_log_or"] = np.log(data["ppmi_or"])
    return data


def corr_stats(x, y):
    from scipy import stats

    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return {
            "n": int(mask.sum()),
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "spearman_r": np.nan,
            "spearman_p": np.nan,
        }
    pearson = stats.pearsonr(x[mask], y[mask])
    spearman = stats.spearmanr(x[mask], y[mask])
    return {
        "n": int(mask.sum()),
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_r": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def direction_summary(df):
    x = df["beta"].to_numpy(dtype="float64")
    y = df["ppmi_log_or"].to_numpy(dtype="float64")
    mask = np.isfinite(x) & np.isfinite(y) & (x != 0) & (y != 0)
    if mask.sum() == 0:
        return {"n_direction_tested": 0, "n_same_direction": 0, "pct_same_direction": np.nan}
    same = np.sign(x[mask]) == np.sign(y[mask])
    return {
        "n_direction_tested": int(mask.sum()),
        "n_same_direction": int(same.sum()),
        "pct_same_direction": float(same.mean() * 100),
    }


def main():
    ppmi = read_ppmi(PPMI_PATH)
    ukb_all = pd.read_csv(UKB_ALL_PATH)
    ukb_sig = pd.read_csv(UKB_SIG_PATH)

    ukb_sig = ukb_sig[ukb_sig["fdr"] < 0.05].copy()
    ukb_sig["ukb_log_hr"] = np.log(ukb_sig["HR"])

    matched = ukb_sig.merge(ppmi, on="protein", how="left", indicator=True)
    matched["has_ppmi"] = matched["_merge"].eq("both") & matched["ppmi_p"].notna() & matched["ppmi_or"].notna()

    matched_with_ppmi = matched[matched["has_ppmi"]].copy()
    matched_with_ppmi["ppmi_fdr_matched_ukb205"] = bh_fdr(matched_with_ppmi["ppmi_p"])
    matched.loc[matched_with_ppmi.index, "ppmi_fdr_matched_ukb205"] = matched_with_ppmi["ppmi_fdr_matched_ukb205"]

    nominal = matched_with_ppmi[matched_with_ppmi["ppmi_p"] < 0.05].copy()
    fdr_all_ppmi = matched_with_ppmi[matched_with_ppmi["ppmi_fdr_all_ppmi"] < 0.05].copy()
    fdr_matched = matched_with_ppmi[matched_with_ppmi["ppmi_fdr_matched_ukb205"] < 0.05].copy()
    bonferroni = matched_with_ppmi[matched_with_ppmi["ppmi_p_bonferroni"] < 0.05].copy()

    corr_15_claim = corr_stats(fdr_all_ppmi["beta"].to_numpy(), fdr_all_ppmi["ppmi_log_or"].to_numpy())
    corr_23_matched_fdr = corr_stats(fdr_matched["beta"].to_numpy(), fdr_matched["ppmi_log_or"].to_numpy())
    corr_nominal = corr_stats(nominal["beta"].to_numpy(), nominal["ppmi_log_or"].to_numpy())
    corr_all_matched = corr_stats(matched_with_ppmi["beta"].to_numpy(), matched_with_ppmi["ppmi_log_or"].to_numpy())
    dir_15 = direction_summary(fdr_all_ppmi)
    dir_nominal = direction_summary(nominal)

    summary = pd.DataFrame(
        [
            {
                "ukb_fdr05_proteins": len(ukb_sig),
                "ukb_fdr05_with_ppmi_result": int(matched["has_ppmi"].sum()),
                "ppmi_nominal_p_lt_0_05_among_matched": len(nominal),
                "ppmi_bh_fdr_lt_0_05_among_all_ppmi_proteins": len(fdr_all_ppmi),
                "ppmi_bh_fdr_lt_0_05_within_ukb205_matched": len(fdr_matched),
                "ppmi_bonferroni_lt_0_05_among_matched": len(bonferroni),
                "corr_set_for_claim": "PPMI BH-FDR<0.05 across all PPMI proteins",
                "corr_n": corr_15_claim["n"],
                "pearson_r_ukb_logHR_vs_ppmi_logOR": corr_15_claim["pearson_r"],
                "pearson_p": corr_15_claim["pearson_p"],
                "spearman_r": corr_15_claim["spearman_r"],
                "spearman_p": corr_15_claim["spearman_p"],
                "fdr15_n_direction_tested": dir_15["n_direction_tested"],
                "fdr15_n_same_direction": dir_15["n_same_direction"],
                "fdr15_pct_same_direction": dir_15["pct_same_direction"],
                "matched_only_fdr23_pearson_r": corr_23_matched_fdr["pearson_r"],
                "matched_only_fdr23_pearson_p": corr_23_matched_fdr["pearson_p"],
                "all_matched_pearson_r": corr_all_matched["pearson_r"],
                "all_matched_pearson_p": corr_all_matched["pearson_p"],
                "nominal_ppmi_pearson_r": corr_nominal["pearson_r"],
                "nominal_ppmi_pearson_p": corr_nominal["pearson_p"],
                "nominal_n_same_direction": dir_nominal["n_same_direction"],
                "nominal_pct_same_direction": dir_nominal["pct_same_direction"],
            }
        ]
    )

    matched_out = OUT_DIR / "ppmi_ukb_pd_protein_comparison_matched.csv"
    summary_out = OUT_DIR / "ppmi_ukb_pd_protein_comparison_summary.csv"
    fdr15_out = OUT_DIR / "ppmi_ukb_pd_protein_comparison_fdr15.csv"
    fdr23_out = OUT_DIR / "ppmi_ukb_pd_protein_comparison_fdr23_matched_only.csv"

    matched.drop(columns=["_merge"], errors="ignore").to_csv(matched_out, index=False)
    summary.to_csv(summary_out, index=False)
    fdr_all_ppmi.to_csv(fdr15_out, index=False)
    fdr_matched.to_csv(fdr23_out, index=False)

    print(summary.to_string(index=False))
    print("\nFDR-significant in both datasets using PPMI BH-FDR across all PPMI proteins:")
    print(fdr_all_ppmi[["protein", "beta", "HR", "ppmi_or", "ppmi_p", "ppmi_fdr_all_ppmi"]].to_string(index=False))
    print("\nOutputs:")
    print(matched_out)
    print(summary_out)
    print(fdr15_out)
    print(fdr23_out)


if __name__ == "__main__":
    main()
