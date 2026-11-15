"""Replacement Figure 4 A/D/E panels; original B/C and author-run F/G are separate."""

import argparse
import hashlib
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.stats import pearsonr


CODE = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE / "results")) / "sensitivity_nonproteomic_weights"
DEFAULT_OUTPUT = CODE.parent / "figures/figure4_independent_weights"


def md5(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "md5").hexdigest()


def configure():
    font_manager.findfont("Arial", fallback_to_default=False)
    matplotlib.rcParams.update({
        "font.family": "Arial", "font.size": 7, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.65, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white", "figure.facecolor": "white",
    })


def save(fig, stem):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for label in fig.findobj(matplotlib.text.Text):
        if label.get_visible() and label.get_text():
            box = label.get_window_extent(renderer)
            if box.width and box.height and (box.x0 < -0.5 or box.y0 < -0.5
                or box.x1 > fig.bbox.x1 + 0.5 or box.y1 > fig.bbox.y1 + 0.5):
                raise ValueError(f"Text outside canvas: {label.get_text()}")
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=600)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(stem.with_name(stem.name + "_preview").with_suffix(".png"), dpi=220)
    plt.close(fig)


def volcano(data, output):
    fig = plt.figure(figsize=(3.55, 2.9), dpi=180)
    ax = fig.add_axes([0.17, 0.18, 0.79, 0.75])
    significant = data.fdr.lt(0.05)
    for selected, color in [(~significant, "#D8D8D8"), (significant, "#28AAC0")]:
        ax.scatter(data.loc[selected, "beta"], -np.log10(data.loc[selected, "p_value"]),
                   s=7, color=color, linewidths=0, rasterized=False)
    ax.axvline(0, color="#777777", linewidth=0.65, linestyle=(0, (4, 4)), zorder=0)
    ax.axhline(-np.log10(0.05), color="#999999", linewidth=0.6, linestyle=(0, (1, 4)))
    ax.set_xlim(-0.026, 0.035)
    ax.set_ylim(0, 94)
    ax.set_xticks([-0.02, -0.01, 0, 0.01, 0.02, 0.03])
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_xlabel(r"$\beta$ coefficient")
    ax.set_ylabel(r"$-\log_{10}(P)$")
    ax.spines[["top", "right"]].set_visible(False)
    for protein, offset in {
        "GAST": (-7, 4), "CHGA": (4, 4), "VSIG2": (-5, 5),
        "ASGR1": (-5, 6), "AMBP": (4, 5), "IL18R1": (-12, 4),
        "GDF15": (-24, 5), "ACE2": (-18, -9),
    }.items():
        row = data.loc[data.protein.eq(protein)].iloc[0]
        ax.annotate(protein, (row.beta, -np.log10(row.p_value)), xytext=offset,
                    textcoords="offset points", fontsize=6.7, fontstyle="italic")
    ax.text(0.025, 0.97, "FDR < 0.05: 834", transform=ax.transAxes, va="top", fontsize=7)
    fig.text(0.018, 0.985, "A", fontsize=11, fontweight="bold", va="top")
    save(fig, output / "Figure4A_independent_pain_score")


def mirrored(shared, output):
    data = shared.sort_values("HR").reset_index(drop=True)
    fig = plt.figure(figsize=(3.7, 4.45), dpi=180)
    left = fig.add_axes([0.19, 0.13, 0.36, 0.81])
    right = fig.add_axes([0.62, 0.13, 0.34, 0.81], sharey=left)
    y = np.arange(len(data))
    left.barh(y, data.HR, color="#F39B85", edgecolor="#777777", linewidth=0.25, height=0.8)
    left.errorbar(data.HR, y, xerr=np.stack([data.HR - data.lower, data.upper - data.HR]),
                   fmt="none", ecolor="#777777", elinewidth=0.45, capsize=0)
    left.axvline(1, color="#333333", linewidth=0.65)
    left.set_yticks(y, labels=data.protein, fontstyle="italic", fontsize=6.4)
    left.set_xlim(0, 1.65)
    left.set_xticks([0, 0.5, 1, 1.5])
    left.set_ylim(len(data) - 0.2, -1.4)
    left.set_xlabel("Hazard ratio (95% CI)", fontsize=7.5)
    left.set_title("PD association", fontsize=7.7, pad=5)
    right.barh(y, data.beta_pain, color="#5DBCC7", edgecolor="#777777", linewidth=0.25, height=0.8)
    right.errorbar(data.beta_pain, y, xerr=data.se, fmt="none", ecolor="#777777", elinewidth=0.45)
    right.axvline(0, color="#333333", linewidth=0.65)
    right.set_xlim(0.014, -0.028)
    right.set_xticks([0.01, 0, -0.01, -0.02])
    right.tick_params(axis="y", left=False, labelleft=False)
    right.set_xlabel(r"$\beta$ coefficient (SE)", fontsize=7.5)
    right.set_title("Pain association", fontsize=7.7, pad=5)
    for ax in [left, right]:
        ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.018, 0.985, "D", fontsize=11, fontweight="bold", va="top")
    save(fig, output / "Figure4D_shared24")


def shared_scatter(shared, output):
    fig = plt.figure(figsize=(3.55, 3.15), dpi=180)
    ax = fig.add_axes([0.19, 0.18, 0.76, 0.76])
    x, y = shared.beta_pain.to_numpy(), shared.HR.to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    grid = np.linspace(x.min(), x.max(), 100)
    ax.plot(grid, intercept + slope * grid, color="#2563A8", linewidth=0.9)
    ax.scatter(x, y, s=18, color="#D55E00", edgecolor="white", linewidth=0.4, zorder=3)
    ax.axhline(1, color="#AAAAAA", linewidth=0.6, linestyle=(0, (3, 3)))
    ax.axvline(0, color="#AAAAAA", linewidth=0.6, linestyle=(0, (3, 3)))
    for protein, offset in {"ITGAV": (3, -9), "EGFR": (3, 5), "NEFL": (-18, 5),
                            "CHGA": (-18, -10), "GPRC5C": (3, 4), "CA14": (-16, 5)}.items():
        row = shared.loc[shared.protein.eq(protein)].iloc[0]
        ax.annotate(protein, (row.beta_pain, row.HR), xytext=offset,
                    textcoords="offset points", fontsize=6.2, fontstyle="italic")
    r, p = pearsonr(x, y)
    ax.text(0.035, 0.96, f"Pearson $r$ = {r:.3f}\n$P$ = 2.26 $\\times$ 10$^{{-7}}$\n24 proteins",
              transform=ax.transAxes, fontsize=7, va="top", linespacing=1.5)
    ax.set_xlabel(r"Pain-score association $\beta$")
    ax.set_ylabel("PD hazard ratio")
    ax.set_xlim(x.min() - 0.003, x.max() + 0.004)
    ax.set_ylim(0.05, 1.55)
    ax.set_yticks([0.25, 0.5, 0.75, 1, 1.25, 1.5])
    ax.set_xticks([-0.02, -0.01, 0, 0.01])
    ax.spines[["top", "right"]].set_visible(False)
    save(fig, output / "FigureS7_shared24_beta_HR")


def mediation(result, joint, paths, joint_paths, output):
    data = result.sort_values("proportion_mediated").reset_index(drop=True)
    fig = plt.figure(figsize=(183 / 25.4, 105 / 25.4), dpi=180)
    diagram = fig.add_axes([0.015, 0.09, 0.65, 0.84])
    forest = fig.add_axes([0.79, 0.13, 0.19, 0.80])
    diagram.set_xlim(0, 1); diagram.set_ylim(0, 1); diagram.axis("off")
    path = paths.set_index("protein").loc[data.protein]
    positions = np.r_[np.linspace(0.94, 0.56, 12), np.linspace(0.44, 0.06, 12)]
    a = path.path_a_standardized_beta.to_numpy()
    b = path.path_b_standardized_beta_logHR_single_protein.to_numpy()
    for i, (row, yy) in enumerate(zip(data.itertuples(), positions)):
        diagram.add_patch(FancyArrowPatch((0.11, 0.5), (0.41, yy), arrowstyle="->",
                         mutation_scale=4, linewidth=0.25 + 0.85 * abs(a[i]) / max(abs(a)),
                         color="#27A4D6", connectionstyle=f"arc3,rad={-0.12 if yy > .5 else .12}"))
        diagram.add_patch(FancyArrowPatch((0.59, yy), (0.89, 0.5), arrowstyle="->",
                         mutation_scale=4, linewidth=0.25 + 0.85 * abs(b[i]) / max(abs(b)),
                         color="#27A4D6", connectionstyle=f"arc3,rad={-0.12 if yy > .5 else .12}"))
        diagram.text(0.5, yy, row.protein, fontsize=6, fontstyle="italic", ha="center", va="center")
    diagram.annotate("", (0.9, 0.5), (0.1, 0.5),
                     arrowprops={"arrowstyle": "->", "color": "#777777", "linewidth": 0.8})
    diagram.text(0.06, 0.5, "Composite\npain score", fontsize=7, ha="center", va="center",
                  bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    diagram.text(0.94, 0.5, "Incident\nPD", fontsize=7, ha="center", va="center",
                  bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    diagram.text(0.16, 0.985, f"Path a\n$\\beta$ = {a.min():.2f}-{a.max():.2f}", fontsize=6.5, va="top", ha="center")
    diagram.text(0.84, 0.985, f"Path b\n$\\beta$ = {b.min():.2f}-{b.max():.2f}", fontsize=6.5, va="top", ha="center")
    total = joint_paths.total_standardized_logHR.iloc[0]
    diagram.text(0.5, 0.51, f"Total $\\beta$ = {total:.2f}", color="#555555", fontsize=7,
                  ha="center", va="bottom", bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    direct = joint_paths.direct_standardized_logHR.iloc[0]
    diagram.text(0.5, 0.485, f"IORW direct $\\beta$ = {direct:.2f}", color="#555555", fontsize=6.3,
                  ha="center", va="top", bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    y = np.arange(len(data))
    forest.hlines(y, data.proportion_lower * 100, data.proportion_upper * 100,
                   color="#333333", linewidth=0.5)
    forest.scatter(data.proportion_mediated * 100, y, color="#2B9DC6", s=12,
                    edgecolors="#333333", linewidths=0.3, zorder=3)
    joint_y = len(data) + 0.2
    row = joint.iloc[0]
    forest.hlines(joint_y, row.proportion_lower * 100, row.proportion_upper * 100,
                   color="#E65A58", linewidth=0.8)
    forest.scatter(row.proportion_mediated * 100, joint_y, color="#D54843", s=17,
                    edgecolors="#333333", linewidths=0.3, zorder=3)
    forest.set_yticks([*y, joint_y], labels=[*data.protein, "Joint 24 proteins"], fontsize=5.7)
    for label in forest.get_yticklabels():
        label.set_fontstyle("italic")
    forest.get_yticklabels()[-1].set_color("#D52D2D")
    forest.get_yticklabels()[-1].set_fontweight("bold")
    lower = min(data.proportion_lower.min(), row.proportion_lower) * 100
    upper = max(data.proportion_upper.max(), row.proportion_upper) * 100
    lo = min(-5, np.floor(lower / 5) * 5)
    hi = np.ceil((upper + 1) / 5) * 5
    forest.set_xlim(lo, hi)
    forest.set_xticks(sorted(set([lo, hi, *np.arange(0, hi + 0.1, 10)])))
    forest.set_ylim(joint_y + 0.8, -0.6)
    forest.axvline(0, color="#777777", linestyle=(0, (4, 3)), linewidth=0.6)
    forest.set_xlabel("Proportion mediated\n(%, 95% bootstrap CI)", fontsize=7)
    forest.spines[["top", "right"]].set_visible(False)
    forest.tick_params(axis="y", length=0, pad=2)
    forest.tick_params(axis="x", labelsize=6)
    fig.text(0.012, 0.985, "E", fontsize=11, fontweight="bold", va="top")
    save(fig, output / "Figure4E_shared24_mediation")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--associations-only", action="store_true")
    args = parser.parse_args()
    configure()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = [args.input_dir / "pain_score_protein_association_results.csv",
             args.input_dir / "shared_proteins_fdr05_concordant.csv"]
    hashes = {path: md5(path) for path in paths}
    pain, shared = [pd.read_csv(path) for path in paths]
    if len(pain) != 2920 or int(pain.fdr.lt(.05).sum()) != 834 or len(shared) != 24 or "COMP" in set(shared.protein):
        raise ValueError("Unexpected new-score protein results")
    volcano(pain, args.output_dir)
    mirrored(shared, args.output_dir)
    shared_scatter(shared, args.output_dir)
    correlation = pearsonr(shared.beta_pain, shared.HR)
    stats = {"n_tested": len(pain), "pain_FDR_significant": 834, "shared_n": 24,
             "shared_lower_PD_risk": int(shared.HR.lt(1).sum()),
             "shared_higher_PD_risk": int(shared.HR.gt(1).sum()),
             "shared_beta_HR_pearson_r": float(correlation.statistic),
             "shared_beta_HR_pearson_p": float(correlation.pvalue)}
    if not args.associations_only:
        mediation_dir = args.input_dir / "mediation24"
        result = pd.read_csv(mediation_dir / "per_protein_mediation_results.csv")
        joint = pd.read_csv(mediation_dir / "joint24_mediation_results.csv")
        path_data = pd.read_csv(mediation_dir / "per_protein_standardized_path_coefficients.csv")
        joint_paths = pd.read_csv(mediation_dir / "joint24_standardized_path_coefficients.csv")
        if not result.n_boot_success.eq(1000).all() or joint.n_boot_success.iloc[0] != 1000:
            raise ValueError("Mediation figure requires all 1,000 bootstrap fits per model")
        mediation(result, joint, path_data, joint_paths, args.output_dir)
        result["indirect_95CI_excludes_zero"] = result.indirect_lower.gt(0) | result.indirect_upper.lt(0)
        table = result[["protein", "n", "events", "proportion_mediated", "proportion_lower",
                        "proportion_upper", "indirect_logHR", "indirect_lower", "indirect_upper",
                        "indirect_95CI_excludes_zero",
                        "n_boot_success"]].copy()
        table = pd.concat([table, pd.DataFrame([{
            "protein": "Joint 24 proteins", "n": joint.n.iloc[0], "events": joint.events.iloc[0],
            "proportion_mediated": joint.proportion_mediated.iloc[0],
            "proportion_lower": joint.proportion_lower.iloc[0], "proportion_upper": joint.proportion_upper.iloc[0],
            "indirect_logHR": joint.indirect_logHR.iloc[0],
            "indirect_lower": joint.indirect_lower.iloc[0], "indirect_upper": joint.indirect_upper.iloc[0],
            "indirect_95CI_excludes_zero": bool(joint.indirect_lower.iloc[0] > 0 or joint.indirect_upper.iloc[0] < 0),
            "n_boot_success": joint.n_boot_success.iloc[0],
        }])], ignore_index=True)
        table.to_csv(args.output_dir / "TableS11_shared24_mediation.csv", index=False)
        warnings = pd.read_csv(mediation_dir / "bootstrap_counts.csv").warnings.sum()
        legend = (
            "Replacement Figure 4 panels based on independently weighted pain scores.\n\n"
            "A, Protein-wise linear regression analyses identified 834 proteins associated with pain score at BH FDR<0.05 "
            "among 2,920 tested proteins (48,586 base complete cases; protein-specific sample sizes vary). The pain score "
            "is the dependent variable and standardized protein is the predictor. Teal points denote FDR-significant associations.\n\n"
            "D, Twenty-four proteins show directionally concordant associations with independently weighted pain score and incident PD. "
            "Orange bars show unchanged protein-PD hazard ratios with 95% CIs; teal bars show the new pain-protein beta estimates "
            "with SEs. Twenty-one proteins are associated with lower PD risk and three with higher risk.\n\n"
            "E, IORW mediation analyses of the 24 proteins use 46,591 participants and 707 incident PD events. Proteins are "
            "direction-aligned for mediation. Diagram a coefficients are from covariate-adjusted standardized mediator-on-pain "
            "linear models; b coefficients are from single-protein Cox models adjusted for standardized pain and covariates. "
            "Line widths reflect absolute association coefficients, not evidence of causality. Total and IORW direct beta estimates "
            "are per one SD of pain score. Forest points show observed proportions mediated and horizontal lines show 95% "
            "percentile CIs from 1,000 participant-level bootstrap resamples per model. The red point denotes the joint24 model. "
            f"For {int(result.indirect_95CI_excludes_zero.sum())} of 24 proteins, the indirect log-HR CI excludes zero; "
            "this is an unadjusted CI criterion, not one-sided or FDR-corrected significance. "
            f"Among 25,000 attempted bootstrap fits, {int(warnings):,} produced warnings about possible infinite nuisance coefficients; "
            "finite estimates were retained following the original implementation and representative warned fits were checked "
            "against the original formula models. These model-stability limitations require review.\n\n"
            "B and C are unchanged. F uses the author-supplied 24-protein PPI network. G is plotted separately by "
            "07_plot_figure4G_top5.py from the supplied 24-protein enrichment exports. "
            "No mixed old/new full Figure 4 has been assembled.\n"
        )
        (args.output_dir / "Figure4_replacement_panel_legends.txt").write_text(legend, encoding="utf-8")
    if any(md5(path) != digest for path, digest in hashes.items()):
        raise ValueError("Protein result inputs changed")
    pain.to_csv(args.output_dir / "Figure4A_source_data.csv", index=False)
    shared.to_csv(args.output_dir / "Figure4D_source_data.csv", index=False)
    (args.output_dir / "shared24_for_PPI_enrichment.txt").write_text("\n".join(sorted(shared.protein)) + "\n", encoding="utf-8")
    (args.output_dir / "tested_protein_background_2920.txt").write_text("\n".join(sorted(pain.protein)) + "\n", encoding="utf-8")
    (args.output_dir / "association_panel_statistics.json").write_text(json.dumps(stats, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    print(f"Replacement panels: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
