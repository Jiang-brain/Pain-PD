"""Plot completed bootstrap score correlations; no bootstrap/model refitting."""

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


CODE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = CODE_ROOT / "results/sensitivity_nonproteomic_weights/bootstrap1000"
DEFAULT_OUTPUT = CODE_ROOT.parent / "figures/bootstrap_spearman_stability"
CONTRACT = {
    "conclusion": (
        "Bootstrap-derived composite pain scores show substantial overall rank "
        "agreement, with variability across pairs of weight estimates."
    ),
    "archetype": "quantitative grid",
    "backend": "Python / matplotlib only",
    "population": "All 49,056 fixed proteomic scoring participants",
    "derivation": "1,000 bootstrap resamples of 408,583 non-proteomic participants",
    "panel_A": "Full Spearman matrix in the original replicate order; diagonal masked",
    "panel_B": "Distribution of the 499,500 unique off-diagonal correlations",
    "size": "183 x 80 mm",
    "export": "600-dpi PNG/TIFF; PDF/SVG with editable text",
    "statistics": "Median and empirical 2.5th/97.5th percentiles; no inferential tests",
    "review_risks": [
        "Pairwise correlations share replicates and are dependent.",
        "The percentile range is not a confidence interval.",
        "The population includes participants with zero scores under every weight set.",
        "Overall rank agreement does not establish stable individual coefficients.",
        "No Pearson correlations or pain-only subgroup results appear in this figure.",
    ],
    "image_integrity": "No clustering, score adjustment, smoothing or model refitting",
}


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load_and_validate(input_dir):
    paths = {
        "matrix": input_dir / "correlation_matrix_all_proteomic_spearman.csv.gz",
        "summary": input_dir / "correlation_summary.csv",
        "diagnostics": input_dir / "bootstrap_fit_diagnostics.csv",
        "checks": input_dir / "bootstrap_validation_checks.csv",
    }
    hashes = {path: sha256(path) for path in paths.values()}
    table = pd.read_csv(paths["matrix"], index_col=0, float_precision="round_trip")
    expected = [f"bootstrap_{index:04d}" for index in range(1, 1001)]
    if list(table.index) != expected or list(table.columns) != expected:
        raise ValueError("Expected exactly 1,000 labelled replicates in saved order")
    matrix = table.to_numpy(dtype=float)
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, atol=1e-12, rtol=0):
        raise ValueError("Non-finite/asymmetric correlation matrix")
    if not np.allclose(np.diag(matrix), 1, atol=1e-12, rtol=0):
        raise ValueError("Invalid self-correlations")
    if matrix.min() < -1 or matrix.max() > 1:
        raise ValueError("Correlation out of range")
    values = matrix[np.tril_indices(1000, k=-1)]
    if len(values) != 499500 or values.min() < 0.5:
        raise ValueError("Unexpected pair count or values outside the declared display range")
    diagnostics = pd.read_csv(paths["diagnostics"])
    if len(diagnostics) != 1000 or not diagnostics.valid.astype(str).str.lower().eq("true").all():
        raise ValueError("Not all 1,000 bootstrap fits were valid")
    if not diagnostics.n_sampled.eq(408583).all():
        raise ValueError("Unexpected bootstrap derivation sample size")
    checks = pd.read_csv(paths["checks"])
    if not checks.passed.astype(str).str.lower().eq("true").all():
        raise ValueError("Existing bootstrap validation has failed checks")
    summary = pd.read_csv(paths["summary"], float_precision="round_trip")
    saved = summary.loc[
        summary.comparison.eq("bootstrap_vs_bootstrap")
        & summary.population.eq("all_proteomic") & summary.method.eq("spearman")
    ]
    if len(saved) != 1 or saved.iloc[0].n_correlations != len(values):
        raise ValueError("Missing/ambiguous saved Spearman summary")
    stats = {
        "n_bootstraps": 1000, "n_scoring_participants": 49056,
        "n_derivation_participants": 408583, "n_correlations": len(values),
        "median": float(np.median(values)), "p025": float(np.quantile(values, 0.025)),
        "p975": float(np.quantile(values, 0.975)), "minimum": float(values.min()),
        "maximum": float(values.max()), "mean": float(values.mean()),
    }
    for field in ["median", "p025", "p975", "minimum", "maximum", "mean"]:
        if abs(stats[field] - float(saved.iloc[0][field])) > 1e-12:
            raise ValueError(f"Calculated {field} does not reproduce the saved summary")
    return matrix, values, stats, hashes


def set_style():
    font_manager.findfont("Arial", fallback_to_default=False)
    matplotlib.rcParams.update({
        "font.family": "Arial", "font.size": 7,
        "axes.labelsize": 7.5, "axes.titlesize": 8,
        "axes.linewidth": 0.7, "axes.grid": False,
        "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.major.width": 0.65, "ytick.major.width": 0.65,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "svg.hashsalt": "bootstrap-spearman-stability",
    })


def make_figure(matrix, values, stats):
    fig = plt.figure(figsize=(183 / 25.4, 80 / 25.4), dpi=180)
    heatmap_ax = fig.add_axes([0.07, 0.18, 0.302, 0.691])
    colorbar_ax = fig.add_axes([0.394, 0.18, 0.018, 0.691])
    histogram_ax = fig.add_axes([0.58, 0.18, 0.39, 0.691])
    fig.text(0.07, 0.982, "1,000 bootstrap scores in 49,056 proteomic participants",
             fontsize=9, va="top")
    fig.text(0.012, 0.919, "A", fontsize=10, fontweight="bold", va="top")
    fig.text(0.522, 0.919, "B", fontsize=10, fontweight="bold", va="top")
    cmap = LinearSegmentedColormap.from_list(
        "pain_rank_agreement", ["#F3F5F5", "#BBDDD5", "#56AEA0", "#145E53"]
    )
    cmap.set_bad("white")
    displayed = np.ma.array(matrix, mask=np.eye(1000, dtype=bool))
    image = heatmap_ax.imshow(displayed, cmap=cmap, vmin=0.5, vmax=1.0,
                              origin="upper", extent=(0.5, 1000.5, 1000.5, 0.5),
                              interpolation="none", rasterized=True)
    heatmap_ax.set_title("Pairwise score correlations", loc="left", pad=7)
    heatmap_ax.set_xlabel("Bootstrap replicate", labelpad=3)
    heatmap_ax.set_ylabel("Bootstrap replicate", labelpad=3)
    heatmap_ax.set_xticks([1, 250, 500, 750, 1000])
    heatmap_ax.set_yticks([1, 250, 500, 750, 1000])
    colorbar = fig.colorbar(image, cax=colorbar_ax, ticks=np.linspace(0.5, 1, 6))
    colorbar.set_label(r"Spearman $\rho$", labelpad=4, fontsize=7.5)
    colorbar.outline.set_linewidth(0.6)
    colorbar.ax.tick_params(labelsize=6.5, length=2)

    bins = np.linspace(0.5, 1, 51)
    counts, edges = np.histogram(values, bins=bins)
    histogram_ax.axvspan(stats["p025"], stats["p975"], color="#E5EBE9", alpha=0.7, zorder=0)
    histogram_ax.hist(values, bins=bins, weights=np.full(len(values), 100 / len(values)),
                      color="#2A9D8F", edgecolor="white", linewidth=0.25, zorder=2)
    histogram_ax.axvline(stats["median"], color="#D55E00", linewidth=1.2, zorder=4)
    for value in [stats["p025"], stats["p975"]]:
        histogram_ax.axvline(value, color="#666666", linewidth=0.8,
                             linestyle=(0, (3, 2)), zorder=3)
    histogram_ax.set_title("Distribution of pairwise correlations", loc="left", pad=7)
    histogram_ax.set_xlabel(r"Spearman $\rho$", labelpad=3)
    histogram_ax.set_ylabel("Pairwise comparisons (%)", labelpad=3)
    histogram_ax.set_xlim(0.5, 1)
    upper = float(np.ceil(counts.max() / len(values) * 100 * 1.14 / 2) * 2)
    histogram_ax.set_ylim(0, upper)
    histogram_ax.set_xticks(np.linspace(0.5, 1, 6))
    histogram_ax.set_yticks(np.arange(0, upper + 0.01, 2))
    histogram_ax.spines[["top", "right"]].set_visible(False)
    histogram_ax.text(0.035, 0.96, f"Median $\\rho$ = {stats['median']:.3f}",
                      transform=histogram_ax.transAxes, fontsize=7.3, va="top", color="#A64A00")
    histogram_ax.text(0.035, 0.855,
                      f"2.5th-97.5th percentiles\n{stats['p025']:.3f}-{stats['p975']:.3f}\n\n499,500 unique pairs",
                      transform=histogram_ax.transAxes, fontsize=6.7, va="top", linespacing=1.45)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if box.width and box.height and (
                box.x0 < -0.5 or box.y0 < -0.5
                or box.x1 > fig.bbox.x1 + 0.5 or box.y1 > fig.bbox.y1 + 0.5
            ):
                raise ValueError(f"Text outside export canvas: {text.get_text()}")
    histogram_data = pd.DataFrame({
        "bin_lower": edges[:-1], "bin_upper": edges[1:], "n_pairs": counts,
        "percentage_of_pairs": counts / len(values) * 100,
    })
    return fig, histogram_data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    matrix, values, stats, hashes = load_and_validate(args.input_dir)
    set_style()
    fig, histogram = make_figure(matrix, values, stats)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / "Bootstrap1000_Spearman_stability"
    for extension in ["png", "pdf", "svg"]:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=600)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600,
                pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(args.output_dir / "Bootstrap1000_Spearman_preview.png", dpi=230)
    plt.close(fig)
    pd.DataFrame([stats]).to_csv(args.output_dir / "Spearman_summary.csv", index=False)
    histogram.to_csv(args.output_dir / "Spearman_histogram_source_data.csv", index=False)
    np.savez_compressed(args.output_dir / "Spearman_matrix_source_data.npz",
                        spearman_correlation=matrix, bootstrap_replicate=np.arange(1, 1001))
    (args.output_dir / "figure_contract.json").write_text(json.dumps(CONTRACT, indent=2), encoding="utf-8")
    report = (
        "RESULTS\n"
        "Across 1,000 bootstrap-derived composite pain scores evaluated in the "
        "same 49,056 proteomic participants, the median pairwise Spearman rank "
        f"correlation was {stats['median']:.3f} "
        f"(2.5th-97.5th percentiles, {stats['p025']:.3f}-{stats['p975']:.3f}; "
        "499,500 unique pairs), indicating substantial overall rank agreement "
        "with variability across bootstrap estimates.\n\n"
        "FIGURE LEGEND\n"
        "Bootstrap stability of composite pain scores. Each of 1,000 bootstrap "
        "samples of the 408,583 non-proteomic derivation participants was used "
        "to re-estimate the pain weights. Negative pain coefficients were "
        "truncated to zero and each replicate's weights were normalized before "
        "being applied to the same 49,056 proteomic participants. A, Spearman "
        "rank correlation matrix among the 1,000 score vectors, in their saved "
        "replicate order; diagonal self-correlations are masked. The colour "
        "scale spans 0.5-1.0. B, Distribution of all 499,500 unique pairwise "
        "Spearman correlations, excluding self-correlations. Histogram heights "
        "show percentages of all pairs. The orange line marks the median; "
        "gray dashed lines bound the shaded empirical 2.5th-97.5th percentile "
        "range, which is not a confidence interval. Pairwise correlations are "
        "dependent because they share bootstrap replicates. The scoring cohort "
        "includes participants with zero scores under all weight sets; this "
        "figure is not restricted to participants with pain. Correlation "
        "agreement does not establish individual-coefficient stability or "
        "external predictive validity.\n"
    )
    (args.output_dir / "Spearman_results_and_legend.txt").write_text(report, encoding="utf-8")
    if any(sha256(path) != expected for path, expected in hashes.items()):
        raise ValueError("Saved bootstrap inputs changed during plotting")
    qa = {
        "all_1000_fits_valid": True, "matrix_symmetric_and_finite": True,
        "unique_pairs": len(values), "saved_summary_reproduced": True,
        "inputs_unchanged": True, "only_spearman_all_proteomic": True,
        "percentiles_not_labelled_confidence_intervals": True,
        "all_text_within_canvas": True, "raster_dpi": 600,
        "input_hashes": {path.name: value for path, value in hashes.items()},
    }
    (args.output_dir / "figure_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    print(json.dumps(stats, indent=2))
    print(f"Validated figure exports: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
