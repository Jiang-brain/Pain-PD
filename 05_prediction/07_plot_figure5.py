"""Redraw Figure 5 from saved OOF risks; never refit or alter prediction models."""

import argparse
import hashlib
import json
import os
import platform
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd
import sklearn
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve


CODE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE_ROOT / "results")) / "sensitivity_nonproteomic_weights/prediction"
DEFAULT_OUTPUT = CODE_ROOT.parent / "figures/figure5_new_prediction"
FEATURES = (
    ("covariates_only", "Covariates only", "#575757"),
    ("new_stable_panel", "370 proteins only", "#2A9D8F"),
    ("covariates_plus_new_stable_panel", "Covariates + 370 proteins", "#D55E00"),
)
PANELS = (
    {
        "letter": "A", "task": "incident_pd_full_followup",
        "title": "PD onset during full follow-up",
        "ylim": (0.65, 0.90), "yticks": (0.7, 0.8, 0.9),
        "window": "After exclusion of prevalent PD and PD diagnosed within 2 years",
    },
    {
        "letter": "B", "task": "incident_pd_2to5y",
        "title": "PD onset within 5 years",
        "ylim": (0.58, 1.0), "yticks": (0.7, 0.8, 0.9, 1.0),
        "window": "PD diagnosed after 2 and within 5 years",
    },
    {
        "letter": "C", "task": "incident_pd_after5y",
        "title": "PD onset after 5 years",
        "ylim": (0.64, 0.88), "yticks": (0.65, 0.75, 0.85),
        "window": "PD diagnosed after 5 years, with no upper 10-year limit",
    },
)
CONTRACT = {
    "conclusion": (
        "The pain-related 370-protein panel carries internally evaluated PD "
        "discrimination information; performance differs by prediction window, "
        "and adding covariates does not improve every panel-only comparison."
    ),
    "archetype": "quantitative grid; three equal ROC panels with inset AUC bars",
    "backend": "Python / matplotlib only",
    "panel_evidence": {
        "A": "Overall incident PD discrimination",
        "B": "Near-term PD discrimination after the two-year exclusion",
        "C": "PD discrimination beyond five years, not a within-ten-year horizon",
    },
    "export": "183 mm wide; white background; 600-dpi PNG/TIFF; editable SVG/PDF",
    "style": "Match the existing gray/teal/orange triptych and inset-bar layout",
    "statistics": "Pooled stratified 10-fold OOF ROC and AUC; no invented intervals",
    "review_risks": [
        "Protein selection preceded PD CV; this is conditional internal evaluation.",
        "The multivariable protein panel is not a single scalar PPS predictor.",
        "Insets use visibly truncated AUC axes; ROC axes remain 0-1 in every panel.",
        "B excludes PD diagnosed within two years; C has no upper ten-year cap.",
    ],
    "image_integrity": "All curves are recomputed from unchanged saved OOF risks; no smoothing.",
}


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load_and_validate(input_dir):
    predictions_path = input_dir / "lightgbm_new/oof_predictions.csv"
    metrics_path = input_dir / "lightgbm_new/cv_metrics.csv"
    summary_path = input_dir / "new_elastic_net_summary.json"
    sources = [predictions_path, metrics_path, summary_path]
    before = {path: digest(path) for path in sources}
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary["n_stable_new"] != 370:
        raise ValueError("This figure is specified for the verified 370-protein panel")
    data = pd.read_csv(predictions_path, float_precision="round_trip")
    metrics = pd.read_csv(metrics_path, float_precision="round_trip")
    required = {"eid", "task", "feature_set", "fold", "outcome", "predicted_risk"}
    if not required.issubset(data.columns):
        raise ValueError(f"Missing prediction fields: {sorted(required - set(data.columns))}")
    results, checks = {}, []
    for panel in PANELS:
        task = panel["task"]
        reference = None
        for feature, label, color in FEATURES:
            sub = data.loc[data.task.eq(task) & data.feature_set.eq(feature)].copy()
            if sub.empty or sub[list(required)].isna().any().any():
                raise ValueError(f"Missing/incomplete OOF predictions: {task}/{feature}")
            if sub.eid.duplicated().any():
                raise ValueError(f"Duplicate participant OOF predictions: {task}/{feature}")
            if set(sub.fold) != set(range(1, 11)) or set(sub.outcome) != {0, 1}:
                raise ValueError(f"Invalid fold/event encoding: {task}/{feature}")
            if not np.isfinite(sub.predicted_risk).all() or not sub.predicted_risk.between(0, 1).all():
                raise ValueError(f"Invalid predicted probabilities: {task}/{feature}")
            cohort = sub.set_index("eid")[["outcome", "fold"]].sort_index()
            if reference is None:
                reference = cohort
            elif not reference.equals(cohort):
                raise ValueError(f"Participant, outcome or CV-fold mismatch: {task}/{feature}")
            y, probability = sub.outcome.to_numpy(), sub.predicted_risk.to_numpy()
            fpr, tpr, _ = roc_curve(y, probability, drop_intermediate=False)
            auc = float(roc_auc_score(y, probability))
            if abs(float(np.trapezoid(tpr, fpr)) - auc) > 1e-12:
                raise ValueError("ROC integral does not equal rank-based AUC")
            saved = metrics.loc[metrics.task.eq(task) & metrics.feature_set.eq(feature)]
            if len(saved) != 1:
                raise ValueError(f"Missing/ambiguous saved metrics: {task}/{feature}")
            saved = saved.iloc[0]
            difference = abs(auc - float(saved.auc))
            if difference > 1e-12 or len(sub) != saved.n_total or y.sum() != saved.events_total:
                raise ValueError(f"OOF metrics do not reproduce saved results: {task}/{feature}")
            key = (task, feature)
            results[key] = {
                "panel": panel["letter"], "task": task, "task_label": panel["title"],
                "feature_set": feature, "model": label, "color": color,
                "n": len(sub), "events": int(y.sum()), "auc": auc,
                "auprc": float(average_precision_score(y, probability)),
                "fpr": fpr, "tpr": tpr,
            }
            checks.append({
                "panel": panel["letter"], "feature_set": feature,
                "auc": auc, "saved_auc": float(saved.auc),
                "absolute_difference": difference, "passed": True,
            })
    return results, checks, before


def configure_style():
    font_manager.findfont("Arial", fallback_to_default=False)
    matplotlib.rcParams.update({
        "font.family": "Arial", "font.size": 6.3,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "axes.linewidth": 0.85, "axes.grid": False,
        "axes.labelsize": 7.7, "xtick.labelsize": 6.3, "ytick.labelsize": 6.3,
        "xtick.major.size": 2.1, "ytick.major.size": 2.1,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.pad": 1.4, "ytick.major.pad": 1.4,
        "legend.frameon": False, "figure.facecolor": "white",
        "savefig.facecolor": "white", "svg.hashsalt": "figure5-new-prediction",
    })


def delta_label(delta):
    return f"{delta:+.3f}" if abs(delta) < 0.01 else f"{delta:+.2f}"


def draw_panel(ax, panel, results):
    rows = [results[(panel["task"], feature)] for feature, _, _ in FEATURES]
    for row in rows:
        ax.plot(row["fpr"], row["tpr"], color=row["color"], linewidth=1.45,
                label=f"{row['model']} (AUC={row['auc']:.3f})", zorder=3)
    ax.plot([0, 1], [0, 1], linestyle=(0, (2.2, 1.5)),
            color="#B8B8B8", linewidth=1.0, zorder=1)
    ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal",
           xlabel="False positive rate", ylabel="True positive rate")
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.xaxis.labelpad = 2
    ax.yaxis.labelpad = 2
    ax.set_title(panel["title"], fontsize=10.5, loc="left", pad=7.2)
    ax.spines[["top", "right"]].set_visible(False)
    legend = ax.legend(loc="lower right", fontsize=5.05, handlelength=0.85,
                       handletextpad=0.28, labelspacing=0.14,
                       borderpad=0, borderaxespad=0.58)
    for line in legend.get_lines():
        line.set_linewidth(1.9)

    inset = ax.inset_axes([0.49, 0.24, 0.43, 0.385])
    aucs = [row["auc"] for row in rows]
    baseline = aucs[0]
    # A visibly truncated AUC axis matches the existing inset design.
    inset.bar([0, 1, 2], aucs, width=0.84,
              color=[row["color"] for row in rows], alpha=0.76, linewidth=0)
    inset.axhline(baseline, color="#7A7A7A", linestyle=(0, (3, 2)), linewidth=0.7)
    inset.set_xlim(-0.75, 2.8)
    inset.set_ylim(*panel["ylim"])
    inset.set_yticks(panel["yticks"])
    inset.set_xticks([])
    inset.tick_params(axis="y", labelsize=6.1, pad=0.5, direction="in", length=1.7)
    inset.yaxis.tick_right()
    inset.spines[["top", "left"]].set_visible(False)
    inset.spines["right"].set_visible(True)
    inset.text(1.02, 1.07, "AUC", transform=inset.transAxes,
               fontsize=6.5, ha="center", va="bottom")
    for x, auc in enumerate(aucs[1:], 1):
        text_x = x + (-0.06 if x == 1 else 0.06)
        inset.text(text_x, auc + (panel["ylim"][1] - panel["ylim"][0]) * 0.035,
                   delta_label(auc - baseline), fontsize=6.0, ha="center", va="bottom")
    return legend, inset


def check_layout(fig, axes, legends):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax, legend in zip(axes, legends):
        if not ax.get_window_extent(renderer).contains(*legend.get_window_extent(renderer).get_points()[0]):
            raise ValueError("Legend begins outside its ROC panel")
        if not ax.get_window_extent(renderer).contains(*legend.get_window_extent(renderer).get_points()[1]):
            raise ValueError("Legend extends outside its ROC panel")
        for inset in ax.child_axes:
            labels = [label.get_window_extent(renderer) for label in inset.texts]
            for index, first in enumerate(labels):
                if any(first.overlaps(second) for second in labels[index + 1:]):
                    raise ValueError("Inset AUC annotations overlap")
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.width and box.height:
            if box.x0 < -0.5 or box.y0 < -0.5 or box.x1 > fig.bbox.x1 + 0.5 or box.y1 > fig.bbox.y1 + 0.5:
                raise ValueError(f"Text outside exported canvas: {text.get_text()}")


def export_figure(fig, stem, dpi, with_tiff=True):
    fig.savefig(stem.with_suffix(".png"), dpi=dpi)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".svg"))
    if with_tiff:
        fig.savefig(stem.with_suffix(".tiff"), dpi=dpi,
                    pil_kwargs={"compression": "tiff_lzw"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()
    if args.dpi < 300:
        raise ValueError("Publication raster exports require at least 300 dpi")
    results, checks, before = load_and_validate(args.input_dir)
    configure_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    width = 183 / 25.4
    height = width * 625 / 1969
    fig = plt.figure(figsize=(width, height), dpi=180)
    axes, legends = [], []
    for panel, left, label_left in zip(PANELS, [0.0508, 0.4078, 0.741], [0.0076, 0.373, 0.706]):
        ax = fig.add_axes([left, 0.1248, 0.243, 0.7664])
        legend, _ = draw_panel(ax, panel, results)
        fig.text(label_left, 0.985, panel["letter"], ha="left", va="top",
                 fontsize=10.8, fontweight="bold", color="#241B18")
        axes.append(ax)
        legends.append(legend)
    check_layout(fig, axes, legends)
    export_figure(fig, args.output_dir / "Figure5_new_prediction", args.dpi)
    fig.savefig(args.output_dir / "Figure5_preview.png", dpi=280)
    plt.close(fig)

    for panel in PANELS:
        single = plt.figure(figsize=(2.65, 2.65), dpi=180)
        ax = single.add_axes([0.14, 0.14, 0.73, 0.73])
        legend, _ = draw_panel(ax, panel, results)
        ax.set_title(panel["title"], fontsize=9.2, loc="left", pad=7.2)
        single.text(0.01, 0.986, panel["letter"], fontsize=10.8,
                    fontweight="bold", ha="left", va="top", color="#241B18")
        check_layout(single, [ax], [legend])
        export_figure(single, args.output_dir / f"Figure5{panel['letter']}", args.dpi, with_tiff=False)
        plt.close(single)

    metrics_rows, curve_frames = [], []
    for panel in PANELS:
        baseline = results[(panel["task"], "covariates_only")]["auc"]
        for feature, _, _ in FEATURES:
            row = results[(panel["task"], feature)]
            metrics_rows.append({key: value for key, value in row.items()
                                 if key not in {"fpr", "tpr", "color"}} |
                                {"delta_auc_vs_covariates": row["auc"] - baseline})
            curve_frames.append(pd.DataFrame({
                "panel": panel["letter"], "task": panel["task"], "feature_set": feature,
                "false_positive_rate": row["fpr"], "true_positive_rate": row["tpr"],
            }))
    pd.DataFrame(metrics_rows).to_csv(args.output_dir / "Figure5_metrics.csv", index=False)
    pd.concat(curve_frames, ignore_index=True).to_csv(args.output_dir / "Figure5_ROC_source_data.csv", index=False)
    pd.DataFrame(checks).to_csv(args.output_dir / "Figure5_validation_checks.csv", index=False)
    (args.output_dir / "Figure5_contract.json").write_text(json.dumps(CONTRACT, indent=2), encoding="utf-8")

    caption = (
        "Figure 5 | Prediction of incident Parkinson's disease using the pain-related "
        "370-protein panel. Receiver operating characteristic (ROC) curves compare "
        "LightGBM models using covariates only (gray), the 370-protein panel only "
        "(teal), and covariates plus the panel (orange). A, Full follow-up "
        "(52,096 participants; 802 events). B, PD diagnosed after 2 and within 5 "
        "years (51,251 participants; 292 events). C, PD diagnosed after 5 years "
        "(50,959 participants; 510 events), without an upper 10-year limit. "
        "Prevalent PD and PD diagnosed within 2 years were excluded. Curves and "
        "AUCs use pooled out-of-fold predicted risks from stratified 10-fold "
        "cross-validation. Insets show AUCs on truncated, labelled axes; "
        "annotations indicate absolute AUC differences from the covariates-only "
        "model. The diagonal dashed line denotes chance discrimination. "
        "The protein panel was selected before PD cross-validation; these "
        "results represent conditional internal evaluation, not external "
        "validation. The protein-panel models are distinct from models using "
        "a single scalar proteomic pain score. No confidence intervals or "
        "formal between-model significance tests are displayed.\n"
    )
    (args.output_dir / "Figure5_legend.txt").write_text(caption, encoding="utf-8")
    if any(digest(path) != expected for path, expected in before.items()):
        raise ValueError("Prediction input files changed during plotting")
    audit = {
        "all_nine_saved_aucs_reproduced": True, "inputs_unchanged": True,
        "cohort_outcome_and_fold_match_within_each_window": True,
        "all_text_within_canvas": True, "all_legends_within_panels": True,
        "input_hashes": {path.name: value for path, value in before.items()},
        "python": platform.python_version(), "matplotlib": matplotlib.__version__,
        "numpy": np.__version__, "pandas": pd.__version__, "sklearn": sklearn.__version__,
        "width_mm": 183, "height_mm": height * 25.4, "raster_dpi": args.dpi,
        "only_aggregate_metrics_and_curves_exported": True,
    }
    (args.output_dir / "Figure5_QA.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(pd.DataFrame(metrics_rows).to_string(index=False))
    print(f"\nValidated figure exports: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
