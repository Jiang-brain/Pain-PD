"""Compact Figure 4E right forest plot, using completed 24-protein bootstrap CIs."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


CODE = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE / "results")) / "sensitivity_nonproteomic_weights/mediation24"
DEFAULT_OUTPUT = CODE.parent / "figures/figure4_independent_weights/E_right_compact_CI"
BLUE = "#2B9DC6"
RED = "#D54843"
MARKER_AREA = 46
JOINT_MARKER_AREA = 56


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def excludes_zero(lower, upper):
    return lower > 0 or upper < 0


def configure():
    font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7,
        "axes.labelsize": 8, "axes.linewidth": 0.65,
        "xtick.labelsize": 6.8, "ytick.labelsize": 7,
        "xtick.major.width": 0.6, "xtick.major.size": 3,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })


def source_rows(single, joint):
    records = []
    for r in single + [joint]:
        group = "joint" if r is joint else "single"
        low = float(r["proportion_lower"]) * 100
        high = float(r["proportion_upper"]) * 100
        estimate = float(r["proportion_mediated"]) * 100
        indirect_low = float(r["indirect_lower"])
        indirect_high = float(r["indirect_upper"])
        if not all(math.isfinite(x) for x in [low,high,estimate,indirect_low,indirect_high]):
            raise ValueError("Nonfinite estimate or interval")
        if not low <= estimate <= high:
            raise ValueError("Point estimate is outside its percentile interval")
        if not (int(r["n"]) == 46591 and int(r["events"]) == 707 and int(r["n_boot_success"]) == 1000):
            raise ValueError("Unexpected cohort or bootstrap count")
        indirect_flag = excludes_zero(indirect_low, indirect_high)
        ratio_flag = excludes_zero(low, high)
        if indirect_flag != ratio_flag:
            raise ValueError("Indirect-effect and proportion CI flags differ; inspect before reporting")
        records.append({
            "protein_or_set": "Joint 24 proteins" if group == "joint" else r["protein"],
            "model": group, "n": int(r["n"]), "events": int(r["events"]),
            "proportion_mediated_percent": estimate,
            "proportion_lower_percent": low, "proportion_upper_percent": high,
            "indirect_logHR": float(r["indirect_logHR"]),
            "indirect_logHR_lower": indirect_low, "indirect_logHR_upper": indirect_high,
            "indirect_95CI_excludes_zero": indirect_flag,
            "proportion_95CI_excludes_zero": ratio_flag,
            "n_bootstrap": int(r["n_boot_success"]),
        })
    return records


def draw(records):
    fig = plt.figure(figsize=(3.2, 3.45), dpi=180)
    ax = fig.add_axes([0.345, 0.14, 0.615, 0.83])
    single = records[:-1]
    joint = records[-1]
    yy = list(range(len(single)))
    joint_y = 24.2
    ax.axvline(0, color="#777777", lw=0.65, linestyle=(0,(4,3)), zorder=0)
    ax.hlines(yy, [r["proportion_lower_percent"] for r in single],
              [r["proportion_upper_percent"] for r in single],
              color="#666666", lw=0.65, zorder=1)
    ax.scatter([r["proportion_mediated_percent"] for r in single], yy,
               s=MARKER_AREA, marker="o", facecolors=BLUE, edgecolors="#333333",
               linewidths=0.4, zorder=3)
    ax.hlines(joint_y, joint["proportion_lower_percent"], joint["proportion_upper_percent"],
              color="#E35B58", lw=0.85, zorder=2)
    ax.scatter(joint["proportion_mediated_percent"], joint_y,
               s=JOINT_MARKER_AREA, marker="o", facecolors=RED,
               edgecolors="#333333", linewidths=0.4, zorder=4)
    ax.set_yticks(yy + [joint_y], labels=[r["protein_or_set"] for r in records])
    for tick in ax.get_yticklabels():
        tick.set_fontstyle("italic")
    ax.get_yticklabels()[-1].set_color("#D72F2F")
    ax.get_yticklabels()[-1].set_fontweight("bold")
    ax.set_ylim(25.0, -0.6)
    ax.set_xlim(-5, 35)
    ticks = list(range(-5,36,5))
    ax.set_xticks(ticks, labels=[f"{v}%" for v in ticks])
    ax.set_xlabel("Proportion mediated (95% CI)", labelpad=4)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.spines[["top","right"]].set_visible(False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.width and box.height and (box.x0 < -0.5 or box.y0 < -0.5 or
           box.x1 > fig.bbox.x1 + 0.5 or box.y1 > fig.bbox.y1 + 0.5):
            raise ValueError(f"Text outside canvas: {text.get_text()}")
    label_boxes = [t.get_window_extent(renderer) for t in ax.get_yticklabels()]
    if any(a.y0 < b.y1 + 0.5 for a,b in zip(label_boxes,label_boxes[1:])):
        raise ValueError("Protein labels overlap vertically")
    pixels_per_row = abs(ax.transData.transform((0,1))[1] - ax.transData.transform((0,0))[1])
    pitch_pt = pixels_per_row * 72 / fig.dpi
    marker_pt = math.sqrt(MARKER_AREA)
    if pitch_pt <= marker_pt + 0.4:
        raise ValueError("Adjacent circles overlap")
    if min(r["proportion_lower_percent"] for r in records) < -5 or max(r["proportion_upper_percent"] for r in records) > 35:
        raise ValueError("CI extends outside x limits")
    return fig, {"row_pitch_pt": pitch_pt, "single_marker_diameter_pt": marker_pt,
                 "gap_between_adjacent_circles_pt": pitch_pt - marker_pt,
                 "labels_and_intervals_inside_canvas": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    files = [args.input / "per_protein_mediation_results.csv",
             args.input / "joint24_mediation_results.csv"]
    hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    single, joint = read(files[0]), read(files[1])[0]
    if len(single) != 24 or len({r["protein"] for r in single}) != 24 or any(r["protein"]=="COMP" for r in single):
        raise ValueError("Expected 24 unique shared proteins excluding COMP")
    single.sort(key=lambda r: float(r["proportion_mediated"]))
    records = source_rows(single, joint)
    ci_count = sum(r["indirect_95CI_excludes_zero"] for r in records if r["model"] == "single")
    configure()
    args.output.mkdir(parents=True, exist_ok=True)
    fig, qa = draw(records)
    stem = args.output / "Figure4E_right_shared24_large_circles"
    exported = []
    for suffix in ("png","pdf","svg","tiff"):
        options = {"dpi": 600}
        if suffix == "png": options["pil_kwargs"] = {"compress_level": 9}
        if suffix == "tiff": options["pil_kwargs"] = {"compression": "tiff_lzw"}
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, **options); exported.append(path)
    preview = stem.with_name(stem.name + "_preview").with_suffix(".png")
    fig.savefig(preview,dpi=240)
    plt.close(fig)
    source_csv = args.output / "Figure4E_right_source_data.csv"
    with source_csv.open("w",encoding="utf-8-sig",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    caption = args.output / "Figure4E_right_legend.txt"
    caption.write_text(
        "Figure 4E, right panel. Observational IORW mediation estimates for the independently weighted pain score, "
        "based on 46,591 participants with 707 incident PD events. Blue circles show the 24 single-protein mediation "
        "proportions; the red circle shows the joint 24-protein estimate. Horizontal lines represent percentile 95% "
        "confidence intervals from 1,000 participant-level bootstrap resamples per model. The dashed line marks zero. "
        f"All single proteins have the same symbol; color does not encode significance. For {ci_count} of 24 proteins, the indirect-effect "
        "95% bootstrap CIs exclude zero; the proportion CIs give the same classification in these data. "
        "Intervals are not adjusted for multiple comparisons, and no one-sided or FDR criterion is used for this figure's "
        "interpretation. The joint proportion is 17.82% (95% CI 9.65%-32.32%). "
        "These observational estimates do not establish causal mediation.\n",encoding="utf-8")
    audit={"core_conclusion":"Display single and joint mediation proportions and uncertainty without causal inference.",
           "single_models":24,"single_indirect_CI_excludes_zero":ci_count,"n":46591,"events":707,
           "bootstrap_resamples_per_model":1000,"x_limits_percent":[-5,35],
           "width_mm":3.2*25.4,"height_mm":3.45*25.4,"dpi":600,
           "joint_percent":records[-1]["proportion_mediated_percent"],
           "joint_CI_percent":[records[-1]["proportion_lower_percent"],records[-1]["proportion_upper_percent"]],
           "source_sha256":{p.name:h for p,h in hashes.items()},"layout":qa}
    assert all(hashlib.sha256(p.read_bytes()).hexdigest()==h for p,h in hashes.items())
    audit_path = args.output / "Figure4E_right_QA.json"
    audit_path.write_text(json.dumps(audit,indent=2),encoding="utf-8")
    archive=args.output / "Figure4E_right_shared24_package.zip"
    with ZipFile(archive,"w",ZIP_DEFLATED,compresslevel=9) as z:
        for path in exported+[source_csv,caption,preview]:
            z.write(path,arcname=path.name)
    print(json.dumps({"output":str(stem),"package_bytes":archive.stat().st_size,
                      "layout":qa,"CI_excludes_zero":ci_count},indent=2))


if __name__=="__main__":
    main()
