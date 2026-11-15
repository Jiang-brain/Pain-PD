"""Redraw the two Figure 4G panels using the manuscript's per-category top-five rule."""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch


CODE = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = CODE / "figure_data/figure4/enrichment_24"
DEFAULT_OUTPUT = CODE.parent / "figures/figure4_independent_weights/G_top5_24"

COLORS = [
    {"bar": "#D1E2E9", "legend": "#4D8BA3", "id": "#39788F"},
    {"bar": "#F4D5CE", "legend": "#EF9D8B", "id": "#B86E5D"},
    {"bar": "#CCDDD6", "legend": "#93AFA5", "id": "#658C80"},
]
GO_GROUPS = [("Process", "GO Biological Process"),
             ("Function", "GO Molecular Function"),
             ("Component", "GO Cellular Component")]
PATHWAY_GROUPS = [("KEGG", "KEGG"), ("RCTM", "Reactome"),
                  ("WikiPathways", "WikiPathways")]
WIDTH, HEIGHT = 3.6, 2.65
TIE_POLICY = "Ascending exported FDR; ties retain the supplied export row order."


def configure():
    font_manager.findfont("Arial", fallback_to_default=False)
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 7,
        "axes.labelsize": 8, "axes.linewidth": 0.6,
        "xtick.labelsize": 6.8, "ytick.labelsize": 6.4,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 0,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "figure.facecolor": "white", "savefig.facecolor": "white",
        "mathtext.fontset": "custom", "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    })


def select_terms(input_dir, groups, shared):
    selected, manifest = [], []
    for group_index, (suffix, category) in enumerate(groups):
        source = input_dir / f"enrichment.{suffix}.tsv"
        with source.open(encoding="utf-8-sig", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        if len({r["#term ID"] for r in records}) != len(records):
            raise ValueError(f"Duplicate IDs within {category}")
        significant = []
        for export_order, r in enumerate(records, start=1):
            fdr = float(r["false discovery rate"])
            if not math.isfinite(fdr) or not 0 < fdr <= 1:
                raise ValueError(f"Invalid FDR for {r['#term ID']}")
            if fdr < 0.05:
                significant.append((fdr, export_order, r))
        significant.sort(key=lambda r: (r[0], r[1]))
        if len(significant) < 5:
            raise ValueError(f"Fewer than five FDR-significant terms in {category}")
        for rank, (fdr, export_order, r) in enumerate(significant[:5], start=1):
            matching = [s.strip() for s in r["matching proteins in your network (labels)"].split(",")]
            if not set(matching).issubset(shared):
                raise ValueError(f"Non-shared protein in {r['#term ID']}")
            if len(set(matching)) != int(r["observed gene count"]):
                raise ValueError(f"Observed/matching count mismatch: {r['#term ID']}")
            selected.append({
                "category": category, "category_index": group_index,
                "rank_within_category": rank,
                "term_id": r["#term ID"], "description": r["term description"],
                "fdr": fdr, "minus_log10_fdr": -math.log10(fdr),
                "observed_gene_count": int(r["observed gene count"]),
                "background_gene_count": int(r["background gene count"]),
                "strength": float(r["strength"]), "signal": float(r["signal"]),
                "matching_proteins": r["matching proteins in your network (labels)"],
                "source_export_row": export_order,
            })
        manifest.append({"file": source.name,
                         "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
                         "supplied_terms": len(records), "significant_terms": len(significant)})
    assert len(selected) == 15
    return selected, manifest


def draw_panel(fig, box, selected, groups, title):
    x0, y0, width, height = box
    ax = fig.add_axes([x0 + 0.215 * width, y0 + 0.195 * height,
                       0.76 * width, 0.715 * height])
    positions = list(range(len(selected)))
    ax.barh(positions, [r["minus_log10_fdr"] for r in selected], height=0.75,
            color=[COLORS[r["category_index"]]["bar"] for r in selected],
            edgecolor="none", zorder=1)
    ax.set_xlim(0, 10)
    ax.set_xticks([0, 2, 4, 6, 8, 10])
    ax.set_ylim(14.65, -0.65)
    ax.set_yticks(positions, labels=[r["term_id"] for r in selected],
                  fontstyle="italic", fontsize=6.4)
    ax.tick_params(axis="y", pad=3)
    for tick, r in zip(ax.get_yticklabels(), selected):
        tick.set_color(COLORS[r["category_index"]]["id"])
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel(r"$-\log_{10}(\mathrm{FDR})$", labelpad=2)
    labels = []
    for y, row in zip(positions, selected):
        labels.append(ax.text(0.10, y, row["description"], fontsize=6.8,
                              va="center", ha="left", color="black", zorder=3))
    handles = [Patch(facecolor=COLORS[i]["legend"], edgecolor="none", label=name)
               for i, (_, name) in enumerate(groups)]
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(x0 + 0.51 * width, y0 + 0.995 * height),
               ncol=3, frameon=False, fontsize=5.8, handlelength=0.9,
               handleheight=0.7, handletextpad=0.25, columnspacing=0.8,
               borderpad=0, borderaxespad=0)
    fig.text(x0 + 0.59 * width, y0 + 0.016 * height, title,
             fontsize=9, fontweight="bold", ha="center", va="bottom")
    return ax, labels


def qa_layout(fig, panels):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    adjusted = []
    for ax, labels in panels:
        right = ax.get_window_extent(renderer).x1 - 3
        for text in labels:
            while text.get_window_extent(renderer).x1 > right and text.get_fontsize() > 6.2:
                text.set_fontsize(round(text.get_fontsize() - 0.1, 2))
            if text.get_window_extent(renderer).x1 > right:
                raise ValueError(f"Description exceeds plot: {text.get_text()}")
            adjusted.append({"description": text.get_text(), "font_size_pt": text.get_fontsize()})
        boxes = [t.get_window_extent(renderer) for t in labels]
        if any(a.y0 < b.y1 + 0.5 for a, b in zip(boxes, boxes[1:])):
            raise ValueError("Adjacent row descriptions overlap")
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        box = text.get_window_extent(renderer)
        if box.width and box.height and (box.x0 < -0.5 or box.y0 < -0.5 or
            box.x1 > fig.bbox.x1 + 0.5 or box.y1 > fig.bbox.y1 + 0.5):
            raise ValueError(f"Text outside canvas: {text.get_text()}")
    return adjusted


def export_panel(selected, groups, title, output_stem):
    fig = plt.figure(figsize=(WIDTH, HEIGHT), dpi=180)
    panel = draw_panel(fig, (0, 0, 1, 1), selected, groups, title)
    typography = qa_layout(fig, [panel])
    exports = []
    for suffix in ("png", "pdf", "svg", "tiff"):
        path = output_stem.with_suffix(f".{suffix}")
        options = {"dpi": 600}
        if suffix == "tiff":
            options["pil_kwargs"] = {"compression": "tiff_lzw"}
        elif suffix == "png":
            options["pil_kwargs"] = {"compress_level": 9}
        fig.savefig(path, **options)
        exports.append(path)
    preview = output_stem.with_name(output_stem.name + "_preview").with_suffix(".png")
    fig.savefig(preview, dpi=240)
    plt.close(fig)
    return exports, typography


def write_source(path, rows):
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    shared_path = args.input.parent / "shared24_for_PPI_enrichment.txt"
    shared = {s.strip() for s in shared_path.read_text().splitlines() if s.strip()}
    if len(shared) != 24 or "COMP" in shared:
        raise ValueError("Expected the fixed 24-protein input set excluding COMP")
    configure()
    args.output.mkdir(parents=True, exist_ok=True)
    go, go_sources = select_terms(args.input, GO_GROUPS, shared)
    pathways, pathway_sources = select_terms(args.input, PATHWAY_GROUPS, shared)
    go_files, go_fonts = export_panel(go, GO_GROUPS, "GO enrichment",
                                     args.output / "Figure4G_GO_top5_shared24")
    path_files, path_fonts = export_panel(pathways, PATHWAY_GROUPS, "Pathway enrichment",
                                         args.output / "Figure4G_pathways_top5_shared24")
    go_csv = args.output / "Figure4G_GO_top5_source_data.csv"
    path_csv = args.output / "Figure4G_pathways_top5_source_data.csv"
    write_source(go_csv, go); write_source(path_csv, pathways)
    pair = plt.figure(figsize=(2 * WIDTH, HEIGHT), dpi=180)
    panels = [draw_panel(pair, (0, 0, 0.5, 1), go, GO_GROUPS, "GO enrichment"),
              draw_panel(pair, (0.5, 0, 0.5, 1), pathways, PATHWAY_GROUPS, "Pathway enrichment")]
    qa_layout(pair, panels)
    pair.savefig(args.output / "Figure4G_two_panels_preview.png", dpi=240)
    plt.close(pair)
    caption = args.output / "Figure4G_legend.txt"
    caption.write_text(
        "Figure 4G. Functional enrichment of the 24 proteins shared between the independently weighted pain score and incident PD. "
        "The left panel shows the five terms with the smallest exported FDR in each GO category "
        "(biological process, molecular function and cellular component). The right panel shows the five pathways "
        "with the smallest exported FDR in each database (KEGG, Reactome and WikiPathways). "
        "All displayed terms have FDR < 0.05. Bar length denotes -log10(FDR); colors identify categories. "
        + TIE_POLICY + " FDR values and term annotations are used directly from the supplied STRING exports; "
        "enrichment tests and multiplicity correction were not rerun for plotting. "
        "Enrichment is descriptive and does not establish causal mechanisms.\n", encoding="utf-8")
    manifest = go_sources + pathway_sources
    assert all(hashlib.sha256((args.input / s["file"]).read_bytes()).hexdigest() == s["sha256"] for s in manifest)
    audit = {"shared_proteins": 24, "terms_per_category": 5,
             "GO_terms": len(go), "pathway_terms": len(pathways),
             "selection_rule": TIE_POLICY, "x_axis_limits": [0,10],
             "panel_width_mm": WIDTH * 25.4, "panel_height_mm": HEIGHT * 25.4,
             "export_dpi": 600, "layout_checks_passed": True,
             "all_input_files_unchanged": True, "sources": manifest,
             "GO_typography": go_fonts, "pathway_typography": path_fonts}
    (args.output / "Figure4G_selection_and_QA.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    archive = args.output / "Figure4G_top5_shared24_package.zip"
    with ZipFile(archive, "w", ZIP_DEFLATED, compresslevel=9) as z:
        for path in go_files + path_files + [go_csv, path_csv, caption,
                    args.output / "Figure4G_two_panels_preview.png"]:
            z.write(path, arcname=path.name)
    print(json.dumps({"GO_selected": [(r["term_id"], r["fdr"]) for r in go],
                      "pathways_selected": [(r["term_id"], r["fdr"]) for r in pathways],
                      "output_dir": str(args.output), "package_bytes": archive.stat().st_size}, indent=2))


if __name__ == "__main__":
    main()
