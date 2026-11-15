"""Final Figure 4 panel entry point; plot saved results without refitting models."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


CODE = Path(__file__).resolve().parents[1]
BRANCH = CODE / "04_proteomics_mediation/03_shared_protein_mediation"


def build_commands(args):
    commands = []
    if "associations" in args.panels or args.path_diagram:
        command = [sys.executable, str(BRANCH / "05_plot_figure4_panels.py"),
                   "--input-dir", str(args.input_dir), "--output-dir", str(args.output_dir)]
        if not args.path_diagram:
            command.append("--associations-only")
        commands.append(command)
    if "mediation" in args.panels:
        commands.append([sys.executable, str(BRANCH / "08_plot_figure4E_forest_compact.py"),
                         "--input", str(args.input_dir / "mediation24"),
                         "--output", str(args.output_dir / "E_right_compact_CI")])
    if "enrichment" in args.panels:
        commands.append([sys.executable, str(BRANCH / "07_plot_figure4G_top5.py"),
                         "--input", str(args.enrichment_dir),
                         "--output", str(args.output_dir / "G_top5_24")])
    return commands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    results = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE / "results"))
    parser.add_argument("--input-dir", type=Path, default=results / "sensitivity_nonproteomic_weights")
    parser.add_argument("--output-dir", type=Path, default=CODE.parent / "figures/figure4_independent_weights")
    parser.add_argument("--enrichment-dir", type=Path, default=CODE / "figure_data/figure4/enrichment_24")
    parser.add_argument("--panels", nargs="+", choices=("associations", "mediation", "enrichment"),
                        default=["associations", "mediation", "enrichment"])
    parser.add_argument("--path-diagram", action="store_true",
                        help="Also export the descriptive path diagram/composite E; use the separate compact right forest for final layout.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only; do not write figures or run models.")
    args = parser.parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.enrichment_dir = args.enrichment_dir.resolve()
    for command in build_commands(args):
        print(subprocess.list2cmdline(command), flush=True)
        if not args.dry_run:
            subprocess.run(command, cwd=CODE, check=True)


if __name__ == "__main__":
    main()
