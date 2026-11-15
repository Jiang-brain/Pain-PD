"""Final Figure 5 entry point: 370-protein LightGBM OOF results, not model refits."""

import argparse
import os
from pathlib import Path
import subprocess
import sys


CODE = Path(__file__).resolve().parents[1]


def build_command(args):
    return [sys.executable, str(CODE / "05_prediction/07_plot_figure5.py"),
            "--input-dir", str(args.input_dir), "--output-dir", str(args.output_dir),
            "--dpi", str(args.dpi)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    results = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE / "results"))
    parser.add_argument("--input-dir", type=Path, default=results / "sensitivity_nonproteomic_weights/prediction")
    parser.add_argument("--output-dir", type=Path, default=CODE.parent / "figures/figure5_new_prediction")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--dry-run", action="store_true", help="Print command only; no model or figure is changed.")
    args = parser.parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    command = build_command(args)
    print(subprocess.list2cmdline(command), flush=True)
    if not args.dry_run:
        subprocess.run(command, cwd=CODE, check=True)


if __name__ == "__main__":
    main()
