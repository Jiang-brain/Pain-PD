"""Run original-stream R bootstrap jobs with bounded local process concurrency."""

import argparse
import os
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", required=True)
    parser.add_argument("--reference-dir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 12:
        raise ValueError("Use 1-12 workers to bound resource consumption")
    code = Path(__file__).resolve().parents[2]
    out = Path(os.environ.get("PD_PAIN_RESULTS_DIR", code / "results")) / "sensitivity_nonproteomic_weights/mediation24"
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PD_PAIN_PREDICTION_SOURCE_DIR=str(args.reference_dir),
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1")
    flags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    jobs = out / "bootstrap_jobs.rds"
    if not jobs.exists():
        subprocess.run([args.rscript, "04_proteomics_mediation/03_shared_protein_mediation/02_bootstrap_mediation.R", "--prepare-jobs"],
                       cwd=code, env=env, creationflags=flags, check=True)
    processes, logs = [], []
    try:
        for worker in range(1, args.workers + 1):
            log = (out / f"bootstrap_worker_{worker:02d}.log").open("a", encoding="utf-8")
            logs.append(log)
            process = subprocess.Popen(
                [args.rscript, "04_proteomics_mediation/03_shared_protein_mediation/02_bootstrap_mediation.R", str(worker), str(args.workers)],
                cwd=code, env=env, stdout=log, stderr=subprocess.STDOUT, creationflags=flags,
            )
            processes.append(process)
        started = time.monotonic()
        while any(process.poll() is None for process in processes):
            failed = [process.returncode for process in processes if process.returncode not in (None, 0)]
            if failed:
                raise RuntimeError(f"Bootstrap worker failed with status {failed}; inspect worker logs")
            paths = list((out / "bootstrap_checkpoints").glob("*.rds"))
            print(f"Mediation bootstrap: {len(paths) * 20}/25000 attempts checkpointed; "
                  f"elapsed {(time.monotonic() - started) / 60:.1f} min", flush=True)
            time.sleep(30)
        if any(process.returncode != 0 for process in processes):
            raise RuntimeError("At least one bootstrap worker failed")
        print("All 25,000 requested bootstrap attempts completed; collecting results.", flush=True)
        subprocess.run([args.rscript, "04_proteomics_mediation/03_shared_protein_mediation/03_collect_and_report.R"],
                       cwd=code, env=env, creationflags=flags, check=True)
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
                process.wait()
        for log in logs:
            log.close()


if __name__ == "__main__":
    main()
