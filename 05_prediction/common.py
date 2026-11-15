"""Paths and read-only source loading for the controlled prediction sensitivity."""
import hashlib
import importlib.util
import os
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parents[1]
BRANCH = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("PD_PAIN_DATA_DIR", CODE_ROOT / "data"))
REFERENCE_DIR = Path(os.environ.get("PD_PAIN_PREDICTION_SOURCE_DIR", DATA_DIR / "prediction_reference"))
RESULTS_ROOT = Path(os.environ.get("PD_PAIN_RESULTS_DIR", CODE_ROOT / "results"))
SENSITIVITY_DIR = RESULTS_ROOT / "sensitivity_nonproteomic_weights"
OUT = SENSITIVITY_DIR / "prediction"
OUT.mkdir(parents=True, exist_ok=True)
NEW_SCORE_PATH = SENSITIVITY_DIR / "composite_pain_scores.csv"


def load_reference(filename, module_name):
    path = BRANCH / "source_reference" / filename
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def md5(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "md5").hexdigest()
