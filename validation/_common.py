from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy

import petal2d


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = Path(__file__).resolve().parent / "results"
FIGURES_DIR = RESULTS_DIR / "figures"


def git_metadata() -> dict:
    def run(*args: str) -> str | None:
        try:
            result = subprocess.run(
                ["git", *args],
                cwd=REPO_ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError):
            return None
        return result.stdout.strip()

    commit = run("rev-parse", "HEAD")
    status = run("status", "--porcelain")
    return {
        "git_commit": commit,
        "git_dirty": None if status is None else bool(status),
    }


def environment_metadata() -> dict:
    return {
        "petal2d_version": getattr(petal2d, "__version__", "unknown"),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        **git_metadata(),
    }


def make_payload(study_name: str, quick: bool, studies: dict) -> dict:
    return {
        "schema_version": 1,
        "study": study_name,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "quick": bool(quick),
        "environment": environment_metadata(),
        "studies": studies,
    }


def write_json(payload: dict, output: str | Path) -> Path:
    path = Path(output)
    if not path.is_absolute():
        path = REPO_ROOT / path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def ensure_figures_dir(path: str | Path | None = None) -> Path:
    out = FIGURES_DIR if path is None else Path(path)
    if not out.is_absolute():
        out = REPO_ROOT / out
    out.mkdir(parents=True, exist_ok=True)
    return out
