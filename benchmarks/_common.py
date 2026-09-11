from __future__ import annotations

import json
import os
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


def _cpu_model() -> str | None:
    try:
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return platform.processor() or None


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
    return {"git_commit": commit, "git_dirty": None if status is None else bool(status)}


def _cpu_governors() -> list[str]:
    governors = set()
    base = Path("/sys/devices/system/cpu")
    try:
        for path in base.glob("cpu[0-9]*/cpufreq/scaling_governor"):
            try:
                governors.add(path.read_text().strip())
            except OSError:
                pass
    except OSError:
        pass
    return sorted(g for g in governors if g)


def _cpu_affinity() -> list[int] | None:
    try:
        return sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        return None


def environment_metadata() -> dict:
    thread_vars = {
        name: os.environ.get(name)
        for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")
    }
    return {
        "petal2d_version": getattr(petal2d, "__version__", "unknown"),
        "python": sys.version,
        "platform": platform.platform(),
        "cpu_model": _cpu_model(),
        "logical_cpu_count": os.cpu_count(),
        "cpu_governors": _cpu_governors(),
        "cpu_affinity": _cpu_affinity(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "thread_environment": thread_vars,
        **git_metadata(),
    }


def make_payload(study_name: str, quick: bool, studies: dict, methodology: dict) -> dict:
    return {
        "schema_version": 3,
        "study": study_name,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "quick": bool(quick),
        "environment": environment_metadata(),
        "methodology": methodology,
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


def summarize_samples(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    q25, q75 = np.percentile(arr, [25.0, 75.0])
    return {
        "n": int(arr.size),
        "median": float(np.median(arr)),
        "iqr": float(q75 - q25),
        "q25": float(q25),
        "q75": float(q75),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }
