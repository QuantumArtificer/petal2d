#!/usr/bin/env python3
"""Fresh-process peak-RSS benchmark for PETAL2D.

Requires the optional benchmark dependency:

    python3 -m pip install -e ".[benchmark]"

Each measurement is performed in a fresh Python process. The worker samples its
own RSS while constructing PETAL2D and reports peak RSS above a post-import
baseline. No claim about hardware memory bandwidth is made.
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import gc
import json
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np

from petal2d import PolarDecomposition

from _common import make_payload, summarize_samples, write_json


def analytic_field(x, y):
    r2 = x * x + y * y
    return np.exp(-0.60 * r2) * (1.0 + 0.03 * (x**3 - 3.0 * x * y**2))


def worker(config: dict) -> dict:
    try:
        import psutil
    except ImportError as exc:
        raise RuntimeError('benchmark_memory.py requires psutil; install with ".[benchmark]"') from exc

    process = psutil.Process(os.getpid())
    post_import_rss = int(process.memory_info().rss)

    # Construct the user-supplied inputs before establishing the PETAL2D
    # baseline. For sampled input, X/Y are benchmark-only construction
    # temporaries: delete them before the baseline so they are not attributed
    # to PETAL2D. The actual field array and coordinate axes remain resident,
    # exactly as they would in a user's process.
    nxy = int(config["nxy"])
    axis = np.linspace(-5.0, 5.0, nxy)
    if config["input_type"] == "sampled":
        X, Y = np.meshgrid(axis, axis, indexing="ij")
        field = analytic_field(X, Y)
        del X, Y
    else:
        field = analytic_field
    gc.collect()

    input_baseline_rss = int(process.memory_info().rss)
    peak = input_baseline_rss
    stop = threading.Event()

    def monitor():
        nonlocal peak
        while not stop.is_set():
            try:
                peak = max(peak, int(process.memory_info().rss))
            except psutil.Error:
                pass
            time.sleep(0.001)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    try:
        dec = PolarDecomposition(
            field,
            axis,
            axis,
            Nr=int(config["nr"]),
            Ntheta=int(config["ntheta"]),
            rmax=5.0,
            origin=(0.0, 0.0),
            recon_err_tol=1.0,
            interp_method="cubic",
        )
        # Touch key retained arrays so lazy allocations cannot escape the
        # monitored region.
        _ = (
            dec.f_polar.nbytes,
            dec.rho_all.nbytes,
            dec.f_recon.nbytes,
        )
        final_rss = int(process.memory_info().rss)
        peak = max(peak, final_rss)
    finally:
        stop.set()
        thread.join(timeout=0.2)

    return {
        **config,
        "cartesian_points": int(config["nxy"] ** 2),
        "polar_points": int(config["nr"] * config["ntheta"]),
        "post_import_rss_bytes": post_import_rss,
        "input_baseline_rss_bytes": input_baseline_rss,
        "input_resident_rss_bytes": max(0, input_baseline_rss - post_import_rss),
        "peak_rss_bytes": peak,
        "petal2d_peak_incremental_rss_bytes": max(0, peak - input_baseline_rss),
        "final_rss_bytes": final_rss,
        "petal2d_final_incremental_rss_bytes": max(0, final_rss - input_baseline_rss),
    }


def run_worker_subprocess(config: dict) -> dict:
    command = [sys.executable, str(Path(__file__).resolve()), "--worker-json", json.dumps(config)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("Memory worker returned no JSON output.")
    return json.loads(lines[-1])


def aggregate_case(config: dict, repeats: int) -> dict:
    raw = [run_worker_subprocess(config) for _ in range(repeats)]
    return {
        **config,
        "cartesian_points": int(config["nxy"] ** 2),
        "polar_points": int(config["nr"] * config["ntheta"]),
        "post_import_rss_bytes": summarize_samples(
            [r["post_import_rss_bytes"] for r in raw]
        ),
        "input_baseline_rss_bytes": summarize_samples(
            [r["input_baseline_rss_bytes"] for r in raw]
        ),
        "input_resident_rss_bytes": summarize_samples(
            [r["input_resident_rss_bytes"] for r in raw]
        ),
        "peak_rss_bytes": summarize_samples([r["peak_rss_bytes"] for r in raw]),
        "petal2d_peak_incremental_rss_bytes": summarize_samples(
            [r["petal2d_peak_incremental_rss_bytes"] for r in raw]
        ),
        "petal2d_final_incremental_rss_bytes": summarize_samples(
            [r["petal2d_final_incremental_rss_bytes"] for r in raw]
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--output", default="benchmarks/results/memory_benchmark.json")
    parser.add_argument("--worker-json", default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        "--cpu",
        type=int,
        default=None,
        help="Optionally pin the benchmark parent/worker processes to one logical CPU on Linux.",
    )
    args = parser.parse_args()

    if args.cpu is not None:
        if args.cpu < 0:
            raise ValueError("--cpu must be >= 0")
        try:
            os.sched_setaffinity(0, {args.cpu})
        except AttributeError as exc:
            raise RuntimeError("--cpu affinity is supported only where os.sched_setaffinity is available") from exc
        except OSError as exc:
            raise RuntimeError(f"Could not pin benchmark to CPU {args.cpu}: {exc}") from exc

    if args.worker_json is not None:
        print(json.dumps(worker(json.loads(args.worker_json))))
        return

    repeats = args.repeats if args.repeats is not None else (3 if args.quick else 7)
    if repeats < 1:
        raise ValueError("repeats must be >=1")

    polar_cases = (
        [(128, 256), (256, 512), (512, 1024)]
        if args.quick
        else [(64, 128), (128, 256), (256, 512), (512, 1024), (768, 1536), (1024, 2048)]
    )
    nxy_values = [256, 512, 1024] if args.quick else [256, 384, 512, 768, 1024, 1536, 2048]
    strategy_nxy_values = [256] if args.quick else [128, 256, 512]

    studies = {
        "callable_polar_grid": [],
        "sampled_cartesian_grid": [],
        "sampled_radial_resolution_strategy_cost": [],
    }
    for nr, ntheta in polar_cases:
        config = {"input_type": "callable", "nxy": 129, "nr": nr, "ntheta": ntheta}
        print(f"Memory case: {config}")
        studies["callable_polar_grid"].append(aggregate_case(config, repeats))

    for nxy in nxy_values:
        config = {"input_type": "sampled", "nxy": nxy, "nr": 256, "ntheta": 512}
        print(f"Memory case: {config}")
        studies["sampled_cartesian_grid"].append(aggregate_case(config, repeats))

    for nxy in strategy_nxy_values:
        for strategy, nr in (
            ("half_cartesian", max(2, nxy // 2)),
            ("match_cartesian", nxy),
            ("double_cartesian", 2 * nxy),
        ):
            config = {
                "input_type": "sampled",
                "nxy": nxy,
                "nr": nr,
                "ntheta": 512,
                "nr_strategy": strategy,
            }
            print(f"Memory case: {config}")
            studies["sampled_radial_resolution_strategy_cost"].append(
                aggregate_case(config, repeats)
            )

    methodology = {
        "measurement": "Fresh Python process per repetition; worker samples its own RSS every ~1 ms using psutil.",
        "input_baseline": "User inputs are constructed first. For sampled input, benchmark-only X/Y construction meshes are deleted and garbage-collected. The baseline then includes the actual field array plus x/y axes, but excludes PETAL2D construction.",
        "reported_incremental_peak": "petal2d_peak_incremental_rss_bytes = peak RSS during PolarDecomposition construction minus the input-resident baseline RSS in the same worker process.",
        "repeats_per_case": repeats,
        "sampled_cartesian_grid_range": nxy_values,
        "sampled_cartesian_fit_policy": "No affine slope is inferred for sampled-input memory unless the measured asymptotic data visibly support one; plots report median and interquartile range only.",
        "cpu_affinity_request": args.cpu,
        "bandwidth_claim": "None. Peak RSS is a memory-footprint measurement, not a hardware memory-bandwidth measurement.",
    }
    path = write_json(make_payload("memory_benchmark", args.quick, studies, methodology), args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
