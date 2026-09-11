#!/usr/bin/env python3
"""Reproducible PETAL2D runtime benchmark.

The script writes JSON only. It reports medians and interquartile ranges after
warm-up runs. Callable and sampled input are benchmarked separately because the
sampled path includes Cartesian-to-polar interpolation.
"""

from __future__ import annotations

import os

# Set before NumPy/SciPy import for reproducibility on BLAS/OpenMP builds.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import time

import numpy as np

from petal2d import PolarDecomposition

from _common import make_payload, summarize_samples, write_json


class TimedPolarDecomposition(PolarDecomposition):
    """Benchmark-only subclass that times the two major internal stages."""

    def _prepare_field(self):
        start = time.perf_counter()
        super()._prepare_field()
        self.benchmark_prepare_field_seconds = time.perf_counter() - start

    def _decompose(self):
        start = time.perf_counter()
        super()._decompose()
        self.benchmark_decompose_seconds = time.perf_counter() - start


def analytic_field(x, y):
    r2 = x * x + y * y
    c3 = x**3 - 3.0 * x * y**2
    return np.exp(-0.60 * r2) * (1.0 + 0.03 * c3)


def make_input(input_type: str, nxy: int):
    axis = np.linspace(-5.0, 5.0, nxy)
    if input_type == "callable":
        return analytic_field, axis, axis
    X, Y = np.meshgrid(axis, axis, indexing="ij")
    return analytic_field(X, Y), axis, axis


def one_run(*, input_type: str, nxy: int, nr: int, ntheta: int, interp_method: str, nr_strategy: str | None = None) -> dict:
    field, x, y = make_input(input_type, nxy)
    start = time.perf_counter()
    dec = TimedPolarDecomposition(
        field,
        x,
        y,
        Nr=nr,
        Ntheta=ntheta,
        rmax=5.0,
        origin=(0.0, 0.0),
        recon_err_tol=1.0,
        interp_method=interp_method,
    )
    total = time.perf_counter() - start
    return {
        "total": total,
        "prepare_field": float(dec.benchmark_prepare_field_seconds),
        "decompose": float(dec.benchmark_decompose_seconds),
    }


def benchmark_case(config: dict, repeats: int, warmups: int) -> dict:
    for _ in range(warmups):
        one_run(**config)

    samples = {"total": [], "prepare_field": [], "decompose": [], "other_overhead": []}
    for _ in range(repeats):
        result = one_run(**config)
        result["other_overhead"] = max(0.0, result["total"] - result["prepare_field"] - result["decompose"])
        for key in samples:
            samples[key].append(float(result[key]))

    return {
        **config,
        "cartesian_points": int(config["nxy"] ** 2),
        "polar_points": int(config["nr"] * config["ntheta"]),
        "timing_seconds": {key: summarize_samples(values) for key, values in samples.items()},
    }


def run_study(configs: list[dict], repeats: int, warmups: int) -> list[dict]:
    rows = []
    for i, config in enumerate(configs, start=1):
        print(f"[{i}/{len(configs)}] {config}")
        rows.append(benchmark_case(config, repeats, warmups))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--output", default="benchmarks/results/runtime_benchmark.json")
    parser.add_argument(
        "--cpu",
        type=int,
        default=None,
        help="Optionally pin the benchmark process to one logical CPU on Linux.",
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

    repeats = args.repeats if args.repeats is not None else (5 if args.quick else 9)
    if repeats < 1 or args.warmups < 0:
        raise ValueError("repeats must be >=1 and warmups must be >=0")

    ntheta_values = [64, 256, 1024] if args.quick else [64, 128, 256, 512, 1024, 2048, 4096]
    nr_values = [64, 256, 1024] if args.quick else [64, 128, 256, 512, 1024, 1536, 2048]
    nxy_values = [64, 128, 256] if args.quick else [64, 96, 128, 192, 256, 384, 512]
    strategy_nxy_values = [128, 256] if args.quick else [64, 128, 256, 384, 512]

    studies = {
        "callable_angular_resolution": run_study(
            [
                {"input_type": "callable", "nxy": 129, "nr": 256, "ntheta": n, "interp_method": "cubic"}
                for n in ntheta_values
            ],
            repeats,
            args.warmups,
        ),
        "callable_radial_resolution": run_study(
            [
                {"input_type": "callable", "nxy": 129, "nr": n, "ntheta": 512, "interp_method": "cubic"}
                for n in nr_values
            ],
            repeats,
            args.warmups,
        ),
        "sampled_cartesian_resolution": run_study(
            [
                {"input_type": "sampled", "nxy": n, "nr": 256, "ntheta": 512, "interp_method": method}
                for method in ("linear", "cubic")
                for n in nxy_values
            ],
            repeats,
            args.warmups,
        ),
        "sampled_radial_resolution_strategy_cost": run_study(
            [
                {
                    "input_type": "sampled",
                    "nxy": nxy,
                    "nr": nr,
                    "ntheta": 512,
                    "interp_method": "cubic",
                    "nr_strategy": strategy,
                }
                for nxy in strategy_nxy_values
                for strategy, nr in (
                    ("half_cartesian", max(2, nxy // 2)),
                    ("match_cartesian", nxy),
                    ("double_cartesian", 2 * nxy),
                )
            ],
            repeats,
            args.warmups,
        ),
    }

    methodology = {
        "threads": "OMP/MKL/OpenBLAS/NumExpr thread counts set to 1 before NumPy/SciPy import when not already specified.",
        "summary": "Wall-clock perf_counter timings; warmups excluded; median and interquartile range reported.",
        "stage_definitions": {
            "prepare_field": "Input evaluation/sampling, Cartesian power/origin preparation, and polar-grid construction/interpolation.",
            "decompose": "Angular FFT, radial powers, adaptive selection, reconstruction, and radial cutoff diagnostics.",
            "total": "Complete PolarDecomposition construction including validation and the two timed internal stages.",
            "other_overhead": "Residual total time outside prepare_field and decompose, including input validation and timing/scheduling noise; reported explicitly rather than silently assigned to either stage.",
        },
        "repeats": repeats,
        "warmups": args.warmups,
        "cpu_affinity_request": args.cpu,
        "scaling_interpretation": "The angular FFT contributes O(Nr*Ntheta*log Ntheta), but the complete implementation also contains O(Nr*Ntheta) work. Benchmark plots therefore show empirical affine fits to the measured total runtime rather than forcing an Ntheta*log(Ntheta) reference curve.",
    }
    path = write_json(make_payload("runtime_benchmark", args.quick, studies, methodology), args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
