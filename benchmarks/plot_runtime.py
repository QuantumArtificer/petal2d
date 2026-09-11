#!/usr/bin/env python3
"""Plot PETAL2D runtime-benchmark JSON without imposing asymptotic fits."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from _common import REPO_ROOT, ensure_figures_dir


def load(path: str) -> dict:
    p = Path(path)
    if not p.is_absolute():
        p = REPO_ROOT / p
    return json.loads(p.read_text())


def med(row: dict, stage="total") -> float:
    return float(row["timing_seconds"][stage]["median"])


def iqr(row: dict, stage="total") -> float:
    return float(row["timing_seconds"][stage]["iqr"])


def affine_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    coeff = np.polyfit(x, y, 1)
    predicted = np.polyval(coeff, x)
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return coeff, r2


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="benchmarks/results/runtime_benchmark.json")
    parser.add_argument("--output-dir", default="benchmarks/results/figures")
    args = parser.parse_args()

    payload = load(args.input)
    out = ensure_figures_dir(args.output_dir)

    rows = sorted(payload["studies"]["callable_angular_resolution"], key=lambda r: r["ntheta"])
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    x = np.array([r["ntheta"] for r in rows], dtype=float)
    y = np.array([med(r) for r in rows], dtype=float)
    ax.errorbar(x, y, yerr=[0.5 * iqr(r) for r in rows], marker="o", capsize=2, label="Measured")
    fit, r2 = affine_fit(x, y)
    ax.plot(x, np.polyval(fit, x), linestyle="--", label=rf"Affine fit in $N_\theta$ ($R^2={r2:.4f}$)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"Angular resolution $N_\theta$")
    ax.set_ylabel("Median runtime (s)")
    ax.set_title("Callable-input angular scaling")
    ax.grid(True, which="both", alpha=0.20)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "runtime_callable_ntheta.png", dpi=220)
    plt.close(fig)

    rows = sorted(payload["studies"]["callable_radial_resolution"], key=lambda r: r["nr"])
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    x = np.array([r["nr"] for r in rows], dtype=float)
    y = np.array([med(r) for r in rows], dtype=float)
    ax.errorbar(x, y, yerr=[0.5 * iqr(r) for r in rows], marker="o", capsize=2, label="Measured")
    fit, r2 = affine_fit(x, y)
    ax.plot(x, np.polyval(fit, x), linestyle="--", label=rf"Affine fit in $N_r$ ($R^2={r2:.4f}$)")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xlabel(r"Radial resolution $N_r$")
    ax.set_ylabel("Median runtime (s)")
    ax.set_title("Callable-input radial scaling")
    ax.grid(True, which="both", alpha=0.20)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "runtime_callable_nr.png", dpi=220)
    plt.close(fig)

    rows_all = payload["studies"]["sampled_cartesian_resolution"]
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    for method in ("linear", "cubic"):
        rows = sorted((r for r in rows_all if r["interp_method"] == method), key=lambda r: r["nxy"])
        ax.loglog(
            [r["nxy"] for r in rows],
            [med(r) for r in rows],
            marker="o",
            label=method.capitalize(),
        )
    ax.set_xlabel(r"Cartesian resolution $N_x=N_y$")
    ax.set_ylabel("Median runtime (s)")
    ax.set_title("Sampled-input Cartesian scaling")
    ax.grid(True, which="both", alpha=0.20)
    ax.legend(frameon=False)
    fig.savefig(out / "runtime_sampled_nxy.png", dpi=220)
    plt.close(fig)

    rows = sorted(payload["studies"]["sampled_cartesian_resolution"], key=lambda r: (r["interp_method"], r["nxy"]))
    fig, ax = plt.subplots(figsize=(5.6, 3.9), layout="constrained")
    for method in ("linear", "cubic"):
        subset = [r for r in rows if r["interp_method"] == method]
        ax.plot(
            [r["nxy"] for r in subset],
            [med(r, "prepare_field") / med(r) for r in subset],
            marker="o",
            label=f"{method.capitalize()} preparation fraction",
        )
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel(r"Cartesian resolution $N_x=N_y$")
    ax.set_ylabel("Fraction of total runtime")
    ax.set_title("Cost of sampled-field preparation and interpolation")
    ax.grid(True, alpha=0.20)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "runtime_sampled_stage_fraction.png", dpi=220)
    plt.close(fig)


    strategy_rows = payload["studies"].get("sampled_radial_resolution_strategy_cost", [])
    if strategy_rows:
        fig, ax = plt.subplots(figsize=(5.8, 4.0), layout="constrained")
        for nxy in sorted({int(r["nxy"]) for r in strategy_rows}):
            subset = [r for r in strategy_rows if int(r["nxy"]) == nxy]
            by_strategy = {r["nr_strategy"]: r for r in subset}
            baseline = med(by_strategy["match_cartesian"])
            x = [0.5, 1.0, 2.0]
            y = [
                med(by_strategy[name]) / baseline
                for name in ("half_cartesian", "match_cartesian", "double_cartesian")
            ]
            ax.plot(x, y, marker="o", label=rf"$N_{{xy}}={nxy}$")
        ax.axhline(1.0, linewidth=0.9, linestyle="--", alpha=0.65)
        ax.set_xticks([0.5, 1.0, 2.0], [r"$N_r=N_{xy}/2$", r"$N_r=N_{xy}$", r"$N_r=2N_{xy}$"])
        ax.set_ylabel(r"Runtime relative to $N_r=N_{xy}$")
        ax.set_title("Runtime cost of radial-resolution strategy")
        ax.grid(True, alpha=0.20)
        ax.legend(frameon=False, fontsize=8, ncols=2)
        fig.savefig(out / "runtime_default_nr_cost.png", dpi=220)
        plt.close(fig)

    print(f"Wrote figures to {out}")


if __name__ == "__main__":
    main()
