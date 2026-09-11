#!/usr/bin/env python3
"""Plot PETAL2D incremental peak-memory benchmark JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from _common import REPO_ROOT, ensure_figures_dir


MIB = 1024.0**2


def load(path: str) -> dict:
    p = Path(path)
    if not p.is_absolute():
        p = REPO_ROOT / p
    return json.loads(p.read_text())


def median_mib(row: dict) -> float:
    return float(row["petal2d_peak_incremental_rss_bytes"]["median"]) / MIB


def affine_fit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float]:
    coeff = np.polyfit(x, y, 1)
    predicted = np.polyval(coeff, x)
    ss_res = float(np.sum((y - predicted) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return coeff, r2


def q25_mib(row: dict) -> float:
    return float(row["petal2d_peak_incremental_rss_bytes"]["q25"]) / MIB


def q75_mib(row: dict) -> float:
    return float(row["petal2d_peak_incremental_rss_bytes"]["q75"]) / MIB


def plot_study(rows, xkey, xlabel, title, output, *, fit_affine: bool):
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    x = np.array([r[xkey] for r in rows], dtype=float)
    y = np.array([median_mib(r) for r in rows], dtype=float)
    lo = y - np.array([q25_mib(r) for r in rows], dtype=float)
    hi = np.array([q75_mib(r) for r in rows], dtype=float) - y
    ax.errorbar(x, y, yerr=np.vstack([lo, hi]), marker="o", capsize=2.5, label="Median ± IQR")
    if fit_affine and len(rows) >= 2:
        fit, r2 = affine_fit(x, y)
        ax.plot(
            x,
            np.polyval(fit, x),
            linestyle="--",
            label=rf"Affine fit ($R^2={r2:.4f}$)",
        )
        bytes_per_point = fit[0] * MIB
        ax.text(
            0.03,
            0.96,
            f"Slope = {bytes_per_point:.1f} bytes/point",
            transform=ax.transAxes,
            va="top",
            fontsize=8,
        )
    ax.set_xlabel(xlabel)
    ax.set_ylabel("PETAL2D peak RSS above input baseline (MiB)")
    ax.set_title(title)
    ax.grid(True, alpha=0.20)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="benchmarks/results/memory_benchmark.json")
    parser.add_argument("--output-dir", default="benchmarks/results/figures")
    args = parser.parse_args()

    payload = load(args.input)
    out = ensure_figures_dir(args.output_dir)

    rows = sorted(payload["studies"]["callable_polar_grid"], key=lambda r: r["polar_points"])
    plot_study(
        rows,
        "polar_points",
        r"Polar grid points $N_rN_\theta$",
        "Incremental memory vs polar-grid size",
        out / "memory_polar_grid.png",
        fit_affine=True,
    )

    rows = sorted(
        payload["studies"]["sampled_cartesian_grid"],
        key=lambda r: r["cartesian_points"],
    )
    plot_study(
        rows,
        "cartesian_points",
        r"Cartesian input points $N_xN_y$",
        "Incremental memory vs sampled-input size",
        out / "memory_cartesian_grid.png",
        fit_affine=False,
    )


    strategy_rows = payload["studies"].get("sampled_radial_resolution_strategy_cost", [])
    if strategy_rows:
        fig, ax = plt.subplots(figsize=(5.8, 4.0), layout="constrained")
        for nxy in sorted({int(r["nxy"]) for r in strategy_rows}):
            subset = [r for r in strategy_rows if int(r["nxy"]) == nxy]
            by_strategy = {r["nr_strategy"]: r for r in subset}
            baseline = median_mib(by_strategy["match_cartesian"])
            x = [0.5, 1.0, 2.0]
            y = [
                median_mib(by_strategy[name]) / baseline
                for name in ("half_cartesian", "match_cartesian", "double_cartesian")
            ]
            ax.plot(x, y, marker="o", label=rf"$N_{{xy}}={nxy}$")
        ax.axhline(1.0, linewidth=0.9, linestyle="--", alpha=0.65)
        ax.set_xticks([0.5, 1.0, 2.0], [r"$N_r=N_{xy}/2$", r"$N_r=N_{xy}$", r"$N_r=2N_{xy}$"])
        ax.set_ylabel(r"Peak PETAL2D RSS relative to $N_r=N_{xy}$")
        ax.set_title("Memory cost of radial-resolution strategy")
        ax.grid(True, alpha=0.20)
        ax.legend(frameon=False, fontsize=8)
        fig.savefig(out / "memory_default_nr_cost.png", dpi=220)
        plt.close(fig)

    print(f"Wrote figures to {out}")


if __name__ == "__main__":
    main()
