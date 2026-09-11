#!/usr/bin/env python3
"""Plot PETAL2D accuracy-validation JSON without recomputing it."""

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


def positive(values):
    tiny = np.finfo(float).tiny
    return [max(float(v), tiny) for v in values]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="validation/results/accuracy_validation.json")
    parser.add_argument("--output-dir", default="validation/results/figures")
    args = parser.parse_args()

    payload = load(args.input)
    out = ensure_figures_dir(args.output_dir)

    sampled = payload["studies"]["sampled_cartesian_resolution"]
    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    for method in sorted({row["interp_method"] for row in sampled}):
        rows = sorted((r for r in sampled if r["interp_method"] == method), key=lambda r: r["Nxy"])
        ax.loglog(
            [r["Nxy"] for r in rows],
            positive([r["interpolation_leakage_fraction"] for r in rows]),
            marker="o",
            label=method.capitalize(),
        )
    ax.set_xlabel(r"Cartesian resolution $N_x=N_y$")
    ax.set_ylabel("Spurious harmonic power fraction")
    ax.set_title("Interpolation-induced spectral leakage")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False)
    fig.savefig(out / "accuracy_cartesian_leakage.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 3.8), layout="constrained")
    for method in sorted({row["interp_method"] for row in sampled}):
        rows = sorted((r for r in sampled if r["interp_method"] == method), key=lambda r: r["Nxy"])
        ax.loglog(
            [r["Nxy"] for r in rows],
            positive([r["max_expected_mode_power_fraction_abs_error"] for r in rows]),
            marker="o",
            label=method.capitalize(),
        )
    ax.set_xlabel(r"Cartesian resolution $N_x=N_y$")
    ax.set_ylabel("Max. expected-mode power-fraction error")
    ax.set_title("Recovery of known harmonic powers")
    ax.grid(True, which="both", alpha=0.22)
    ax.legend(frameon=False)
    fig.savefig(out / "accuracy_mode_power_error.png", dpi=220)
    plt.close(fig)

    # End-to-end study used to decide the default Nr rule.
    strategies = payload["studies"].get("default_radial_resolution_strategy", [])
    if strategies:
        labels = {
            "half_cartesian": r"$N_r=\lfloor N_{xy}/2\rfloor$",
            "match_cartesian": r"$N_r=N_{xy}$ (default)",
            "double_cartesian": r"$N_r=2N_{xy}$",
        }
        fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.7), layout="constrained")
        for strategy in ("half_cartesian", "match_cartesian", "double_cartesian"):
            rows = sorted(
                (r for r in strategies if r["nr_strategy"] == strategy),
                key=lambda r: r["Nxy"],
            )
            x = [r["Nxy"] for r in rows]
            axs[0].loglog(
                x,
                positive([r["max_expected_mode_power_fraction_abs_error"] for r in rows]),
                marker="o",
                label=labels[strategy],
            )
            axs[1].loglog(
                x,
                positive([r["expected_mode_reconstruction_error_percent"] for r in rows]),
                marker="o",
                label=labels[strategy],
            )
        axs[0].set_xlabel(r"Cartesian resolution $N_x=N_y$")
        axs[0].set_ylabel("Max. harmonic power-fraction error")
        axs[0].set_title("Known-spectrum recovery")
        axs[1].set_xlabel(r"Cartesian resolution $N_x=N_y$")
        axs[1].set_ylabel(r"Expected-mode reconstruction error (\%)")
        axs[1].set_title("End-to-end reconstruction")
        for ax in axs:
            ax.grid(True, which="both", alpha=0.22)
        axs[0].legend(frameon=False, fontsize=8)
        fig.suptitle("Effect of the radial-resolution strategy", fontsize=11)
        fig.savefig(out / "accuracy_default_nr_strategy.png", dpi=220)
        plt.close(fig)

    # Keep dimensionless and dimensional convergence errors on separate axes.
    radial = sorted(payload["studies"]["radial_resolution"], key=lambda r: r["Nr"])
    fig, axs = plt.subplots(1, 2, figsize=(9.2, 3.7), layout="constrained")
    axs[0].loglog(
        [r["Nr"] for r in radial],
        positive([r["relative_mode_power_error"] for r in radial]),
        marker="o",
    )
    axs[0].set_xlabel(r"Radial resolution $N_r$")
    axs[0].set_ylabel("Relative mode-power error")
    axs[0].set_title("Radial quadrature")

    axs[1].loglog(
        [r["Nr"] for r in radial],
        positive([r["cutoff_radius_abs_error"] for r in radial]),
        marker="s",
    )
    axs[1].set_xlabel(r"Radial resolution $N_r$")
    axs[1].set_ylabel("Absolute cutoff-radius error")
    axs[1].set_title("Cutoff-radius convergence")
    for ax in axs:
        ax.grid(True, which="both", alpha=0.22)
    fig.suptitle("Radial-resolution convergence", fontsize=11)
    fig.savefig(out / "accuracy_radial_convergence.png", dpi=220)
    plt.close(fig)

    print(f"Wrote figures to {out}")


if __name__ == "__main__":
    main()
