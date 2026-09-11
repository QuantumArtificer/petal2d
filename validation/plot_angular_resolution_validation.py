#!/usr/bin/env python3
"""Plot angular-resolution validation JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt

from _common import REPO_ROOT, ensure_figures_dir


def load(path: str) -> dict:
    p = Path(path)
    if not p.is_absolute():
        p = REPO_ROOT / p
    return json.loads(p.read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="validation/results/angular_resolution_validation.json")
    parser.add_argument("--output-dir", default="validation/results/figures")
    args = parser.parse_args()

    payload = load(args.input)
    rows = payload["studies"]["complex_aliasing"]
    out = ensure_figures_dir(args.output_dir)

    ntheta_values = sorted({r["Ntheta"] for r in rows})
    fig, ax = plt.subplots(figsize=(6.2, 4.4), layout="constrained")
    for ntheta in ntheta_values:
        subset = sorted((r for r in rows if r["Ntheta"] == ntheta), key=lambda r: r["true_m"])
        ax.plot(
            [r["true_m"] for r in subset],
            [r["recovered_dominant_m"] for r in subset],
            marker=".",
            linewidth=1.0,
            label=rf"$N_\theta={ntheta}$",
        )
    max_m = max(r["true_m"] for r in rows)
    ax.plot([0, max_m], [0, max_m], linestyle="--", linewidth=1.0, color="0.35", label="No aliasing")
    ax.set_xlabel(r"True angular harmonic $m$")
    ax.set_ylabel(r"Recovered discrete harmonic $m$")
    ax.set_title("Angular aliasing follows the discrete FFT spectrum")
    ax.grid(True, alpha=0.20)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.savefig(out / "angular_aliasing_map.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.8, 4.0), layout="constrained")
    for ntheta in ntheta_values:
        subset = sorted((r for r in rows if r["Ntheta"] == ntheta and r["unaliased"]), key=lambda r: r["true_m"])
        if not subset:
            continue
        ax.semilogy(
            [r["true_m"] for r in subset],
            [max(r["spectral_leakage_fraction"], 1e-18) for r in subset],
            marker="o",
            linewidth=1.0,
            label=rf"$N_\theta={ntheta}$",
        )
    ax.set_xlabel(r"Resolved angular harmonic $m$")
    ax.set_ylabel("Spectral leakage fraction")
    ax.set_title("Pure-harmonic recovery below the Nyquist limit")
    ax.grid(True, which="both", alpha=0.20)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(out / "angular_unaliased_leakage.png", dpi=220)
    plt.close(fig)

    print(f"Wrote figures to {out}")


if __name__ == "__main__":
    main()
