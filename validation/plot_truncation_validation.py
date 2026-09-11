#!/usr/bin/env python3
"""Plot adaptive-truncation validation JSON."""

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
    parser.add_argument("--input", default="validation/results/truncation_validation.json")
    parser.add_argument("--output-dir", default="validation/results/figures")
    args = parser.parse_args()

    payload = load(args.input)
    out = ensure_figures_dir(args.output_dir)

    fig, ax = plt.subplots(figsize=(5.8, 4.1), layout="constrained")
    for key, label, marker in [
        ("complex_independent_modes", "Complex modes", "o"),
        ("real_conjugate_pairs", "Real conjugate pairs", "s"),
    ]:
        rows = sorted(payload["studies"][key], key=lambda r: r["recon_err_tol_percent"], reverse=True)
        ax.plot(
            [r["expected_error_percent"] for r in rows],
            [r["measured_error_percent"] for r in rows],
            marker=marker,
            linestyle="none",
            label=label,
        )
    limits = ax.get_xlim()
    lo = max(0.0, limits[0])
    hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lo, hi], [lo, hi], linestyle="--", color="0.35", linewidth=1.0, label="Exact agreement")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Theoretical retained-spectrum error (%)")
    ax.set_ylabel("Measured weighted $L^2$ error (%)")
    ax.set_title("Adaptive truncation error certification")
    ax.grid(True, alpha=0.20)
    ax.legend(frameon=False)
    fig.savefig(out / "truncation_error_certification.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.8, 4.1), layout="constrained")
    for key, label, marker in [
        ("complex_independent_modes", "Complex modes", "o"),
        ("real_conjugate_pairs", "Real conjugate pairs", "s"),
    ]:
        rows = sorted(payload["studies"][key], key=lambda r: r["recon_err_tol_percent"], reverse=True)
        ax.step(
            [r["recon_err_tol_percent"] for r in rows],
            [r["retained_mode_count"] for r in rows],
            where="mid",
            marker=marker,
            label=label,
        )
    ax.invert_xaxis()
    ax.set_xlabel("Requested reconstruction-error tolerance (%)")
    ax.set_ylabel("Retained angular modes")
    ax.set_title("Adaptive spectral compression")
    ax.grid(True, alpha=0.20)
    ax.legend(frameon=False)
    fig.savefig(out / "truncation_retained_modes.png", dpi=220)
    plt.close(fig)

    print(f"Wrote figures to {out}")


if __name__ == "__main__":
    main()
