#!/usr/bin/env python3
"""Validate PETAL2D adaptive truncation against known harmonic powers."""

from __future__ import annotations

import argparse
from functools import lru_cache
from math import sqrt

import numpy as np
from scipy.integrate import quad

from petal2d import PolarDecomposition

from _common import make_payload, write_json


BETA = 0.75
RMAX = 5.0


@lru_cache(maxsize=None)
def radial_norm(m_abs: int) -> float:
    return quad(
        lambda r: r ** (2 * m_abs + 1) * np.exp(-2.0 * BETA * r * r),
        0.0,
        RMAX,
        epsabs=1e-14,
    )[0]


def radial_basis(m_abs: int, r: np.ndarray) -> np.ndarray:
    return r**m_abs * np.exp(-BETA * r * r) / sqrt(radial_norm(m_abs))


def complex_field_factory(power_fractions: dict[int, float]):
    phases = {m: np.exp(1j * (0.19 * m + 0.07)) for m in power_fractions}

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        out = np.zeros_like(r, dtype=complex)
        for m, fraction in power_fractions.items():
            out += sqrt(fraction) * phases[m] * radial_basis(abs(m), r) * np.exp(1j * m * theta)
        return out

    return field


def real_field_factory(pair_power_fractions: dict[int, float]):
    # Keys are nonnegative |m|. For m>0, the listed weight is the combined
    # power of the conjugate (+m,-m) pair.
    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        out = np.zeros_like(r, dtype=float)
        for m_abs, fraction in pair_power_fractions.items():
            if m_abs == 0:
                out += sqrt(fraction) * radial_basis(0, r)
            else:
                out += sqrt(2.0 * fraction) * radial_basis(m_abs, r) * np.cos(m_abs * theta)
        return out

    return field


def expected_selection(weights: list[tuple[tuple[int, ...], float]], tolerance_percent: float) -> tuple[list[tuple[int, ...]], float]:
    cumulative = 0.0
    selected: list[tuple[int, ...]] = []
    for modes, weight in weights:
        selected.append(modes)
        cumulative += weight
        error = 100.0 * sqrt(max(0.0, 1.0 - cumulative))
        if error <= tolerance_percent:
            return selected, error
    return selected, 100.0 * sqrt(max(0.0, 1.0 - cumulative))


def complex_study(quick: bool) -> list[dict]:
    powers = {0: 0.55, 1: 0.20, -2: 0.12, 4: 0.08, -6: 0.05}
    ranked = sorted([((m,), p) for m, p in powers.items()], key=lambda item: -item[1])
    tolerances = [55, 45, 30, 10, 1] if quick else [65, 55, 45, 40, 35, 30, 25, 15, 5, 1, 0.1]
    axis = np.linspace(-RMAX, RMAX, 161)
    rows = []
    for tol in tolerances:
        expected_pairs, expected_error = expected_selection(ranked, tol)
        dec = PolarDecomposition(
            complex_field_factory(powers),
            axis,
            axis,
            Nr=192,
            Ntheta=256,
            rmax=RMAX,
            origin=(0.0, 0.0),
            recon_err_tol=float(tol),
        )
        rows.append(
            {
                "recon_err_tol_percent": float(tol),
                "known_mode_power_fractions": {str(m): float(p) for m, p in powers.items()},
                "expected_selected_pairs": [list(p) for p in expected_pairs],
                "actual_selected_pairs": [list(map(int, p)) for p in dec.selected_pairs],
                "expected_error_percent": float(expected_error),
                "predicted_error_percent": float(dec.recon_error),
                "measured_error_percent": float(dec.recon_error_measured),
                "predicted_minus_expected_percent": float(dec.recon_error - expected_error),
                "measured_minus_predicted_percent": float(dec.recon_error_measured - dec.recon_error),
                "retained_power_fraction": float(dec.retained_power_fraction),
                "retained_mode_count": int(dec.Nm),
                "selection_matches_theory": bool(dec.selected_pairs == expected_pairs),
            }
        )
    return rows


def real_pair_study(quick: bool) -> list[dict]:
    pair_powers = {0: 0.50, 1: 0.25, 3: 0.15, 5: 0.10}
    ranked = sorted(
        [((0,), pair_powers[0])] + [((m, -m), p) for m, p in pair_powers.items() if m > 0],
        key=lambda item: -item[1],
    )
    tolerances = [55, 45, 30, 10, 1] if quick else [65, 55, 45, 40, 35, 30, 25, 15, 5, 1, 0.1]
    axis = np.linspace(-RMAX, RMAX, 161)
    rows = []
    for tol in tolerances:
        expected_pairs, expected_error = expected_selection(ranked, tol)
        dec = PolarDecomposition(
            real_field_factory(pair_powers),
            axis,
            axis,
            Nr=192,
            Ntheta=256,
            rmax=RMAX,
            origin=(0.0, 0.0),
            recon_err_tol=float(tol),
        )
        rows.append(
            {
                "recon_err_tol_percent": float(tol),
                "known_pair_power_fractions": {str(m): float(p) for m, p in pair_powers.items()},
                "expected_selected_pairs": [list(p) for p in expected_pairs],
                "actual_selected_pairs": [list(map(int, p)) for p in dec.selected_pairs],
                "expected_error_percent": float(expected_error),
                "predicted_error_percent": float(dec.recon_error),
                "measured_error_percent": float(dec.recon_error_measured),
                "predicted_minus_expected_percent": float(dec.recon_error - expected_error),
                "measured_minus_predicted_percent": float(dec.recon_error_measured - dec.recon_error),
                "retained_power_fraction": float(dec.retained_power_fraction),
                "retained_mode_count": int(dec.Nm),
                "selection_matches_theory": bool(dec.selected_pairs == expected_pairs),
                "reconstruction_is_real": bool(np.isrealobj(dec.f_recon)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", default="validation/results/truncation_validation.json")
    args = parser.parse_args()

    studies = {
        "complex_independent_modes": complex_study(args.quick),
        "real_conjugate_pairs": real_pair_study(args.quick),
    }
    path = write_json(make_payload("truncation_validation", args.quick, studies), args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
