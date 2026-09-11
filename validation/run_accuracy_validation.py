#!/usr/bin/env python3
"""Numerical-accuracy validation for PETAL2D.

This script writes JSON only. Plotting is handled by
``validation/plot_accuracy_validation.py``.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from petal2d import PolarDecomposition

from _common import make_payload, write_json


def mode_index(dec: PolarDecomposition, m: int) -> int:
    matches = np.flatnonzero(dec.m_all == int(m))
    if matches.size != 1:
        raise RuntimeError(f"Could not identify unique FFT bin for m={m}.")
    return int(matches[0])


def smooth_c3_exact_powers(beta: float, alpha: float, rmax: float) -> dict[int, float]:
    # f = exp(-beta r^2) [1 + alpha r^3 cos(3 theta)]
    # rho_0 = exp(-beta r^2)
    # rho_{+/-3} = alpha/2 * r^3 exp(-beta r^2)
    p0 = quad(lambda r: np.exp(-2.0 * beta * r * r) * r, 0.0, rmax, epsabs=1e-14)[0]
    p3 = quad(
        lambda r: (alpha * alpha / 4.0) * r**7 * np.exp(-2.0 * beta * r * r),
        0.0,
        rmax,
        epsabs=1e-14,
    )[0]
    return {0: p0, 3: p3, -3: p3}


def sampled_cartesian_resolution_study(quick: bool) -> list[dict]:
    beta = 0.60
    alpha = 0.20
    extent = 5.0
    rmax = extent
    nr = 192
    ntheta = 256
    nxy_values = [51, 101, 161] if quick else [41, 51, 61, 81, 101, 141, 201, 281]
    interpolation_methods = ["linear", "cubic"]

    exact_powers = smooth_c3_exact_powers(beta, alpha, rmax)
    exact_total = sum(exact_powers.values())
    exact_fractions = {m: p / exact_total for m, p in exact_powers.items()}

    rows: list[dict] = []
    for method in interpolation_methods:
        for nxy in nxy_values:
            axis = np.linspace(-extent, extent, nxy)
            X, Y = np.meshgrid(axis, axis, indexing="ij")
            r2 = X * X + Y * Y
            c3 = X**3 - 3.0 * X * Y**2
            sampled = np.exp(-beta * r2) * (1.0 + alpha * c3)

            dec = PolarDecomposition(
                sampled,
                axis,
                axis,
                Nr=nr,
                Ntheta=ntheta,
                rmax=rmax,
                origin=(0.0, 0.0),
                recon_err_tol=1.0,
                interp_method=method,
            )

            recovered = {
                m: float(dec.power_fracs_all[mode_index(dec, m)])
                for m in exact_fractions
            }
            expected_fraction_recovered = float(sum(recovered.values()))
            max_fraction_abs_error = float(
                max(abs(recovered[m] - exact_fractions[m]) for m in exact_fractions)
            )

            rows.append(
                {
                    "Nxy": int(nxy),
                    "cartesian_points": int(nxy * nxy),
                    "Nr": nr,
                    "Ntheta": ntheta,
                    "interp_method": method,
                    "expected_modes": [0, 3, -3],
                    "exact_power_fractions": {str(m): float(v) for m, v in exact_fractions.items()},
                    "recovered_power_fractions": {str(m): float(v) for m, v in recovered.items()},
                    "max_expected_mode_power_fraction_abs_error": max_fraction_abs_error,
                    "interpolation_leakage_fraction": max(0.0, 1.0 - expected_fraction_recovered),
                    "expected_mode_reconstruction_error_percent": float(
                        dec.reconstruction_error([0, 3, -3])
                    ),
                    "domain_consistency": float(dec.domain_consistency),
                    "parseval_relative_error": float(dec.parseval_relative_error),
                }
            )
    return rows


def default_radial_resolution_strategy_study(quick: bool) -> list[dict]:
    """End-to-end sampled-field accuracy for three Nr-vs-Nxy strategies.

    This study documents the accuracy tradeoff behind PETAL2D's default
    ``Nr = min(Nx, Ny)`` rule. The same smooth C3 field and exact harmonic
    powers used by the Cartesian-resolution study are used here, while only
    the radial-resolution strategy is changed. The half-resolution rule is
    retained as a lower-cost comparison and ``2*Nxy`` as a higher-resolution
    reference.
    """
    beta = 0.60
    alpha = 0.20
    extent = 5.0
    rmax = extent
    ntheta = 256
    nxy_values = [51, 101, 161] if quick else [41, 61, 81, 101, 141, 201, 281]
    strategies = (
        ("half_cartesian", lambda n: max(2, n // 2)),
        ("match_cartesian", lambda n: n),
        ("double_cartesian", lambda n: 2 * n),
    )

    exact_powers = smooth_c3_exact_powers(beta, alpha, rmax)
    exact_total = sum(exact_powers.values())
    exact_fractions = {m: p / exact_total for m, p in exact_powers.items()}

    rows: list[dict] = []
    for nxy in nxy_values:
        axis = np.linspace(-extent, extent, nxy)
        X, Y = np.meshgrid(axis, axis, indexing="ij")
        r2 = X * X + Y * Y
        c3 = X**3 - 3.0 * X * Y**2
        sampled = np.exp(-beta * r2) * (1.0 + alpha * c3)

        for strategy_name, nr_rule in strategies:
            nr = int(nr_rule(nxy))
            dec = PolarDecomposition(
                sampled,
                axis,
                axis,
                Nr=nr,
                Ntheta=ntheta,
                rmax=rmax,
                origin=(0.0, 0.0),
                recon_err_tol=1.0,
                interp_method="cubic",
            )
            recovered = {
                m: float(dec.power_fracs_all[mode_index(dec, m)])
                for m in exact_fractions
            }
            expected_fraction_recovered = float(sum(recovered.values()))
            rows.append(
                {
                    "Nxy": int(nxy),
                    "Nr": nr,
                    "Ntheta": ntheta,
                    "nr_strategy": strategy_name,
                    "nr_over_nxy": float(nr / nxy),
                    "max_expected_mode_power_fraction_abs_error": float(
                        max(abs(recovered[m] - exact_fractions[m]) for m in exact_fractions)
                    ),
                    "interpolation_leakage_fraction": max(
                        0.0, 1.0 - expected_fraction_recovered
                    ),
                    "expected_mode_reconstruction_error_percent": float(
                        dec.reconstruction_error([0, 3, -3])
                    ),
                    "domain_consistency": float(dec.domain_consistency),
                    "domain_consistency_abs_error_from_unity": float(
                        abs(dec.domain_consistency - 1.0)
                    ),
                }
            )
    return rows


def gaussian_exact_support(beta: float, rmax: float, power_tail: float, amplitude_threshold: float) -> tuple[float, float, float]:
    finite_norm_factor = 1.0 - np.exp(-2.0 * beta * rmax * rmax)
    target = 1.0 - power_tail
    rhs = max(np.finfo(float).tiny, 1.0 - target * finite_norm_factor)
    r_power = float(np.sqrt(-np.log(rhs) / (2.0 * beta)))
    r_amplitude = float(min(rmax, np.sqrt(-np.log(amplitude_threshold) / beta)))
    return r_power, r_amplitude, max(r_power, r_amplitude)


def radial_resolution_study(quick: bool) -> list[dict]:
    beta = 0.70
    extent = 6.0
    rmax = extent
    ntheta = 192
    nr_values = [48, 96, 192] if quick else [24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
    axis = np.linspace(-extent, extent, 401)

    exact_power = quad(lambda r: np.exp(-2.0 * beta * r * r) * r, 0.0, rmax, epsabs=1e-14)[0]

    rows: list[dict] = []
    for nr in nr_values:
        dec = PolarDecomposition(
            lambda x, y: np.exp(-beta * (x * x + y * y)),
            axis,
            axis,
            Nr=nr,
            Ntheta=ntheta,
            rmax=rmax,
            origin=(0.0, 0.0),
            recon_err_tol=1e-10,
        )
        rp_exact, ra_exact, rc_exact = gaussian_exact_support(
            beta,
            rmax,
            dec.radial_power_tail_fraction,
            dec.radial_relative_amplitude_threshold,
        )
        rows.append(
            {
                "Nr": int(nr),
                "Ntheta": ntheta,
                "radial_step": float(dec.r[1] - dec.r[0]),
                "mode_power": float(dec.powers[0]),
                "exact_mode_power": float(exact_power),
                "relative_mode_power_error": float(abs(dec.powers[0] - exact_power) / exact_power),
                "domain_consistency": float(dec.domain_consistency),
                "parseval_relative_error": float(dec.parseval_relative_error),
                "radial_power_support_radius": float(dec.radial_power_support_radius[0]),
                "exact_radial_power_support_radius": rp_exact,
                "radial_power_support_radius_abs_error": float(abs(dec.radial_power_support_radius[0] - rp_exact)),
                "radial_amplitude_support_radius": float(dec.radial_amplitude_support_radius[0]),
                "exact_radial_amplitude_support_radius": ra_exact,
                "radial_amplitude_support_radius_abs_error": float(abs(dec.radial_amplitude_support_radius[0] - ra_exact)),
                "cutoff_radius": float(dec.cutoff_radius[0]),
                "exact_cutoff_radius": rc_exact,
                "cutoff_radius_abs_error": float(abs(dec.cutoff_radius[0] - rc_exact)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true", help="Run a reduced validation grid.")
    parser.add_argument(
        "--output",
        default="validation/results/accuracy_validation.json",
        help="JSON output path, relative to the repository root unless absolute.",
    )
    args = parser.parse_args()

    studies = {
        "sampled_cartesian_resolution": sampled_cartesian_resolution_study(args.quick),
        "default_radial_resolution_strategy": default_radial_resolution_strategy_study(args.quick),
        "radial_resolution": radial_resolution_study(args.quick),
    }
    path = write_json(make_payload("accuracy_validation", args.quick, studies), args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
