#!/usr/bin/env python3
"""Validate angular sampling, Nyquist behavior, and Fourier aliasing."""

from __future__ import annotations

import argparse

import numpy as np

from petal2d import PolarDecomposition

from _common import make_payload, write_json


def localized_pure_harmonic(m: int, r0: float = 1.3):
    m = int(m)

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        if m == 0:
            amplitude = np.exp(-0.5 * (r / r0) ** 2)
        else:
            s = r / r0
            amplitude = np.zeros_like(r, dtype=float)
            mask = s > 0.0
            log_amp = m * np.log(s[mask]) - 0.5 * m * (s[mask] ** 2 - 1.0)
            amplitude[mask] = np.exp(log_amp)
        return amplitude * np.exp(1j * m * theta)

    return field


def expected_fft_alias(m: int, ntheta: int) -> int:
    bins = np.rint(np.fft.fftfreq(ntheta, d=1.0 / ntheta)).astype(int)
    return int(bins[m % ntheta])


def complex_aliasing_study(quick: bool) -> list[dict]:
    ntheta_values = [16, 32, 64, 128, 256] if quick else [16, 24, 32, 48, 64, 96, 128, 192, 256]
    m_values = [0, 1, 4, 7, 8, 16, 31, 32, 63, 64, 100] if quick else list(range(0, 101))
    axis = np.linspace(-4.0, 4.0, 81)

    rows: list[dict] = []
    for ntheta in ntheta_values:
        for m in m_values:
            dec = PolarDecomposition(
                localized_pure_harmonic(m),
                axis,
                axis,
                Nr=64,
                Ntheta=ntheta,
                rmax=4.0,
                origin=(0.0, 0.0),
                recon_err_tol=1e-8,
            )
            dominant_index = int(np.argmax(dec.powers_all))
            dominant_m = int(dec.m_all[dominant_index])
            expected_alias = expected_fft_alias(m, ntheta)
            dominant_fraction = float(dec.power_fracs_all[dominant_index])
            rows.append(
                {
                    "true_m": int(m),
                    "Ntheta": int(ntheta),
                    "m_nyquist": int(dec.m_nyquist),
                    "expected_discrete_alias_m": expected_alias,
                    "recovered_dominant_m": dominant_m,
                    "alias_prediction_correct": bool(dominant_m == expected_alias),
                    "unaliased": bool(expected_alias == m),
                    "dominant_power_fraction": dominant_fraction,
                    "spectral_leakage_fraction": max(0.0, 1.0 - dominant_fraction),
                    "parseval_relative_error": float(dec.parseval_relative_error),
                }
            )
    return rows


def real_pair_boundary_study(quick: bool) -> list[dict]:
    ntheta_values = [16, 32, 64] if quick else [16, 24, 32, 48, 64, 96, 128]
    axis = np.linspace(-4.0, 4.0, 81)
    rows: list[dict] = []

    for ntheta in ntheta_values:
        candidates = sorted({1, 3, max(1, ntheta // 2 - 1), ntheta // 2})
        for m in candidates:
            def field(x, y, m=m):
                r = np.hypot(x, y)
                theta = np.arctan2(y, x)
                s = r / 1.3
                amplitude = np.zeros_like(r)
                mask = s > 0
                log_amp = m * np.log(s[mask]) - 0.5 * m * (s[mask] ** 2 - 1.0)
                amplitude[mask] = np.exp(log_amp)
                return amplitude * np.cos(m * theta)

            dec = PolarDecomposition(
                field,
                axis,
                axis,
                Nr=64,
                Ntheta=ntheta,
                rmax=4.0,
                origin=(0.0, 0.0),
                recon_err_tol=1e-8,
            )
            rows.append(
                {
                    "m": int(m),
                    "Ntheta": int(ntheta),
                    "m_nyquist": int(dec.m_nyquist),
                    "nyquist_mode": None if dec.nyquist_mode is None else int(dec.nyquist_mode),
                    "selected_pairs": [list(map(int, p)) for p in dec.selected_pairs],
                    "retained_modes": [int(v) for v in dec.m_sorted],
                    "is_singleton_nyquist_channel": bool(
                        ntheta % 2 == 0 and m == ntheta // 2 and dec.selected_pairs == [(-ntheta // 2,)]
                    ),
                    "reconstruction_error_percent": float(dec.recon_error_measured),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--output", default="validation/results/angular_resolution_validation.json")
    args = parser.parse_args()

    studies = {
        "complex_aliasing": complex_aliasing_study(args.quick),
        "real_pair_boundary": real_pair_boundary_study(args.quick),
    }
    path = write_json(make_payload("angular_resolution_validation", args.quick, studies), args.output)
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
