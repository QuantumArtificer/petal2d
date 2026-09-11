import numpy as np
import pytest
from scipy.integrate import trapezoid

from petal2d import PolarDecomposition


def _grid(n=181, extent=5.0):
    axis = np.linspace(-extent, extent, n)
    return axis, axis.copy()


def test_rotation_covariance_of_harmonic_coefficients():
    x, y = _grid()
    m = 4
    phi = 0.37

    def reference(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m * np.exp(-r**2) * np.cos(m * theta)

    def rotated(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m * np.exp(-r**2) * np.cos(m * (theta - phi))

    a = PolarDecomposition(reference, x, y, Nr=120, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=1e-8)
    b = PolarDecomposition(rotated, x, y, Nr=120, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=1e-8)

    mask = np.abs(a.rho[m]) > 1e-8 * np.max(np.abs(a.rho[m]))
    ratio_plus = b.rho[m][mask] / a.rho[m][mask]
    ratio_minus = b.rho[-m][mask] / a.rho[-m][mask]
    assert np.allclose(ratio_plus, np.exp(-1j * m * phi), atol=1e-11, rtol=1e-11)
    assert np.allclose(ratio_minus, np.exp(1j * m * phi), atol=1e-11, rtol=1e-11)


def test_translation_is_removed_by_centroid_origin():
    x = np.linspace(-6.0, 6.0, 201)
    y = np.linspace(-6.0, 6.0, 201)

    def shifted(x, y):
        return np.exp(-((x - 1.2) ** 2 + (y + 0.8) ** 2))

    dec = PolarDecomposition(shifted, x, y, Nr=100, Ntheta=192, recon_err_tol=1e-7)
    assert dec.origin == pytest.approx((1.2, -0.8), abs=5e-5)
    assert dec.m_sorted == [0]


def test_conjugacy_of_complete_spectrum_for_real_field():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.2 * np.cos(3 * theta) - 0.1 * np.sin(5 * theta))

    dec = PolarDecomposition(field, x, y, Nr=100, Ntheta=255, origin=(0.0, 0.0), recon_err_tol=1e-9)
    lookup = {int(m): i for i, m in enumerate(dec.m_all)}
    for m in range(1, dec.m_nyquist + 1):
        assert np.allclose(
            dec.rho_all[:, lookup[-m]],
            np.conjugate(dec.rho_all[:, lookup[m]]),
            atol=1e-12,
            rtol=1e-12,
        )


def test_tighter_tolerance_never_reduces_retained_information():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        g = np.exp(-r**2)
        return g * (
            1.0
            + 0.5 * np.exp(1j * theta)
            + 0.25 * np.exp(-2j * theta)
            + 0.1 * np.exp(4j * theta)
        )

    tolerances = [50.0, 30.0, 15.0, 5.0]
    results = [
        PolarDecomposition(field, x, y, Nr=100, Ntheta=192, origin=(0.0, 0.0), recon_err_tol=tol)
        for tol in tolerances
    ]
    counts = [d.Nm for d in results]
    errors = [d.recon_error for d in results]
    fractions = [d.retained_power_fraction for d in results]
    assert counts == sorted(counts)
    assert errors == sorted(errors, reverse=True)
    assert fractions == sorted(fractions)


def test_parseval_full_spectrum_matches_direct_polar_norm():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-0.8 * r**2) * (1.0 + 0.3j * np.exp(2j * theta))

    dec = PolarDecomposition(field, x, y, Nr=111, Ntheta=193, origin=(0.0, 0.0))
    assert dec.parseval_relative_error < 1e-13
    assert dec.polar_power == pytest.approx(2.0 * np.pi * dec.total_power, rel=1e-13)


def test_default_real_reconstruction_is_real_but_asymmetric_mode_is_complex():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.3 * np.cos(3.0 * theta))

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=100,
        Ntheta=192,
        origin=(0.0, 0.0),
        recon_err_tol=1e-9,
    )
    assert np.isrealobj(dec.f_recon)
    assert np.isrealobj(dec.reconstruct())
    assert np.isrealobj(dec.reconstruct("all"))
    assert np.iscomplexobj(dec.reconstruct([3]))


def test_reconstruction_error_public_api_matches_stored_diagnostic():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (
            1.0 + 0.35 * np.exp(2j * theta) + 0.12 * np.exp(-5j * theta)
        )

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=100,
        Ntheta=192,
        origin=(0.0, 0.0),
        recon_err_tol=15.0,
    )
    assert dec.reconstruction_error() == pytest.approx(
        dec.recon_error_measured, rel=1e-13, abs=1e-13
    )
    assert dec.reconstruction_error("all") < 1e-11
    assert dec.reconstruction_error([0]) > dec.reconstruction_error()


def test_selected_pair_powers_are_consistent_with_selected_modes():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.4 * np.cos(3.0 * theta))

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=100,
        Ntheta=192,
        origin=(0.0, 0.0),
        recon_err_tol=20.0,
    )
    assert dec.selected_pairs == [(0,), (3, -3)]
    assert len(dec.selected_pair_powers) == len(dec.selected_pairs)
    assert np.sum(dec.selected_pair_powers) == pytest.approx(np.sum(dec.powers))
    assert np.sum(dec.selected_pair_power_fracs) == pytest.approx(
        dec.retained_power_fraction
    )


def test_normalization_connects_domain_consistency_and_polar_power():
    x, y = _grid(n=241, extent=6.0)

    def field(x, y):
        return np.exp(-0.5 * (x**2 + y**2))

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=180,
        Ntheta=256,
        origin=(0.0, 0.0),
        normalize=True,
        rmax=1.25,
    )
    assert dec.cartesian_power == pytest.approx(1.0, rel=1e-12)
    assert dec.polar_power == pytest.approx(dec.domain_consistency, rel=1e-13)


def test_px_orbital_selects_only_physical_conjugate_pair():
    x = np.linspace(-6.0, 6.0, 121)
    y = np.linspace(-6.0, 6.0, 121)

    def px(x, y):
        return x * np.exp(-0.5 * (x**2 + y**2))

    dec = PolarDecomposition(px, x, y, Nr=120, Ntheta=192, recon_err_tol=1e-8)
    assert dec.selected_pairs == [(1, -1)]
    assert set(dec.m_sorted) == {1, -1}
    assert dec.reconstruction_error() < 1e-10


def test_spectrum_summary_defaults_to_retained_modes(capsys):
    x = np.linspace(-4.0, 4.0, 81)
    y = np.linspace(-4.0, 4.0, 81)

    def field(x, y):
        return np.exp(-(x**2 + y**2))

    dec = PolarDecomposition(field, x, y, Nr=80, Ntheta=128, recon_err_tol=1e-8)
    rows = dec.spectrum_table()
    assert [row["m"] for row in rows] == [0]
    dec.print_spectrum()
    output = capsys.readouterr().out
    assert "     0" in output
    assert "    64" not in output


def test_default_radial_support_uses_both_power_and_amplitude_criteria():
    x = np.linspace(-6.0, 6.0, 161)
    y = np.linspace(-6.0, 6.0, 161)

    def field(x, y):
        r2 = x**2 + y**2
        return x * np.exp(-0.5 * r2)

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=240,
        Ntheta=256,
        origin=(0.0, 0.0),
        recon_err_tol=1e-8,
    )

    assert dec.radial_power_tail_fraction == pytest.approx(1e-6)
    assert dec.radial_relative_amplitude_threshold == pytest.approx(1e-3)

    for m in (-1, 1):
        power_radius = dec.radial_power_support_radius[m]
        amplitude_radius = dec.radial_amplitude_support_radius[m]
        support_radius = dec.cutoff_radius[m]

        # The final diagnostic radius is explicitly the stricter criterion.
        assert support_radius == pytest.approx(max(power_radius, amplitude_radius))

        rho = np.abs(dec.rho[m])
        radial_power_density = rho**2 * dec.r
        total_power = trapezoid(radial_power_density, dec.r)
        mask = dec.r <= power_radius
        enclosed_power = trapezoid(radial_power_density[mask], dec.r[mask])
        # The power radius is interpolated between radial samples, so use a
        # small discretization allowance while checking the 1e-6 tail target.
        assert enclosed_power / total_power > 1.0 - 3e-6

        peak = np.max(rho)
        significant = np.flatnonzero(
            rho >= dec.radial_relative_amplitude_threshold * peak
        )
        assert significant.size > 0
        assert dec.r[significant[-1]] <= amplitude_radius + (dec.r[1] - dec.r[0])


def test_radial_support_parameter_validation_is_explicit():
    x = np.linspace(-2.0, 2.0, 41)
    y = np.linspace(-2.0, 2.0, 41)

    def field(x, y):
        return np.exp(-(x**2 + y**2))

    with pytest.raises(ValueError, match="radial_power_tail_fraction"):
        PolarDecomposition(field, x, y, radial_power_tail_fraction=1.0)
    with pytest.raises(ValueError, match="radial_relative_amplitude_threshold"):
        PolarDecomposition(field, x, y, radial_relative_amplitude_threshold=-1e-3)


def test_smooth_c3_example_is_regular_at_origin_and_has_expected_modes():
    x = np.linspace(-5.0, 5.0, 161)
    y = np.linspace(-5.0, 5.0, 161)

    def field(x, y):
        r2 = x**2 + y**2
        c3_harmonic = x**3 - 3.0 * x * y**2
        return np.exp(-0.45 * r2) * (1.0 + 0.08 * c3_harmonic)

    theta = np.linspace(0.0, 2.0 * np.pi, 121, endpoint=False)
    r_small = 1e-3
    values_at_origin = field(np.zeros_like(theta), np.zeros_like(theta))
    values_near_origin = field(
        r_small * np.cos(theta), r_small * np.sin(theta)
    )
    assert np.ptp(values_at_origin) == pytest.approx(0.0, abs=1e-15)
    # The angular variation is O(r^3), so at r=1e-3 it is O(1e-10).
    assert np.ptp(values_near_origin) < 2e-10

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=160,
        Ntheta=240,
        origin=(0.0, 0.0),
        recon_err_tol=1e-8,
    )
    assert dec.selected_pairs == [(0,), (3, -3)]
    assert set(dec.m_sorted) == {0, 3, -3}
