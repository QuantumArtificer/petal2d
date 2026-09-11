import numpy as np
import pytest

from petal2d import PolarDecomposition
import petal2d.polar as polar_module


def _grid(n=161, extent=6.0):
    x = np.linspace(-extent, extent, n)
    y = np.linspace(-extent, extent, n)
    return x, y


def test_isotropic_gaussian_is_pure_m0():
    x, y = _grid()

    def field(x, y):
        return np.exp(-(x**2 + y**2))

    dec = PolarDecomposition(
        field, x, y, Nr=128, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=1e-8
    )
    assert dec.m_sorted == [0]
    assert dec.target_reached
    assert dec.recon_error < 1e-8
    assert dec.recon_error_measured < 1e-8
    assert dec.domain_consistency == pytest.approx(1.0, rel=1e-3)


def test_single_complex_harmonic_is_identified_exactly():
    x, y = _grid()
    m0 = 4

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m0 * np.exp(-r**2) * np.exp(1j * m0 * theta)

    dec = PolarDecomposition(
        field, x, y, Nr=128, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=1e-8
    )
    assert not dec.is_real_input
    assert dec.m_sorted == [m0]
    assert dec.power_fracs[m0] == pytest.approx(1.0, abs=1e-12)
    assert dec.recon_error_measured < 1e-8


def test_real_cosine_has_equal_plus_minus_m_power():
    x, y = _grid()
    m0 = 3

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m0 * np.exp(-r**2) * np.cos(m0 * theta)

    dec = PolarDecomposition(
        field, x, y, Nr=128, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=1e-8
    )
    assert dec.is_real_input
    assert set(dec.m_sorted) == {-m0, m0}
    assert dec.power_fracs[m0] == pytest.approx(0.5, abs=1e-12)
    assert dec.power_fracs[-m0] == pytest.approx(0.5, abs=1e-12)
    assert dec.recon_error_measured < 1e-8
    assert dec.real_reconstruction_residual < 1e-13


def test_parseval_selection_matches_measured_weighted_error():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        g = np.exp(-0.7 * r**2)
        return g * (1.0 + 0.4 * np.exp(2j * theta) + 0.2 * np.exp(-3j * theta))

    dec = PolarDecomposition(
        field, x, y, Nr=160, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=20.0
    )
    assert dec.m_sorted == [0, 2]
    assert dec.power_fracs[0] == pytest.approx(1.0 / 1.2, rel=1e-10)
    assert dec.power_fracs[2] == pytest.approx(0.16 / 1.2, rel=1e-10)
    assert dec.recon_error == pytest.approx(100.0 * np.sqrt(0.04 / 1.2), rel=1e-10)
    assert dec.recon_error_measured == pytest.approx(dec.recon_error, rel=1e-10, abs=1e-10)
    assert dec.parseval_relative_error < 1e-12


def test_real_adaptive_selection_preserves_conjugate_pairs():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.4 * np.cos(3.0 * theta))

    dec = PolarDecomposition(
        field, x, y, Nr=128, Ntheta=256, origin=(0.0, 0.0), recon_err_tol=20.0
    )
    # Keeping only one of +/-3 would satisfy 20%, but PETAL2D deliberately
    # selects the conjugate pair together for a real input field.
    assert dec.selected_pairs == [(0,), (3, -3)]
    assert set(dec.m_sorted) == {0, 3, -3}
    assert dec.recon_error_measured < 1e-10
    assert dec.real_reconstruction_residual < 1e-13


def test_m_abs_max_reports_unreachable_target():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        g = np.exp(-r**2)
        return g * (0.1 + np.exp(5j * theta))

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=128,
        Ntheta=256,
        origin=(0.0, 0.0),
        m_abs_max=2,
        recon_err_tol=1.0,
    )
    assert set(dec.m_sorted) == {-2, -1, 0, 1, 2}
    assert not dec.target_reached
    assert dec.recon_error > 99.0
    assert dec.recon_error_measured == pytest.approx(dec.recon_error, rel=1e-10)


def test_centroid_default_and_explicit_origin():
    x = np.linspace(-5.0, 7.0, 241)
    y = np.linspace(-4.0, 8.0, 241)
    x0, y0 = 1.5, 2.0

    def field(x, y):
        return np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / 0.5)

    auto = PolarDecomposition(field, x, y, Nr=100, Ntheta=128, recon_err_tol=1e-6)
    explicit = PolarDecomposition(
        field, x, y, Nr=100, Ntheta=128, origin=np.array([x0, y0]), recon_err_tol=1e-6
    )

    assert auto.origin_mode == "centroid"
    assert auto.origin[0] == pytest.approx(x0, abs=2e-4)
    assert auto.origin[1] == pytest.approx(y0, abs=2e-4)
    assert explicit.origin_mode == "explicit"
    assert explicit.origin == pytest.approx((x0, y0))
    expected_safe = min(x0 - x[0], x[-1] - x0, y0 - y[0], y[-1] - y0)
    assert auto.safe_rmax == pytest.approx(expected_safe, abs=2e-4)
    assert auto.r[-1] == pytest.approx(auto.safe_rmax)
    assert auto.m_sorted == [0]


def test_domain_consistency_quantifies_excluded_weight():
    x, y = _grid(241, 6.0)

    def field(x, y):
        return np.exp(-0.5 * (x**2 + y**2))

    full = PolarDecomposition(field, x, y, Nr=180, Ntheta=256, origin=(0.0, 0.0))
    small = PolarDecomposition(
        field, x, y, Nr=180, Ntheta=256, origin=(0.0, 0.0), rmax=1.0
    )
    assert full.domain_consistency == pytest.approx(1.0, rel=3e-4)
    assert small.domain_consistency == pytest.approx(1.0 - np.exp(-1.0), rel=2e-3)
    assert small.domain_consistency < full.domain_consistency


def test_callable_path_bypasses_cartesian_interpolation(monkeypatch):
    x, y = _grid(81, 4.0)

    def fail_if_called(*args, **kwargs):
        raise AssertionError("map_coordinates should not be used for callable input")

    monkeypatch.setattr(polar_module, "map_coordinates", fail_if_called)

    def field(x, y):
        return np.exp(-(x**2 + y**2))

    dec = PolarDecomposition(field, x, y, Nr=64, Ntheta=128, origin=(0.0, 0.0))
    assert dec.m_sorted == [0]


def test_sampled_field_recovers_expected_modes_with_small_interpolation_leakage():
    x, y = _grid(101, 5.0)
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X, Y)
    theta = np.arctan2(Y, X)
    sampled = np.exp(-0.6 * r**2) * (1.0 + 0.2 * np.cos(3.0 * theta))

    dec = PolarDecomposition(
        sampled,
        x,
        y,
        Nr=128,
        Ntheta=256,
        origin=(0.0, 0.0),
        recon_err_tol=1.0,
        interp_method="cubic",
    )
    assert set(dec.m_sorted) == {0, -3, 3}
    assert dec.power_fracs[0] == pytest.approx(1.0 / 1.02, abs=5e-4)
    assert dec.power_fracs[3] == pytest.approx(0.01 / 1.02, abs=5e-4)
    assert dec.power_fracs[-3] == pytest.approx(0.01 / 1.02, abs=5e-4)
    assert 1.0 - dec.retained_power_fraction < 1e-4
    assert dec.recon_error_measured < 1.0


def test_reconstruct_and_full_spectrum_api():
    x, y = _grid()

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.3 * np.exp(2j * theta))

    dec = PolarDecomposition(
        field, x, y, Nr=100, Ntheta=128, origin=(0.0, 0.0), recon_err_tol=5.0
    )
    assert dec.m_all.shape == (128,)
    assert dec.rho_all.shape == (100, 128)
    assert dec.powers_all.shape == (128,)
    assert dec.power_fracs_all.sum() == pytest.approx(1.0)
    assert np.allclose(dec.reconstruct(), dec.f_recon)
    assert dec._weighted_relative_error_percent(dec.f_polar, dec.reconstruct("all")) < 1e-10

    m0_only = dec.reconstruct([0])
    assert m0_only.shape == dec.f_polar.shape
    with pytest.raises(ValueError, match="outside"):
        dec.reconstruct([1000])


def test_nyquist_limit_even_and_odd():
    x, y = _grid(121, 4.0)

    def even_field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * np.cos(8.0 * theta)

    even = PolarDecomposition(
        even_field, x, y, Nr=80, Ntheta=16, origin=(0.0, 0.0), m_abs_max=8, recon_err_tol=1e-8
    )
    assert even.m_nyquist == 8
    assert even.nyquist_mode == -8
    assert -8 in even.m_all and 8 not in even.m_all
    assert even.m_sorted == [-8]

    def odd_field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * np.cos(7.0 * theta)

    odd = PolarDecomposition(
        odd_field, x, y, Nr=80, Ntheta=15, origin=(0.0, 0.0), m_abs_max=7, recon_err_tol=1e-8
    )
    assert odd.m_nyquist == 7
    assert odd.nyquist_mode is None
    assert set(odd.m_sorted) == {-7, 7}

    with pytest.raises(ValueError, match="Nyquist"):
        PolarDecomposition(even_field, x, y, Nr=40, Ntheta=16, m_abs_max=9)


def test_normalization_and_input_validation():
    x, y = _grid(81, 3.0)
    dec = PolarDecomposition(
        lambda x, y: np.exp(-(x**2 + y**2)),
        x,
        y,
        Nr=64,
        normalize=True,
        origin=(0.0, 0.0),
    )
    assert dec.cartesian_power == pytest.approx(1.0, rel=1e-12)

    nonuniform = np.array([-2.0, -1.0, 0.2, 2.0])
    y2 = np.linspace(-2.0, 2.0, 5)
    with pytest.raises(ValueError, match="uniformly spaced"):
        PolarDecomposition(np.ones((4, 5)), nonuniform, y2, Nr=4)

    x2 = np.linspace(-2.0, 2.0, 5)
    with pytest.raises(ValueError, match="shape"):
        PolarDecomposition(np.ones((4, 5)), x2, y2, Nr=4)
    with pytest.raises(ValueError, match="non-zero finite L2 power"):
        PolarDecomposition(np.zeros((5, 5)), x2, y2, Nr=4)
    with pytest.raises(ValueError, match="length-two"):
        PolarDecomposition(np.ones((5, 5)), x2, y2, Nr=4, origin=(0.0,))
    with pytest.raises(ValueError, match="inside"):
        PolarDecomposition(np.ones((5, 5)), x2, y2, Nr=4, origin=(9.0, 0.0))


def test_spectrum_table_and_print_spectrum_are_dependency_free(capsys):
    x, y = _grid(101, 4.0)

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return np.exp(-r**2) * (1.0 + 0.2 * np.cos(2.0 * theta))

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=80,
        Ntheta=128,
        origin=(0.0, 0.0),
        recon_err_tol=1e-8,
    )
    rows = dec.spectrum_table(retained_only=True)
    assert isinstance(rows, list)
    assert {row["m"] for row in rows} == {0, -2, 2}
    expected_keys = {
        "m",
        "power",
        "power_fraction",
        "radial_power_support_radius",
        "radial_amplitude_support_radius",
        "cutoff_radius",
        "retained",
    }
    assert all(set(row) == expected_keys for row in rows)
    assert all(row["retained"] for row in rows)

    dec.print_spectrum(retained_only=True)
    out = capsys.readouterr().out
    assert "power" in out
    assert "fraction" in out
    assert "retained" in out
    assert "r_cutoff" in out
    assert "True" in out



def test_ranked_power_plot_uses_categorical_harmonic_ticks_for_sparse_high_m():
    """A sparse high-m spectrum must not create a huge linear m-axis."""
    import matplotlib.pyplot as plt

    x, y = _grid(121, 5.0)
    m0 = 40

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m0 * np.exp(-2.0 * r**2) * np.exp(1j * m0 * theta)

    dec = PolarDecomposition(
        field,
        x,
        y,
        Nr=96,
        Ntheta=128,
        origin=(0.0, 0.0),
        recon_err_tol=1e-8,
    )
    assert dec.m_sorted == [m0]

    fig, axs = dec.plot_harmonics_with_hist()
    labels = [tick.get_text().replace("$", "") for tick in axs[1].get_xticklabels()]
    assert labels == [str(m0)]
    assert axs[1].get_xlim()[1] - axs[1].get_xlim()[0] < 2.0
    plt.close(fig)
