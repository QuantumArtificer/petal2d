import numpy as np

from petal2d import PolarDecomposition


def test_preparation_coordinate_meshes_are_not_retained():
    axis = np.linspace(-4.0, 4.0, 81)

    def field(x, y):
        return np.exp(-(x * x + y * y))

    dec = PolarDecomposition(
        field,
        axis,
        axis,
        Nr=96,
        Ntheta=128,
        origin=(0.0, 0.0),
    )

    # These were formerly large persistent work arrays. PETAL2D needs only the
    # sampled polar field and spectral arrays after preparation.
    for name in ("_X", "_Y", "R", "TH", "Xp", "Yp", "_f"):
        assert not hasattr(dec, name)

    assert dec.f_polar.shape == (96, 128)
    assert dec.rho_all.shape == (96, 128)


def test_sampled_input_still_works_without_retained_work_meshes():
    axis = np.linspace(-4.0, 4.0, 81)
    X, Y = np.meshgrid(axis, axis, indexing="ij")
    sampled = np.exp(-(X * X + Y * Y)) * (1.0 + 0.1 * (X**3 - 3 * X * Y**2))

    dec = PolarDecomposition(
        sampled,
        axis,
        axis,
        Nr=96,
        Ntheta=128,
        origin=(0.0, 0.0),
        interp_method="cubic",
        recon_err_tol=1.0,
    )

    for name in ("_X", "_Y", "R", "TH", "Xp", "Yp", "_f"):
        assert not hasattr(dec, name)

    assert dec.domain_consistency > 0.99
    assert dec.reconstruction_error() <= 1.0
