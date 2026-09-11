import numpy as np

from petal2d import PolarDecomposition


def _field(x, y):
    return np.exp(-(x * x + y * y))


def test_default_radial_resolution_matches_smaller_cartesian_dimension():
    x = np.linspace(-4.0, 4.0, 81)
    y = np.linspace(-3.0, 3.0, 61)
    dec = PolarDecomposition(_field, x, y, Ntheta=96, origin=(0.0, 0.0))
    assert dec.r.size == 61


def test_explicit_radial_resolution_overrides_default():
    axis = np.linspace(-4.0, 4.0, 81)
    dec = PolarDecomposition(
        _field, axis, axis, Nr=37, Ntheta=96, origin=(0.0, 0.0)
    )
    assert dec.r.size == 37
