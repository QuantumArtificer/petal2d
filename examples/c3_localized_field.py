"""Localized real field with C3 angular modulation: m=0 and m=+/-3."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-7.0, 7.0, 181)
y = np.linspace(-7.0, 7.0, 181)


def c3_field(x, y):
    """Smooth C3-symmetric field with m=0 and m=+/-3 content.

    The anisotropic term is proportional to
    r^3 cos(3 theta) = x^3 - 3 x y^2, so it vanishes smoothly at r=0.
    """
    r2 = x**2 + y**2
    c3_harmonic = x**3 - 3.0 * x * y**2
    return np.exp(-0.45 * r2) * (1.0 + 0.08 * c3_harmonic)


dec = PolarDecomposition(
    c3_field, x, y, Nr=180, Ntheta=300, origin=(0.0, 0.0), recon_err_tol=1e-8
)
report_and_plot(dec, "Smooth C3-symmetric localized field")
