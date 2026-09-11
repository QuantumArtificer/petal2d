"""Shifted localized field demonstrating automatic power-centroid centering."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x0, y0 = 1.4, -0.9
x = np.linspace(-6.0, 6.0, 181)
y = np.linspace(-6.0, 6.0, 181)


def shifted(x, y):
    return np.exp(-0.8 * ((x - x0) ** 2 + (y - y0) ** 2))


dec = PolarDecomposition(
    shifted, x, y, Nr=150, Ntheta=240, recon_err_tol=1e-8
)
report_and_plot(
    dec,
    "Shifted Gaussian with automatic power-centroid origin",
    extra={"expected_origin": (x0, y0)},
)
