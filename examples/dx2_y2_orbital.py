"""Real d_(x^2-y^2)-like orbital with a conjugate m=+/-2 pair."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)


def d_x2_y2(x, y):
    return (x**2 - y**2) * np.exp(-0.5 * (x**2 + y**2))


dec = PolarDecomposition(
    d_x2_y2, x, y, Nr=160, Ntheta=256, recon_err_tol=1e-8
)
report_and_plot(dec, r"d_(x^2-y^2)-like orbital")
