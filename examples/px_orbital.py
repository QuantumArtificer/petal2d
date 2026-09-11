"""Real p_x-like orbital with a conjugate m=+/-1 pair."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)


def px(x, y):
    return x * np.exp(-0.5 * (x**2 + y**2))


dec = PolarDecomposition(
    px, x, y, Nr=160, Ntheta=256, recon_err_tol=1e-8
)
report_and_plot(dec, "p_x-like orbital")
