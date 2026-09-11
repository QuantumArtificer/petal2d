"""Complex localized angular-momentum state with a single m=4 harmonic."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

m0 = 4
x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)


def orbital(x, y):
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return r**m0 * np.exp(-0.7 * r**2) * np.exp(1j * m0 * theta)


dec = PolarDecomposition(
    orbital, x, y, Nr=160, Ntheta=256, recon_err_tol=1e-8
)
report_and_plot(dec, f"Complex localized angular-momentum state (m={m0})")
