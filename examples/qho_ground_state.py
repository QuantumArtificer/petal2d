"""Isotropic 2D harmonic-oscillator ground-state probability density."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-8.0, 8.0, 161)
y = np.linspace(-8.0, 8.0, 161)


def density(x, y):
    return np.exp(-(x**2 + y**2)) / np.pi


dec = PolarDecomposition(
    density, x, y, Nr=192, Ntheta=256, recon_err_tol=1e-8
)
report_and_plot(dec, "2D QHO ground-state probability density")
