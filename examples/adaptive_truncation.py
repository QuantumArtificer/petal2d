"""Smooth mixed-harmonic field illustrating Parseval-controlled truncation."""

import numpy as np

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)


def field(x, y):
    """Smooth field with exact angular channels m=0, +2, and -4."""
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    g = r**4 * np.exp(-0.7 * r**2)
    return g * (
        1.0
        + 0.4 * np.exp(2j * theta)
        + 0.2 * np.exp(-4j * theta)
    )


dec = PolarDecomposition(
    field,
    x,
    y,
    Ntheta=256,
    recon_err_tol=20.0,
    origin=(0.0, 0.0),
)

report_and_plot(
    dec,
    "Adaptive truncation of a smooth known mixed spectrum",
    extra={
        "Parseval-predicted error (%)": f"{dec.recon_error:.6f}",
        "directly measured error (%)": f"{dec.recon_error_measured:.6f}",
    },
)
