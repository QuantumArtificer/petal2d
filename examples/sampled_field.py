"""Sampled Cartesian C3 field illustrating interpolation and leakage diagnostics."""

import numpy as np
from matplotlib import pyplot as plt

from petal2d import PolarDecomposition
from _common import report_and_plot

x = np.linspace(-5.0, 5.0, 121)
y = np.linspace(-5.0, 5.0, 121)
X, Y = np.meshgrid(x, y, indexing="ij")
r2 = X**2 + Y**2
sampled = np.exp(-0.6 * r2) * (1.0 + 0.05 * (X**3 - 3.0 * X * Y**2))

dec = PolarDecomposition(
    sampled,
    x,
    y,
    interp_method="cubic",
    recon_err_tol=0.1,
)

rows = dec.spectrum_table(retained_only=False)
leakage = sum(
    row["power_fraction"] for row in rows if row["m"] not in {0, 3, -3}
)

report_and_plot(
    dec,
    "Sampled smooth C3 field",
    extra={"power outside expected m={0,+/-3}": f"{leakage:.3e}"},
)
