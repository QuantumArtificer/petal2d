import matplotlib.pyplot as plt
import numpy as np

from petal2d import PolarDecomposition


def test_ranked_angular_power_panel_has_concise_label_and_title():
    axis = np.linspace(-5.0, 5.0, 121)
    m0 = 40

    def field(x, y):
        r = np.hypot(x, y)
        theta = np.arctan2(y, x)
        return r**m0 * np.exp(-2.0 * r**2) * np.exp(1j * m0 * theta)

    dec = PolarDecomposition(
        field,
        axis,
        axis,
        Nr=96,
        Ntheta=128,
        origin=(0.0, 0.0),
        recon_err_tol=1e-8,
    )
    fig, axs = dec.plot_harmonics_with_hist()
    assert axs[1].get_xlabel() == r"Angular harmonic $m$"
    assert axs[1].get_title() == "Retained harmonics (power ranked)"
    labels = [tick.get_text().replace("$", "") for tick in axs[1].get_xticklabels()]
    assert labels == [str(m0)]
    plt.close(fig)
