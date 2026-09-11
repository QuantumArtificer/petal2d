"""Shared presentation helpers for PETAL2D examples."""

from matplotlib import pyplot as plt


def report_and_plot(dec, title: str, extra: dict | None = None) -> None:
    """Print a consistent summary and show both standard PETAL2D figures."""
    print(f"\n=== {title} ===")
    print(f"origin: {dec.origin}")
    print(f"domain_consistency: {dec.domain_consistency:.8f}")
    print(f"selected_pairs: {dec.selected_pairs}")
    print(f"retained_modes: {dec.m_sorted}")
    print(f"radial_power_tail_fraction: {dec.radial_power_tail_fraction:.3e}")
    print(
        "radial_relative_amplitude_threshold: "
        f"{dec.radial_relative_amplitude_threshold:.3e}"
    )
    print(f"reconstruction_error_percent: {dec.recon_error_measured:.6g}")
    print(f"target_reached: {dec.target_reached}")
    if extra:
        for key, value in extra.items():
            print(f"{key}: {value}")
    dec.print_spectrum()
    dec.plot_harmonics_with_hist(title)
    dec.plot_original_vs_reconstructions(title)
    plt.show()
