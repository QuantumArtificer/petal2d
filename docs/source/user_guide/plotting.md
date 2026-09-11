# Plotting

PETAL2D provides two plotting helpers. They use a local Matplotlib `rc_context`, so they do not modify the user's global plotting configuration.


## Harmonic spectrum and radial profiles

```python
fig, axes = dec.plot_harmonics_with_hist("My localized state")
```

The left panel shows $|\rho_m(r)|$ for retained harmonics, with color-matched dashed lines marking `cutoff_radius[m]`. The right panel is categorical and power-ranked: only retained harmonics appear, ordered from greatest to least fractional power. A pure $m=100$ state therefore produces one bar labelled `100` rather than an empty linear axis spanning hundreds of unused integers.

The right-panel title is Retained harmonics (power ranked) and the x-axis is Angular harmonic $m$.

```{figure} ../_static/examples/m4_spectrum.png
:width: 88%
:alt: PETAL2D radial profiles and ranked angular-power plot

Pure complex $m=4$ example.
```

## Original versus reconstruction

```python
fig, axes = dec.plot_original_vs_reconstructions("My localized state")
```

For real sign-changing fields PETAL2D uses a zero-centered diverging color scale and preserves the sign. For nonnegative real fields it uses a sequential scale. For complex fields it displays magnitudes. Original and reconstruction always share one color normalization.

```{figure} ../_static/examples/px_reconstruction.png
:width: 88%
:alt: Signed px original and reconstruction

Real $p_x$ example with the signed field preserved.
```

By default the displayed window is based on the largest retained cutoff radius, with modest padding. This changes the viewport only. It does not crop the analyzed field or alter reconstruction. An explicit display radius can be requested:

```python
dec.plot_original_vs_reconstructions(
    "My localized state",
    display_radius=5.0,
)
```

## Saving figures

The helpers return ordinary Matplotlib figure and axes objects:

```python
fig, _ = dec.plot_harmonics_with_hist("My state")
fig.savefig("spectrum.pdf", dpi=300)
```

For publication, vector output such as PDF or SVG is preferable for line and bar graphics. The field panels are rasterized internally within vector figures so that large polar meshes do not create unnecessarily large PDF files.
