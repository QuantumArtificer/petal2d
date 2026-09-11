# Choosing parameters

PETAL2D's defaults are chosen for general-purpose accuracy without requiring users to tune every grid. This page explains what each numerical parameter controls and when it should be changed.

## Radial resolution: `Nr`

If omitted,

```python
Nr = min(len(x), len(y))
```

The validation suite found this to be a good accuracy/cost compromise. Halving the radial resolution increases radial-quadrature errors by about a factor of four, consistent with second-order trapezoidal convergence, while doubling the default gives another roughly factor-of-four improvement at increased runtime and memory cost.

Increase `Nr` when:

- radial nodes or narrow shells are important
- cutoff radii must be located more accurately
- interpolation error is already small and radial quadrature is the dominant error.

Decrease `Nr` only when speed or memory matters more than precise radial powers and cutoff diagnostics.

## Angular resolution: `Ntheta`

The angular grid is uniform on $[0,2\pi)$. Its discrete Fourier spectrum contains integer modes up to the FFT Nyquist limit.

A practical rule is to choose `Ntheta` comfortably larger than twice the largest $|m|$ you expect. The Nyquist condition sets the representability limit. Additional angular samples provide margin for mixtures of modes and sampled data.

PETAL2D exposes

```python
print(dec.m_nyquist)
print(dec.nyquist_mode)
```

For `Ntheta=512`:

```text
256
-256
```

For even `Ntheta`, NumPy represents the Nyquist class by the single index `-Ntheta/2`. There is no independent `+Ntheta/2` bin. PETAL2D rejects `m_abs_max` values above the representable absolute Nyquist limit.

## Angular truncation: `recon_err_tol`

`recon_err_tol` is the requested maximum relative weighted-$L^2$ reconstruction error in percent. PETAL2D ranks admissible angular selection units by power and retains the smallest number needed to meet that target whenever the allowed spectrum makes it possible.

For real fields, nonzero selection units are conjugate pairs $(m,-m)$, except for $m=0$ and the even-grid Nyquist singleton. For complex fields, channels are selected independently.

In {doc}`../examples/adaptive_truncation`, `recon_err_tol=20.0` retains two of three known complex channels and gives a measured error of approximately 18.2574%.

If `target_reached` is false, the requested target cannot be met with the allowed spectrum. A restrictive `m_abs_max` is the most common reason.

## Maximum angular index: `m_abs_max`

Use `m_abs_max` when physics or a controlled low-angular-momentum model tells you that only a bounded set of harmonics should be considered. Do not use it merely to hide weak numerical leakage: leakage is diagnostically useful because convergence with grid resolution can distinguish numerical artifacts from physical angular content.

## Interpolation: `interp_method`

For sampled arrays:

- `"cubic"` is the default and is strongly preferred for accuracy
- `"linear"` is faster and can be useful for exploratory scans or exceptionally large inputs.

Validation on smooth analytic fields showed orders-of-magnitude lower harmonic leakage for cubic interpolation at a modest runtime penalty. Callable inputs bypass Cartesian-to-polar interpolation because they are evaluated directly on the polar grid.

## Radial cutoff diagnostics

Two independent parameters define when each retained radial harmonic is negligible:

- `radial_power_tail_fraction=1e-6`: at most $10^{-6}$ of the integrated mode power may remain outside the power-support radius
- `radial_relative_amplitude_threshold=1e-3`: the amplitude-support radius is the outermost sampled radius where $|\rho_m|$ remains at least $10^{-3}$ of its peak.

The final

```python
dec.cutoff_radius[m]
```

is the larger of these two radii. It is a diagnostic of radial extent and does not truncate the reconstruction. See {doc}`../theory/radial_cutoff` for the exact definitions and proof.
