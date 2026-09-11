# Numerical considerations and failure modes

## Exact mathematics versus numerical implementation

The Fourier identities in the theory section are exact statements for square-integrable fields. A PETAL2D calculation additionally contains finite-domain, radial-quadrature, angular-sampling, and, when arrays are supplied, Cartesian interpolation errors.

PETAL2D exposes diagnostics so these effects are visible rather than hidden.

## Radial quadrature

Mode powers are evaluated using trapezoidal integration on the radial grid. Validation on smooth analytic profiles shows the expected approximately second-order convergence with radial spacing. If mode powers or cutoff radii are more important than runtime, increase `Nr`.

## Angular aliasing

Sampling $N_\theta$ equally spaced angles identifies harmonics only modulo $N_\theta$. A continuum harmonic above the Nyquist range aliases into a representable discrete FFT index. PETAL2D reports the discrete spectrum correctly. It cannot infer the unsampled continuum frequency from aliased data.

See {doc}`../theory/sampling`.

## Interpolation leakage

A sampled Cartesian field with exact angular content can acquire small spurious angular channels when interpolated to polar coordinates. This is numerical leakage, not necessarily physical symmetry breaking. Increasing Cartesian resolution, using cubic interpolation, and checking convergence help distinguish the two.

Do not automatically suppress weak modes with `m_abs_max` before checking whether they converge away.

## Centering sensitivity

Angular decompositions are origin dependent. Even a perfectly symmetric field develops multiple angular harmonics when expanded about the wrong point. `origin="centroid"` minimizes user effort and is translation covariant for localized single-center fields, but a known physical center is preferable when available.

## Multiple localized objects

If a field contains two separated objects, the global power centroid may lie between them. PETAL2D will then describe angular structure about that midpoint. For site-resolved interpretation, window the field around one object or provide the corresponding explicit origin.

## Complex fields

PETAL2D preserves independent $m$ channels for complex input. Real-field pair selection should not be expected for a wavefunction with genuine phase winding.

## What PETAL2D does not currently support

- nonuniform Cartesian coordinate grids
- vector or tensor fields as a single decomposition object
- 3D spherical harmonics
- automatic point-group identification
- automatic separation of multiple localization centers
- inference of an unaliased continuum harmonic above the angular Nyquist limit.
