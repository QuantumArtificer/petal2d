# Scope and limitations

PETAL2D is designed for localized two-dimensional scalar fields around a single expansion center.

## Supported problem class

PETAL2D is designed for

- two-dimensional real or complex scalar fields
- localized structure around one meaningful expansion center
- uniformly spaced Cartesian coordinates for sampled arrays
- polar angular-Fourier analysis with retained radial profiles.

## Not currently provided

PETAL2D does not currently provide

- 3D spherical harmonics
- coupled vector/tensor harmonic decomposition
- interpolation from nonuniform Cartesian grids
- automatic crystallographic point-group classification
- automated multi-center decomposition
- boundary-condition-aware periodic decomposition of an entire unit cell
- inference of an unaliased continuum harmonic beyond the angular Nyquist limit
- a guarantee that every tiny nonzero FFT bin is physical.

## Interpretation limits

Angular-harmonic labels are basis labels about the chosen origin. For a complex wavefunction they coincide with eigenchannels of the planar $L_z$ operator about that origin, but in a crystal a dominant $m$ should not automatically be interpreted as an exact isolated-atom orbital quantum number: the lattice can mix angular channels allowed by the local symmetry.

For a density such as $|\psi|^2$, PETAL2D analyzes the angular morphology of the density. It cannot recover wavefunction phase or angular-momentum information that was removed before the density was supplied.

Likewise, a spectrum compatible with a $C_n$ rotational selection rule is evidence about angular content around the chosen center. By itself it is not a complete proof of the full point-group symmetry of the underlying physical system.
