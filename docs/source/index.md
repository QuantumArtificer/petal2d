# PETAL2D

Polar Expansion Toolkit for Atomic Orbitals and Localized Fields in 2D

PETAL2D is a lightweight, general-purpose toolkit for decomposing and projecting localized two-dimensional scalar fields onto polar angular harmonics, the natural circular-harmonic basis for angular structure in 2D:

$$
f(r,\theta)=\sum_{m\in\mathbb Z}\rho_m(r)e^{im\theta}.
$$

PETAL2D accepts either analytical callables or data sampled on a uniform Cartesian grid. It preserves the radial profiles $\rho_m(r)$, ranks angular channels by their integrated spectral weight, and can construct a controlled adaptive reconstruction from the smallest retained set that satisfies a requested weighted-$L^2$ error.

::::{grid} 2
:::{grid-item-card} Start here
:link: getting_started
:link-type: doc
Install PETAL2D, decompose your first field, inspect the spectrum, and reconstruct it.
:::
:::{grid-item-card} User guide
:link: user_guide/index
:link-type: doc
Learn the input conventions, parameters, diagnostics, numerical behavior, and recommended workflow.
:::
:::{grid-item-card} API reference
:link: reference/index
:link-type: doc
Detailed documentation for the public class, methods, parameters, and returned quantities.
:::
:::{grid-item-card} Guided examples
:link: examples/index
:link-type: doc
Move from elementary angular harmonics to controlled compression, sampled data, real orbitals, symmetry, and centering.
:::
:::{grid-item-card} Validation and benchmarks
:link: validation/index
:link-type: doc
See how interpolation, radial quadrature, angular sampling, runtime, and memory behave in controlled tests.
:::
:::{grid-item-card} Theory and full proofs
:link: theory/index
:link-type: doc
Derive the expansion, Parseval identities, error certificate, real-field pairing, symmetry rules, sampling limits, and cutoff diagnostics.
:::
::::

## What PETAL2D is for

PETAL2D is useful whenever a localized real or complex 2D scalar field has meaningful angular structure about a chosen center. Typical applications include atomic-orbital-like wavefunctions, Wannier and other localized crystal orbitals, impurity or defect states, charge-density anisotropies, complex order parameters and vortices, excitonic relative-coordinate wavefunctions, localized photonic or phononic modes, and angular descriptors of scientific images or simulation fields.

The same decomposition is also useful as a compact angular basis in problems with approximate rotational structure: for example, diagnosing local $C_n$ selection rules, separating symmetry-compatible angular channels, constructing reduced descriptions of localized PDE eigenmodes, or identifying the dominant angular content of a field before a problem-specific multipole or symmetry analysis.

PETAL2D does not replace a Cartesian FFT, a 3D spherical-harmonic expansion, a vector or tensor harmonic decomposition, or a full crystallographic point-group analysis. See {doc}`reference/limitations` for the package scope.

## Release archive

PETAL2D 0.1.0 is archived on Zenodo with DOI [10.5281/zenodo.22715741](https://doi.org/10.5281/zenodo.22715741).

```{toctree}
:maxdepth: 2
:hidden:

getting_started
user_guide/index
reference/index
examples/index
validation/index
theory/index
development/index
```
