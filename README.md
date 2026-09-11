# PETAL2D

**PETAL2D** — **P**olar **E**xpansion **T**oolkit for **A**tomic Orbitals and **L**ocalized Fields in **2D** — decomposes localized two-dimensional scalar fields into angular Fourier harmonics while retaining their radial structure.

```text
2D field f(x,y)
      ↓ choose/estimate origin
polar field f(r,θ)
      ↓ angular FFT at every r
radial harmonics ρₘ(r)
      ↓ Parseval-weighted power ranking
interpretable spectrum + controlled reconstruction
```

For

$$
f(r,\theta)=\sum_m\rho_m(r)e^{im\theta},
$$

PETAL2D computes the radial harmonic profiles $\rho_m(r)$, their polar-area-weighted powers

$$
P_m=\int|\rho_m(r)|^2r\,dr,
$$

and an adaptive reconstruction with an exact discrete Parseval error certificate.

## Why PETAL2D?

A Cartesian 2D FFT decomposes a field into plane waves. PETAL2D answers a different question: **what angular structure lives around a chosen localization center, and where does each angular channel live radially?**

That makes the representation natural for localized states whose shape is better described as $s$-, $p$-, $d$-, vortex-, or crystal-field-like than as a collection of Cartesian plane waves.

### Highlights

- analytic callables or sampled Cartesian arrays;
- automatic $L^2$ power-centroid origin or an explicit physical center;
- direct callable evaluation on the polar grid;
- cubic or linear interpolation for sampled arrays;
- complete discrete angular spectrum plus adaptive retained spectrum;
- real-valued adaptive truncation that preserves $(+m,-m)$ conjugate pairs;
- exact Parseval-based reconstruction-error prediction for the computed spectrum;
- public reconstruction of retained, complete, or user-selected harmonics;
- angular Nyquist and aliasing diagnostics;
- `domain_consistency` to quantify whether the polar disk captures the Cartesian-domain field weight;
- explicit radial power-support, amplitude-support, and final `cutoff_radius` diagnostics;
- publication-oriented plotting, mathematical tests, JSON validation, and reproducible benchmarks.

## Installation

From PyPI after the first release:

```bash
python -m pip install petal2d
```

For development:

```bash
python -m pip install -e ".[test,docs]"
```

## 60-second start

```python
import numpy as np
from petal2d import PolarDecomposition

x = np.linspace(-6, 6, 161)
y = np.linspace(-6, 6, 161)

def field(x, y):
    return x * np.exp(-0.5 * (x**2 + y**2))   # p_x-like orbital

dec = PolarDecomposition(field, x, y, recon_err_tol=0.1)

dec.print_spectrum()
print("origin:", dec.origin)
print("selected pairs:", dec.selected_pairs)
print("domain consistency:", dec.domain_consistency)
print("reconstruction error (%):", dec.recon_error_measured)
```

For a real $p_x$-like field, PETAL2D identifies the conjugate pair $m=\pm1$ and keeps the adaptive reconstruction real.

Plot the radial profiles, power-ranked retained harmonics, and reconstruction:

```python
import matplotlib.pyplot as plt

dec.plot_harmonics_with_hist("p_x-like orbital")
dec.plot_original_vs_reconstructions("p_x-like orbital")
plt.show()
```

## Inputs

A callable:

```python
def psi(x, y):
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return r**4 * np.exp(-0.7*r**2) * np.exp(4j*theta)

dec = PolarDecomposition(psi, x, y)
```

or a sampled array:

```python
X, Y = np.meshgrid(x, y, indexing="ij")
f_xy = np.exp(-(X**2 + Y**2))
dec = PolarDecomposition(f_xy, x, y, interp_method="cubic")
```

The sampled array must have shape `(len(x), len(y))` on strictly increasing, uniformly spaced Cartesian axes.

## What can PETAL2D be used for?

PETAL2D is useful whenever a **localized 2D scalar or complex scalar field has meaningful angular structure about a center**. Potential applications include:

- **Atomic and molecular orbital analysis in planar models.** Resolve $s$-, $p$-, $d$-like and mixed angular content while keeping radial profiles.
- **Wannier and localized crystal orbitals.** Quantify local angular character, crystal-field mixing, and symmetry-compatible harmonics around a chosen atomic or moiré site.
- **Defect, impurity, and bound-state wavefunctions.** Compare angular channels and localization radii between states or control parameters.
- **Moiré electronic states.** Characterize local orbital textures around high-symmetry stacking regions without reducing the result to a single shape label.
- **Complex order parameters and vortices.** Analyze phase-winding channels of a complex scalar field; a pure $e^{im\theta}$ vortex appears as a single complex harmonic.
- **Charge, probability, or scalar spin-density anisotropy.** Measure spatial angular structure of real observables. For vector spin textures, decompose components separately.
- **Excitonic relative-coordinate wavefunctions.** Resolve angular channels in an effectively two-dimensional electron-hole relative coordinate.
- **Localized photonic, phononic, acoustic, or PDE eigenmodes.** Decompose scalar mode profiles centered on a resonator, defect, or localization center.
- **Scientific imaging and simulation descriptors.** Compress or compare localized angular morphology when the field can be interpreted as a scalar intensity or signed scalar map.

The harmonic index must be interpreted according to the input. For a complex wavefunction, $m$ is an $L_z/\hbar$ angular basis label. For a density $|\psi|^2$, it describes density anisotropy and cannot recover phase information discarded by the modulus square.

## Key diagnostics

### Expansion origin

The default is the $L^2$ power centroid,

$$
\mathbf r_c=\frac{\int \mathbf r|f|^2d^2r}{\int|f|^2d^2r},
$$

or supply a known center directly:

```python
dec = PolarDecomposition(f_xy, x, y, origin=(x0, y0))
```

### Domain consistency

$$
\mathrm{domain\_consistency}=
\frac{\int_{\Omega_{\rm polar}}|f|^2d^2r}
{\int_{\Omega_{\rm Cartesian}}|f|^2d^2r}.
$$

A value near one means the analyzed polar disk contains essentially the same $L^2$ weight as the supplied Cartesian domain. It measures captured weight, not pointwise interpolation accuracy.

### Adaptive reconstruction

```python
f_adaptive = dec.reconstruct()
f_all = dec.reconstruct("all")
f_selected = dec.reconstruct([0, 3, -3])

print(dec.reconstruction_error())
```

For real input, automatic selection keeps conjugate pairs together. Explicit user-selected subsets are reconstructed exactly as requested.

### Radial cutoff

For every retained $m$:

```python
dec.radial_power_support_radius[m]
dec.radial_amplitude_support_radius[m]
dec.cutoff_radius[m]
```

The final cutoff radius is the larger of an integrated-power criterion and an outermost relative-amplitude criterion. It is a diagnostic of radial extent; **PETAL2D never truncates the reconstruction at this radius.**

## Documentation

The documentation is organized like a Scientific Python project rather than a long README:

- [**Getting started**](docs/source/getting_started.md) — first decomposition and core workflow;
- [**User guide**](docs/source/user_guide/index.md) — inputs, centering, parameter choice, diagnostics, applications, plotting, and numerical failure modes;
- [**Theory**](docs/source/theory/index.md) — Fourier conventions, physical interpretation, Nyquist structure, and radial cutoff;
- [**Full proofs**](docs/source/theory/proofs.md) — Parseval, adaptive error certificate, real-field conjugacy, rotation/translation covariance, $C_n$ selection rules, aliasing, and domain consistency;
- [**Guided examples**](docs/source/examples/index.md) — QHO, pure angular momentum, real $p/d$ orbitals, $C_3$ field, automatic centering, sampled data, adaptive compression, and a Wannier workflow;
- [**Validation and benchmarks**](docs/source/validation/index.md) — accuracy, convergence, runtime, and memory methodology;
- [**API reference**](docs/source/reference/index.rst) — exact public class and method documentation.

Build locally with:

```bash
python -m pip install -e ".[docs]"
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

## Mathematical guarantees and numerical caveats

For the computed discrete angular spectrum, orthogonality gives the adaptive reconstruction error directly from omitted power. PETAL2D also evaluates the weighted error independently and the validation suite checks their agreement.

Finite Cartesian domains, finite `Nr`, finite `Ntheta`, and Cartesian-to-polar interpolation remain numerical approximations. Weak harmonics in sampled data should therefore be convergence-tested before being interpreted as physical symmetry breaking.

## Validation and benchmarks

Run mathematical tests:

```bash
python -m pytest -v
```

Reproduce numerical validation and performance studies:

```bash
./tools/run_full_validation_and_benchmarks.sh
```

The compute stages write JSON first; plotting scripts read the saved JSON. Benchmark metadata records package/dependency versions and relevant hardware/runtime information.

## Examples

```bash
python examples/qho_ground_state.py
python examples/angular_momentum_state.py
python examples/px_orbital.py
python examples/dx2_y2_orbital.py
python examples/c3_localized_field.py
python examples/shifted_orbital_centroid.py
```

## Project scope

PETAL2D currently targets single-center, two-dimensional scalar fields on uniform Cartesian grids. It does not currently provide 3D spherical harmonics, coupled vector/tensor decomposition, nonuniform-grid interpolation, or automatic crystallographic point-group classification.

## Citation

PETAL2D includes `CITATION.cff`. A versioned Zenodo DOI will be added with the first archived release. If PETAL2D contributes materially to a publication, cite the archived software version used for the calculation.

## Contributing

Bug reports, mathematical test cases, documentation improvements, physically motivated examples, and carefully validated numerical improvements are welcome. See `docs/source/development/contributing.md` and `CONTRIBUTING.md`.

## License

MIT License.

Developed by **Alex Santacruz**, 2DQMAT Research @ IF-UNAM.
