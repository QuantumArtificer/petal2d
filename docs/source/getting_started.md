# Getting started


## Installation

For a release from PyPI:

```bash
python -m pip install petal2d
```

For a development checkout:

```bash
python -m pip install -e ".[test,docs]"
```

Verify the import:

```python
import petal2d
print(petal2d.__version__)
```


```text
0.1.0
```

## The central idea

Suppose $f(x,y)$ is localized around a point $(x_0,y_0)$. Relative to that point,

$$
x-x_0=r\cos\theta,\qquad
y-y_0=r\sin\theta.
$$

At every radius, PETAL2D expands the angular dependence as

$$
f(r,\theta)=\sum_m\rho_m(r)e^{im\theta}.
$$

The integer $m$ labels an angular harmonic and $\rho_m(r)$ tells you where radially that angular structure lives. If the analyzed field is a complex wavefunction, the $e^{im\theta}$ channels are eigenchannels of the planar angular-momentum operator $L_z$ about the chosen origin. Their integrated power fractions therefore resolve the wavefunction's $L_z$-sector weights on the analyzed disk. For a density or another real scalar observable, the same $m$ values describe angular shape rather than hidden wavefunction phase.

## First decomposition: an isotropic Gaussian

```python
import numpy as np
from petal2d import PolarDecomposition

x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)

def field(x, y):
    return np.exp(-(x**2 + y**2))

dec = PolarDecomposition(field, x, y)
```

With the defaults PETAL2D will:

1. evaluate the field on the supplied Cartesian domain.
2. determine the $L^2$ power centroid.
3. choose the largest centered polar disk fully contained in the rectangle.
4. evaluate a callable directly on the polar grid, or interpolate sampled arrays.
5. perform the angular FFT.
6. compute radial mode powers.
7. retain the smallest power-ranked set meeting the requested reconstruction error.
8. compute radial cutoff diagnostics for every retained harmonic.

Inspect the result:

```python
print("origin:", dec.origin)
print("domain consistency:", dec.domain_consistency)
print("selected pairs:", dec.selected_pairs)
print("reconstruction error (%):", dec.recon_error_measured)
dec.print_spectrum()
```


```text
origin: (3.533949646070574e-17, 5.478458685715792e-19)
domain consistency: 0.9995311180756926
selected pairs: [(0,)]
reconstruction error (%): 1.266726935300768e-14
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    2.498828e-01  1.000000e+00           2.63027               2.62854       2.63027      True
```

The tiny nonzero coordinates are floating-point noise around the exact origin $(0,0)$. The physically meaningful result is that the isotropic field is entirely in the $m=0$ channel to numerical precision.

## Plot the decomposition

```python
import matplotlib.pyplot as plt

dec.plot_harmonics_with_hist("Isotropic Gaussian")
dec.plot_original_vs_reconstructions("Isotropic Gaussian")
plt.show()
```


```{figure} _static/examples/getting_started_spectrum.png
:width: 90%
:alt: Angular spectrum and radial profile for the getting-started Gaussian

The retained radial profile and power-ranked angular spectrum. The dashed line marks `cutoff_radius` for the retained mode.
```

```{figure} _static/examples/getting_started_reconstruction.png
:width: 90%
:alt: Original and reconstructed getting-started Gaussian

Original field and retained-mode reconstruction using the same color normalization.
```

The cutoff radius is a diagnostic, not a numerical truncation of the reconstruction. It is the larger of two independently computed radii: one that leaves no more than the requested radial power tail, and one that lies beyond the last sampled radius whose amplitude exceeds the configured relative-amplitude threshold. See {doc}`theory/radial_cutoff` for the exact definitions.

## Reconstruct the field

```python
f_adaptive = dec.reconstruct()
f_all = dec.reconstruct("all")
f_m0 = dec.reconstruct([0])

print(dec.reconstruction_error())
print(dec.reconstruction_error("all"))
print(dec.reconstruction_error([0]))
```


```text
1.266726935300768e-14
1.406744239744267e-14
1.266726935300768e-14
```

The returned error is the polar-area-weighted relative $L^2$ error in percent. Here all three reconstructions are effectively exact because the field contains only $m=0$.

## Sampled arrays instead of callables

PETAL2D also accepts a two-dimensional array sampled on the Cartesian grid:

```python
X, Y = np.meshgrid(x, y, indexing="ij")
sampled = np.exp(-(X**2 + Y**2))

dec_sampled = PolarDecomposition(
    sampled,
    x,
    y,
    interp_method="cubic",
)

print(dec_sampled.selected_pairs)
print(dec_sampled.recon_error_measured)
```


```text
[(0,)]
4.6017946300227716e-05
```

The sampled version retains the same physically expected $m=0$ channel. Its small nonzero reconstruction error is controlled by interpolation and finite-grid effects rather than direct callable evaluation.

The array shape must be `(len(x), len(y))`, and `x` and `y` must be strictly increasing and uniformly spaced. Cubic interpolation is the default because validation shows a large accuracy advantage for a modest runtime cost. A more demanding sampled-data example with measured interpolation leakage is developed in {doc}`examples/sampled_field`.

## What to check after every decomposition

A useful first-pass checklist is:

- `target_reached`: did the admissible angular spectrum satisfy `recon_err_tol`?
- `selected_pairs`: which adaptive angular units were retained?
- `recon_error_measured`: what weighted reconstruction error was measured directly?
- `domain_consistency`: does the analyzed polar disk contain essentially the same $L^2$ weight as the supplied Cartesian domain?
- `print_spectrum()`: what are the retained mode powers and radial cutoff diagnostics?

For the Gaussian run above:

```python
print(dec.target_reached)
print(dec.selected_pairs)
print(dec.recon_error_measured)
print(dec.domain_consistency)
```


```text
True
[(0,)]
1.266726935300768e-14
0.9995311180756926
```

Continue with {doc}`user_guide/choosing_parameters` for parameter selection, {doc}`examples/index` for progressively more demanding examples, and {doc}`theory/index` for the derivation.
