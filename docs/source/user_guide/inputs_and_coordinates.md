# Inputs, coordinates, and centering

PETAL2D accepts the same physical field in two forms: an analytical callable or a sampled Cartesian array. The distinction matters numerically because callables can be evaluated directly on the polar grid, while sampled arrays require interpolation.


## Callable fields

A callable must accept array-valued `x` and `y` and return either a scalar or an array with the broadcasted mesh shape. PETAL2D evaluates it on the Cartesian mesh to determine normalization and centering, then evaluates it directly on the polar grid.

```python
import numpy as np
from petal2d import PolarDecomposition

x = np.linspace(-5.0, 5.0, 121)
y = np.linspace(-5.0, 5.0, 121)

def psi(x, y):
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    return r**2 * np.exp(-r**2) * np.exp(2j * theta)

dec_callable = PolarDecomposition(
    psi,
    x,
    y,
    origin=(0.0, 0.0),
    Ntheta=256,
)

print("callable selected:", dec_callable.selected_pairs)
```


```text
callable selected: [(2,)]
```

Because the callable is evaluated directly at $(r\cos\theta,r\sin\theta)$, this path avoids Cartesian-to-polar interpolation error.

## Sampled Cartesian fields

For sampled data, pass an array with shape `(len(x), len(y))`. PETAL2D uses `scipy.ndimage.map_coordinates` to interpolate the array onto the polar grid. The coordinates must be uniformly spaced because the interpolation indices are constructed from constant Cartesian spacings.

Using samples of the same $m=2$ field:

```python
X, Y = np.meshgrid(x, y, indexing="ij")
f_xy = psi(X, Y)

dec_sampled = PolarDecomposition(
    f_xy,
    x,
    y,
    origin=(0.0, 0.0),
    Ntheta=256,
    interp_method="cubic",
)

print("sampled selected:", dec_sampled.selected_pairs)
print(
    "sampled reconstruction error (%):",
    dec_sampled.recon_error_measured,
)
```


```text
sampled selected: [(2,)]
sampled reconstruction error (%): 0.0002737698248805121
```

The physically correct channel is still recovered, but sampled input carries a finite interpolation error. See {doc}`../examples/sampled_field` and {doc}`numerical_considerations` for convergence guidance.

## Expansion origin

By default,

```python
origin="centroid"
```

uses the $L^2$ power centroid

$$
\mathbf r_c=
\frac{\int \mathbf r\,|f(\mathbf r)|^2\,d^2r}
{\int |f(\mathbf r)|^2\,d^2r}.
$$

This definition is phase independent and applies to real, sign-changing, and complex fields. Its translation covariance is derived in {doc}`../theory/proofs`.

You may instead specify a physical center directly:

```python
dec_explicit = PolarDecomposition(
    f_xy,
    x,
    y,
    origin=(0.25, -0.10),
    Ntheta=256,
)

print(dec_explicit.origin)
```


```text
(0.25, -0.1)
```

Use an explicit origin when the physics supplies a distinguished site or symmetry center, or when the field contains several separated localized objects whose global power centroid would lie between them.

## `rmax` and `safe_rmax`

If `rmax=None`, PETAL2D uses the largest disk centered on the chosen origin that is fully contained in the Cartesian rectangle:

$$
r_{\mathrm{safe}}
=
\min(
x_0-x_{\min},
x_{\max}-x_0,
y_0-y_{\min},
y_{\max}-y_0
).
$$

This avoids angular sectors leaving the sampled domain.

For sampled arrays, requesting a larger `rmax` causes points outside the rectangle to be filled with zero by the interpolation rule. That discontinuity can create artificial angular structure near the boundary, so enlarging `rmax` beyond `safe_rmax` should be deliberate and convergence-tested.

For callables, an explicit larger `rmax` is mathematically well defined because the function can be evaluated outside the original Cartesian rectangle. In that case `domain_consistency` may legitimately exceed one because the polar domain contains physical area not included in the Cartesian reference domain.
