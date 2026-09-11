# Example 3: controlled adaptive compression

To test adaptive truncation cleanly, use a field with a known, smooth angular spectrum:

$$
f(r,\theta)=g(r)\left[
1+0.4e^{2i\theta}+0.2e^{-4i\theta}
\right],
\qquad
g(r)=r^4e^{-0.7r^2}.
$$

The common $r^4$ factor makes every angular contribution regular at the origin. In Cartesian language, each term is a smooth polynomial times a Gaussian:

- $r^4$ is $(x^2+y^2)^2$
- $r^4e^{2i\theta}=r^2(x+iy)^2$
- $r^4e^{-4i\theta}=(x-iy)^4$.

Because all three channels share the same radial envelope, their powers are proportional to

$$
1^2=1,\qquad 0.4^2=0.16,\qquad 0.2^2=0.04.
$$

The normalized power fractions are therefore

$$
\frac{1}{1.20},\qquad
\frac{0.16}{1.20},\qquad
\frac{0.04}{1.20}.
$$

If only the first two channels are kept, the exact relative $L^2$ error is

$$
\sqrt{\frac{0.04}{1.20}}\approx18.2574\%.
$$

A complete minimal version is

```python
import numpy as np
from petal2d import PolarDecomposition

x = np.linspace(-6.0, 6.0, 161)
y = np.linspace(-6.0, 6.0, 161)

def field(x, y):
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    g = r**4 * np.exp(-0.7 * r**2)
    return g * (
        1.0
        + 0.4 * np.exp(2j * theta)
        + 0.2 * np.exp(-4j * theta)
    )

dec = PolarDecomposition(
    field,
    x,
    y,
    Ntheta=256,
    recon_err_tol=20.0,
    origin=(0.0, 0.0),
)

print(dec.m_sorted)
print(dec.recon_error)
print(dec.recon_error_measured)
```


```text
[0, 2]
18.257418583505537
18.257418583505537
```

Run the complete example:

```bash
python examples/adaptive_truncation.py
```


```text
=== Adaptive truncation of a smooth known mixed spectrum ===
domain_consistency: 1.00000000
selected_pairs: [(0,), (2,)]
retained_modes: [0, 2]
reconstruction_error_percent: 18.2574
target_reached: True
Parseval-predicted error (%): 18.257419
directly measured error (%): 18.257419
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    2.231213e+00  8.333333e-01           4.09288               4.24065       4.24065      True
     2    3.569941e-01  1.333333e-01           4.09288               4.24065       4.24065      True
```

The $m=-4$ channel is omitted because the first two channels already meet the requested 20% error. The predicted and directly measured errors agree to numerical precision. This is not a heuristic compression score: for the computed discrete spectrum it follows directly from orthogonality and Parseval's identity.

```{figure} ../_static/examples/adaptive_spectrum.png
:width: 88%
:alt: Adaptively retained smooth mixed angular spectrum

Only the two channels required by the requested error are retained.
```

```{figure} ../_static/examples/adaptive_reconstruction.png
:width: 88%
:alt: Adaptive reconstruction of the smooth mixed spectrum

The visible difference is the omitted $m=-4$ contribution.
```
