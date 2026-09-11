# Example 4: sampled numerical data

Real applications often provide arrays rather than analytic functions. PETAL2D therefore accepts scalar data sampled on a uniform Cartesian grid and interpolates that data to the polar analysis grid before performing the angular decomposition.

The example uses a smooth $C_3$ field whose exact analytic angular content is $m=0,\pm3$, but supplies PETAL2D only with the sampled array:

```python
import numpy as np
from petal2d import PolarDecomposition

x = np.linspace(-5, 5, 121)
y = np.linspace(-5, 5, 121)
X, Y = np.meshgrid(x, y, indexing="ij")

r2 = X**2 + Y**2
field = np.exp(-0.6*r2) * (1 + 0.05*(X**3 - 3*X*Y**2))

dec = PolarDecomposition(
    field,
    x,
    y,
    interp_method="cubic",
    recon_err_tol=0.1,
)

dec.print_spectrum()
print("domain consistency:", dec.domain_consistency)
```

The complete repository example additionally measures power outside the analytically expected set:

```bash
python examples/sampled_field.py
```


```text
=== Sampled smooth C3 field ===
domain_consistency: 0.99965390
selected_pairs: [(0,), (3, -3)]
retained_modes: [0, 3, -3]
reconstruction_error_percent: 2.70168e-05
target_reached: True
power outside expected m={0,+/-3}: 7.299e-14
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    4.165218e-01  9.956770e-01           3.39594               3.39388       3.39594      True
    -3    9.042219e-04  2.161502e-03           4.22068               4.37027       4.37027      True
     3    9.042219e-04  2.161502e-03           4.22068               4.37027       4.37027      True
```

The interpolation leakage in this run is

```python
rows = dec.spectrum_table(retained_only=False)
leakage = sum(
    row["power_fraction"]
    for row in rows
    if row["m"] not in {0, 3, -3}
)
print(leakage)
```

with output

```text
7.299081726523838e-14
```

Weak symmetry-incompatible channels in sampled data should be tested for convergence before they are interpreted physically. The complete spectrum exposes interpolation-level leakage.

```{figure} ../_static/examples/sampled_spectrum.png
:width: 88%
:alt: Spectrum of sampled smooth C3 field

The retained spectrum contains only the expected $m=0,\pm3$ channels.
```

```{figure} ../_static/examples/sampled_reconstruction.png
:width: 88%
:alt: Original and reconstructed sampled C3 field

Original sampled field and retained-mode reconstruction.
```
