# Example 6: a smooth $C_3$-symmetric localized field

PETAL2D is naturally suited to analyzing planar rotational symmetry around a chosen center. Consider

$$
f(x,y)=e^{-\alpha(x^2+y^2)}
\left[1+\beta(x^3-3xy^2)\right].
$$

Because

$$
x^3-3xy^2=r^3\cos3\theta,
$$

the anisotropic contribution is smooth at the origin and contains only the conjugate angular pair $m=\pm3$.

Run:

```bash
python examples/c3_localized_field.py
```


```text
=== Smooth C3-symmetric localized field ===
domain_consistency: 0.99977646
selected_pairs: [(0,), (3, -3)]
retained_modes: [0, 3, -3]
reconstruction_error_percent: 2.4606e-14
target_reached: True
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    5.554281e-01  9.743327e-01           3.91976               3.91837       3.91976      True
     3    7.315958e-03  1.283366e-02           4.87289               5.04603       5.04603      True
    -3    7.315958e-03  1.283366e-02           4.87289               5.04603       5.04603      True
```

This is the simplest explicit illustration of the $C_n$ selection rule: an exactly $C_3$-invariant scalar field can contain only harmonics with $m=0\pmod 3$.

In a crystal calculation, weak $m=\pm6,\pm9,\ldots$ channels may encode additional $C_3$-compatible angular structure rather than symmetry breaking. By contrast, robust converged $m=\pm1$ or $\pm2$ content would be incompatible with exact $C_3$ invariance about the chosen center.

```{figure} ../_static/examples/c3_spectrum.png
:width: 88%
:alt: Smooth C3 harmonic spectrum

The spectrum contains $m=0$ and the real conjugate pair $m=\pm3$.
```

```{figure} ../_static/examples/c3_reconstruction.png
:width: 88%
:alt: Smooth C3 original and reconstruction

Original field and retained-mode reconstruction.
```
