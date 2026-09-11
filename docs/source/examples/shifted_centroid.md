# Example 7: automatic centering

An angular decomposition is origin dependent: translating a field away from the expansion center generally produces additional angular harmonics. PETAL2D therefore defaults to the $L^2$ power centroid when no physical center is supplied explicitly.

The example translates a Gaussian to

$$
(x_0,y_0)=(1.4,-0.9).
$$

Run:

```bash
python examples/shifted_orbital_centroid.py
```


```text
=== Shifted Gaussian with automatic power-centroid origin ===
origin: (1.3999999999999995, -0.9)
domain_consistency: 0.99974580
selected_pairs: [(0,)]
retained_modes: [0]
reconstruction_error_percent: 4.45883e-14
target_reached: True
expected_origin: (1.4, -0.9)
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    3.124206e-01  1.000000e+00           2.93995                2.9388       2.93995      True
```

The recovered origin agrees with the imposed displacement to numerical precision. Once recentered, the translated Gaussian is again purely $m=0$.

```{figure} ../_static/examples/shifted_reconstruction.png
:width: 88%
:alt: Shifted Gaussian original and reconstruction

Original field and retained-mode reconstruction after power-centroid centering.
```

For atomistic or other site-resolved data, prefer an explicit center when one is physically known. Automatic centering is most useful when the localization center is unknown or has drifted numerically.
