# Example 5: real $p_x$ and $d_{x^2-y^2}$ orbitals

Real orbitals illustrate why PETAL2D treats real and complex fields differently. A real field obeys

$$
\rho_{-m}(r)=\rho_m(r)^*,
$$

so nonzero harmonics occur in conjugate pairs.

## $p_x$: an $m=\pm1$ pair

Using

$$
x=r\cos\theta
=\frac{r}{2}\left(e^{i\theta}+e^{-i\theta}\right),
$$

a localized $p_x$-like orbital contains equal $m=+1$ and $m=-1$ power. PETAL2D therefore treats $(1,-1)$ as one adaptive selection unit rather than keeping either member alone and creating an artificial complex reconstruction.

Run:

```bash
python examples/px_orbital.py
```


```text
=== p_x-like orbital ===
domain_consistency: 1.00000003
selected_pairs: [(1, -1)]
retained_modes: [1, -1]
reconstruction_error_percent: 2.70394e-14
target_reached: True
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
    -1    1.250000e-01  5.000000e-01           4.08711                4.2065        4.2065      True
     1    1.250000e-01  5.000000e-01           4.08711                4.2065        4.2065      True
```

```{figure} ../_static/examples/px_spectrum.png
:width: 88%
:alt: px conjugate harmonic pair

The real $p_x$ orbital has equal $m=\pm1$ power.
```

```{figure} ../_static/examples/px_reconstruction.png
:width: 88%
:alt: px original and reconstruction

Conjugate-pair selection preserves the signed real orbital.
```

## $d_{x^2-y^2}$: an $m=\pm2$ pair

Similarly,

$$
x^2-y^2=r^2\cos2\theta
=\frac{r^2}{2}\left(e^{i2\theta}+e^{-i2\theta}\right),
$$

so a real $d_{x^2-y^2}$-like orbital contains equal $m=\pm2$ power.

```bash
python examples/dx2_y2_orbital.py
```


```text
=== d_(x^2-y^2)-like orbital ===
domain_consistency: 1.00000000
selected_pairs: [(2, -2)]
retained_modes: [2, -2]
reconstruction_error_percent: 4.43516e-14
target_reached: True
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
    -2    2.500000e-01  5.000000e-01           4.37494               4.52437       4.52437      True
     2    2.500000e-01  5.000000e-01           4.37494               4.52437       4.52437      True
```

```{figure} ../_static/examples/d2_spectrum.png
:width: 88%
:alt: d x2-y2 conjugate harmonic pair

The $d_{x^2-y^2}$ orbital has equal $m=\pm2$ power.
```

```{figure} ../_static/examples/d2_reconstruction.png
:width: 88%
:alt: d x2-y2 original and reconstruction

Signed original field and retained-mode reconstruction.
```


