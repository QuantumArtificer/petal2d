# Interpreting results

This page uses the sampled $C_3$ field from {doc}`../examples/sampled_field` as a concrete reference.

## Full spectrum versus retained spectrum

PETAL2D computes the complete discrete angular FFT before adaptive truncation. The complete spectrum is available as

```python
dec.m_all
dec.rho_all
dec.powers_all
dec.power_fracs_all
```

For the sampled $C_3$ example:

```python
print("m_all shape:", dec.m_all.shape)
print("rho_all shape:", dec.rho_all.shape)
print("powers_all shape:", dec.powers_all.shape)
print("power_fracs_all shape:", dec.power_fracs_all.shape)
```


```text
m_all shape: (512,)
rho_all shape: (121, 512)
powers_all shape: (512,)
power_fracs_all shape: (512,)
```

The adaptively retained subset is available through

```python
dec.m_sorted
dec.rho
dec.power_fracs
dec.selected_pairs
```

and the same run gives

```text
m_sorted: [0, 3, -3]
selected_pairs: [(0,), (3, -3)]
rho keys: [0, 3, -3]
power_fracs: {0: 0.9956769952656367, 3: 0.0021615023671450488, -3: 0.0021615023671450496}
```

The full spectrum includes floating-point and interpolation-level leakage. The retained spectrum is the compressed representation chosen by `recon_err_tol`.

## Radial coefficients `rho[m]`

`dec.rho[m]` is $\rho_m(r)$ sampled on `dec.r`. It is generally complex even when the original field is real. For a real field,

$$
\rho_{-m}(r)=\rho_m(r)^*,
$$

so the two channels carry equal integrated power.

## Mode power

PETAL2D defines

$$
P_m=\int_0^{r_{\max}}|\rho_m(r)|^2r\,dr.
$$

`power_fracs_all` is `powers_all` normalized by $\sum_mP_m$. The factor $2\pi$ relating this convention to the physical polar $L^2$ norm cancels in every fraction.

A large power fraction means that a harmonic accounts for a large fraction of the field's polar $L^2$ weight. Radial location is described by `rho[m]` and the cutoff diagnostics.

## `selected_pairs`

For real input, `selected_pairs` lists the adaptive units actually retained. The sampled $C_3$ run gives

```text
[(0,), (3, -3)]
```

The singleton `(0,)` is selected independently. The nonzero real-field harmonics are selected as the conjugate pair `(3, -3)` so that automatic adaptive reconstruction preserves real-valuedness. For genuinely complex input, each selected channel is a singleton.

## Reconstruction error

`recon_error` is the Parseval-predicted error from omitted spectral power. `recon_error_measured` independently evaluates the reconstructed field in the weighted polar norm.

```python
print(dec.recon_error)
print(dec.recon_error_measured)
```


```text
2.7069303579009995e-05
2.7016812777500244e-05
```

Both values are percentages. Their agreement is an internal consistency check. Small differences can appear for sampled data because the two errors are evaluated through different numerical routes.

## `domain_consistency`

PETAL2D defines

$$
C_D=
\frac{\int_{\Omega_{\rm polar}}|f|^2d^2r}
{\int_{\Omega_{\rm Cartesian}}|f|^2d^2r}.
$$

For the sampled $C_3$ run,

```text
domain_consistency = 0.9996539035424497
```

Interpretation:

- $C_D\approx1$: the polar disk and Cartesian rectangle contain essentially the same field weight
- $C_D<1$: appreciable weight lies in parts of the rectangle outside the analyzed disk, or numerical discretization makes the polar estimate smaller
- $C_D>1$: possible for callable input with an explicitly enlarged polar domain, or from small quadrature/interpolation overshoot.

A useful limiting example is a uniform field on a square: the largest inscribed disk contains area fraction $\pi/4$, so its exact domain consistency is $\pi/4$, not one. A localized Gaussian can have `domain_consistency` very close to one even though the disk occupies less geometric area, because the omitted corners contain negligible field weight.

`domain_consistency` is a weight diagnostic, not a pointwise interpolation-error metric.

## Cutoff radii

For every retained mode PETAL2D reports

```python
dec.radial_power_support_radius[m]
dec.radial_amplitude_support_radius[m]
dec.cutoff_radius[m]
```

with

$$
r_{\rm cutoff,m}
=
\max(r_{\rm power,m},r_{\rm amplitude,m}).
$$

For the sampled $C_3$ example, the actual values are

```text
m= 0: r_power=3.395936536178625, r_amplitude=3.3938767167529695, r_cutoff=3.395936536178625
m= 3: r_power=4.220679790622508, r_amplitude=4.370270544732886,  r_cutoff=4.370270544732886
m=-3: r_power=4.22067979060771,  r_amplitude=4.370270544732886,  r_cutoff=4.370270544732886
```

`cutoff_radius` describes where the remaining sampled radial tail is negligible according to both configured criteria. It never changes `f_recon`. Reconstruction is performed on the complete analyzed polar grid.
