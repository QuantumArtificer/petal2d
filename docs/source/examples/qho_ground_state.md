# Example 1: isotropic QHO ground-state density

The two-dimensional isotropic harmonic-oscillator ground-state wavefunction is rotationally symmetric. Its probability density may be written, up to units and normalization, as

$$
n(r)=\frac{1}{\pi}e^{-r^2}.
$$

There is no angular dependence, so theory predicts

$$
\rho_0(r)=n(r),\qquad \rho_{m\ne0}(r)=0.
$$

Run:

```bash
python examples/qho_ground_state.py
```


```text
=== 2D QHO ground-state probability density ===
domain_consistency: 0.99941502
selected_pairs: [(0,)]
retained_modes: [0]
reconstruction_error_percent: 1.40994e-14
target_reached: True
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     0    2.531548e-02  1.000000e+00           2.63121               2.62905       2.63121      True
```

The $m=0$ result also checks radial integration, reconstruction, and the cutoff-radius diagnostic for a smooth monotone tail.

:::{note}
Because the input here is a density, the spectrum describes density anisotropy. An excited wavefunction $\psi\propto e^{im\theta}$ can carry nonzero angular momentum while $|\psi|^2$ remains isotropic.
:::

```{figure} ../_static/examples/qho_spectrum.png
:width: 88%
:alt: QHO radial harmonic spectrum

The isotropic density is entirely $m=0$ within numerical precision.
```

```{figure} ../_static/examples/qho_reconstruction.png
:width: 88%
:alt: QHO original and reconstruction

Original field and retained-mode reconstruction on a common color scale.
```
