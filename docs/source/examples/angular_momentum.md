# Example 2: a pure angular-momentum state

Consider a separable localized wavefunction

$$
\psi_m(r,\theta)=r^{|m|}e^{-\alpha r^2}e^{im\theta}.
$$

The factor $r^{|m|}$ makes the field regular at the origin, while the angular factor is an exact eigenfunction of the planar angular-momentum operator:

$$
L_z\psi_m=m\hbar\psi_m.
$$

The executable example uses $m=4$:

```bash
python examples/angular_momentum_state.py
```


```text
=== Complex localized angular-momentum state (m=4) ===
domain_consistency: 1.00000000
selected_pairs: [(4,)]
retained_modes: [4]
reconstruction_error_percent: 8.55582e-14
target_reached: True
     m           power      fraction   r_power_support   r_amplitude_support      r_cutoff  retained
----------------------------------------------------------------------------------------------------
     4    2.231213e+00  1.000000e+00           4.09379               4.24121       4.24121      True
```

PETAL2D therefore returns a single complex channel `(4,)`. The power-ranked harmonic panel contains one categorical bar labelled `4`. It does not create an empty linear axis from $-4$ to $4$. The same presentation remains compact even for a pure high-order harmonic such as $m=100$.

Rotating the state by $\phi$ multiplies $\rho_4$ by $e^{-4i\phi}$ but leaves $P_4$ unchanged. Coefficient phase therefore carries orientation information, whereas harmonic power does not.

```{figure} ../_static/examples/m4_spectrum.png
:width: 88%
:alt: Pure m=4 harmonic spectrum

The $m=4$ channel carries all angular power.
```

```{figure} ../_static/examples/m4_reconstruction.png
:width: 88%
:alt: Pure m=4 original and reconstruction

Magnitude of the original complex field and its reconstruction.
```
