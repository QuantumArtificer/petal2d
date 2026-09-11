# Radial cutoff diagnostics

The cutoff radius is designed to answer a practical question:

> Beyond what radius is a retained harmonic negligible according to both its integrated power and its pointwise radial amplitude?

It is not used to truncate the PETAL2D reconstruction.

## Power-support radius

For a retained harmonic define

$$
P_m(r)=\int_0^r|\rho_m(r')|^2r'\,dr'.
$$

Given `radial_power_tail_fraction = ε_P`, the power-support radius is the smallest radius satisfying

$$
P_m(r_{m}^{\rm power})\ge(1-\varepsilon_P)P_m(R).
$$

This controls integrated omitted power.

## Amplitude-support radius

Integrated power alone can place a radius inside a visually or physically meaningful low-amplitude outer lobe. Therefore PETAL2D also uses `radial_relative_amplitude_threshold = ε_A` and defines

$$
r_m^{\rm amp}=\sup\left\{r:\;|\rho_m(r)|\ge
\varepsilon_A\max_{r'}|\rho_m(r')|\right\}.
$$

The outermost crossing is used so oscillatory profiles cannot discard a later significant lobe.

## Final cutoff radius

$$
r_{m}^{\rm cutoff}=\max(r_m^{\rm power},r_m^{\rm amp}).
$$

Consequently, beyond the cutoff radius both conditions hold: the remaining integrated power is sufficiently small and no sampled outer region exceeds the chosen relative amplitude threshold.

The defaults are

```python
radial_power_tail_fraction = 1e-6
radial_relative_amplitude_threshold = 1e-3
```

corresponding to 99.9999% enclosed radial mode power and a 0.1% relative amplitude threshold.


:::{important}
`cutoff_radius` is a radial-information diagnostic, not a claim that the field becomes exactly zero there. If a downstream calculation needs a finite radial domain, choosing a domain at or beyond the cutoff ensures that both configured tail criteria have been met. If the field itself must be forced to zero, use a smooth taper/window outside the physically relevant region rather than multiplying by a discontinuous hard step. A sharp truncation can introduce artificial kinks and high-frequency content.
:::

## Why two criteria are necessary

For a normalized Gaussian-like tail, a radius containing 99% of the integrated power can still occur where the pointwise amplitude is around ten percent of its peak. Increasing the power percentage alone is not a universal solution because the relation between integrated weight and visible amplitude depends on radial shape. The combined criterion is explicitly designed to protect both notions of support.
