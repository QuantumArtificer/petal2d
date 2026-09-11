# Polar harmonic expansion

## Function space and convention

Let $f(r,\theta)$ be a complex scalar field on a disk $0\le r\le R$, $0\le\theta<2\pi$, with finite polar $L^2$ norm

$$
\|f\|_{L^2(D_R)}^2=
\int_0^R\int_0^{2\pi}|f(r,\theta)|^2\,r\,d\theta\,dr<\infty.
$$

For almost every fixed $r$, the angular function belongs to $L^2(S^1)$ and therefore has a Fourier expansion

$$
f(r,\theta)=\sum_{m\in\mathbb Z}\rho_m(r)e^{im\theta},
$$

with PETAL2D's continuum normalization

$$
\rho_m(r)=\frac{1}{2\pi}\int_0^{2\pi}
f(r,\theta)e^{-im\theta}\,d\theta.
$$

This is the continuum convention approximated by `fft(f_polar)/Ntheta`.

## Physical meaning of $m$

The angular momentum operator in two-dimensional polar coordinates is

$$
L_z=-i\hbar\frac{\partial}{\partial\theta}.
$$

Because

$$
L_z e^{im\theta}=m\hbar e^{im\theta},
$$

the $e^{im\theta}$ channels are angular-momentum eigenfunctions. If $f$ is a complex wavefunction, $\rho_m(r)$ resolves its $L_z$ content while retaining radial information.

For a density or other real scalar observable, the same basis is still mathematically natural, but $m$ labels angular shape harmonics, not necessarily a quantum number of an underlying state.

For a normalized complex wavefunction, the fraction $P_m/\sum_nP_n$ is the probability weight in the $L_z=m\hbar$ angular sector after the radial degree of freedom is summed over. This interpretation is proved explicitly in {doc}`proofs`.

## Mode power

Define the radial power of mode $m$ as

$$
P_m=\int_0^R|\rho_m(r)|^2r\,dr.
$$

Angular Parseval gives

$$
\frac{1}{2\pi}\int_0^{2\pi}|f(r,\theta)|^2d\theta
=\sum_m|\rho_m(r)|^2.
$$

Integrating over $r$ yields

$$
\|f\|_{L^2(D_R)}^2=2\pi\sum_mP_m.
$$

Thus the fractional power

$$
w_m=\frac{P_m}{\sum_n P_n}
$$

is a natural dimensionless measure of angular content.

## Why retain radial profiles?

A single number $P_m$ says how important a harmonic is globally. The function $\rho_m(r)$ says where it lives. Two fields can have identical angular power fractions but very different radial structures, nodes, shell locations, or localization lengths. PETAL2D keeps both levels of information.

## Point-group fingerprints

If a scalar field is exactly invariant under a $2\pi/n$ rotation,

$$
f(r,\theta+2\pi/n)=f(r,\theta),
$$

then only harmonics with $m$ divisible by $n$ may be nonzero. More generally, if a complex field transforms with a one-dimensional rotational character

$$
f(r,\theta+2\pi/n)=e^{i\ell 2\pi/n}f(r,\theta),
$$

then nonzero harmonics satisfy

$$
m\equiv \ell \pmod n.
$$

This makes PETAL2D useful for diagnosing local crystal-field mixing and rotational symmetry around a site. The formal proof is in {doc}`proofs`.

:::{warning}
A spectrum compatible with a point group is evidence of angular structure, not by itself a proof of full spatial symmetry. Radial dependence, the chosen origin, additional reflections, layer/orbital labels, and other degrees of freedom may matter.
:::
