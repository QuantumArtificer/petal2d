# Discrete angular sampling and Nyquist structure

PETAL2D samples

$$
\theta_j=\frac{2\pi j}{N_\theta},\qquad j=0,\ldots,N_\theta-1,
$$

and computes a discrete Fourier transform.

## Representable FFT indices

For odd $N_\theta$, the canonical integer indices are

$$
-\frac{N_\theta-1}{2},\ldots,0,\ldots,+\frac{N_\theta-1}{2}.
$$

For even $N_\theta$, NumPy's FFT convention contains

$$
0,1,\ldots,\frac{N_\theta}{2}-1,
-\frac{N_\theta}{2},\ldots,-1.
$$

The Nyquist harmonic is represented by the single index $-N_\theta/2$. The index $+N_\theta/2$ gives the same sampled sequence and is not a second independent bin.

## Aliasing

For sampled angles,

$$
e^{i(m+kN_\theta)\theta_j}=e^{im\theta_j}
$$

for any integer $k$. Thus the data determine $m$ only modulo $N_\theta$. PETAL2D's `m_all` reports the canonical FFT representative of each discrete class.

## Real fields

For real angular samples the DFT satisfies discrete conjugate symmetry. Away from $m=0$ and the even-grid Nyquist singleton, bins occur in $(+m,-m)$ pairs. PETAL2D's adaptive selector preserves those pairs automatically.

## Choosing `Ntheta`

The strict Nyquist condition is only the beginning. In practical sampled-data problems, interpolation and finite Cartesian resolution can generate weak higher harmonics. Use enough angular samples that the physically relevant range is comfortably below Nyquist, then verify convergence by increasing `Ntheta`.
