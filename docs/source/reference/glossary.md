# Glossary

```{glossary}
Angular harmonic
  An integer Fourier channel $m$ multiplying $e^{im\theta}$.

Mode power
  $P_m=\int|\rho_m(r)|^2r\,dr$. PETAL2D power fractions are normalized by $\sum_mP_m$.

Power centroid
  The first spatial moment of $|f|^2$, used by `origin="centroid"`.

Domain consistency
  Ratio of physical $L^2$ weight in the analyzed polar domain to that in the supplied Cartesian domain.

Selected pair
  An adaptive real-field selection unit $(m,-m)$. $m=0$ and the even-grid Nyquist mode are singletons. Complex fields use singleton channels.

Nyquist mode
  The highest discrete angular frequency class. For even `Ntheta`, NumPy represents it by `-Ntheta/2`.

Power-support radius
  Smallest radius enclosing the configured fraction of integrated radial mode power.

Amplitude-support radius
  Outermost radius where the mode amplitude remains above the configured fraction of its peak.

Cutoff radius
  The larger of the power-support and amplitude-support radii. It marks a radius beyond which both configured radial-tail criteria are satisfied. It does not truncate the reconstruction.

Interpolation leakage
  Spurious angular spectral weight introduced when a sampled Cartesian field is interpolated to a polar grid.
```
