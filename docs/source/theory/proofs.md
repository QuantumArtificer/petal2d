# Proofs

## Polar Fourier expansion and Parseval identity

Let

$$
D_R=\{(r,\theta):0\le r\le R,\ 0\le\theta<2\pi\}
$$

and let $f\in L^2(D_R,r\,dr\,d\theta)$. Define

$$
\rho_m(r)=\frac{1}{2\pi}\int_0^{2\pi}
f(r,\theta)e^{-im\theta}\,d\theta.
$$

For almost every fixed $r$, the function $f(r,\cdot)$ belongs to $L^2(S^1)$. The Fourier basis is complete, so

$$
f(r,\theta)=\sum_{m\in\mathbb Z}\rho_m(r)e^{im\theta}
$$

in angular $L^2$, and

$$
\frac{1}{2\pi}\int_0^{2\pi}|f(r,\theta)|^2\,d\theta
=
\sum_m|\rho_m(r)|^2.
$$

### Proof

Since

$$
\int_0^R
\left[
\int_0^{2\pi}|f(r,\theta)|^2\,d\theta
\right]r\,dr<\infty,
$$

Fubini's theorem implies that the inner integral is finite for almost every $r$. Thus $f(r,\cdot)\in L^2(S^1)$ for almost every $r$.

The functions

$$
\frac{e^{im\theta}}{\sqrt{2\pi}},
\qquad m\in\mathbb Z,
$$

form a complete orthonormal basis of $L^2(S^1)$. The coefficients $\rho_m(r)$ are the corresponding Fourier coefficients. Completeness gives angular $L^2$ convergence, and Parseval's identity gives the pointwise-in-$r$ norm relation above.

Multiplying by $2\pi r$ and integrating radially gives

$$
\int_0^R\int_0^{2\pi}|f|^2r\,d\theta\,dr
=
2\pi\sum_m
\int_0^R|\rho_m(r)|^2r\,dr.
$$

Define

$$
P_m=\int_0^R|\rho_m(r)|^2r\,dr.
$$

Then the total polar power is $2\pi\sum_mP_m$.

## Exact reconstruction error

Let $S\subset\mathbb Z$ be a set of retained harmonics and define

$$
f_S(r,\theta)=
\sum_{m\in S}\rho_m(r)e^{im\theta}.
$$

Then

$$
\frac{\|f-f_S\|_{L^2(D_R)}}{\|f\|_{L^2(D_R)}}
=
\sqrt{
1-
\frac{\sum_{m\in S}P_m}
{\sum_mP_m}
}.
$$

### Proof

The omitted field is

$$
f-f_S=
\sum_{m\notin S}\rho_m(r)e^{im\theta}.
$$

Orthogonality of distinct angular harmonics gives

$$
\|f-f_S\|^2
=
2\pi\sum_{m\notin S}P_m,
$$

while Parseval gives

$$
\|f\|^2
=
2\pi\sum_mP_m.
$$

Using

$$
\sum_{m\notin S}P_m
=
\sum_mP_m-\sum_{m\in S}P_m
$$

and taking the positive square root gives the result.

PETAL2D reports this quantity as a percentage.

## Optimal power-ranked truncation

Suppose admissible selection units have powers

$$
Q_1\ge Q_2\ge\cdots\ge0.
$$

For a complex field, a unit can be one angular channel. For a real field, the admissible units are $m=0$, conjugate pairs $(m,-m)$, and the even-grid Nyquist singleton when present.

For any fixed number $K$ of units, retaining the first $K$ units minimizes the reconstruction error. The first $K$ that satisfies the requested error tolerance also uses the minimum number of admissible units.

### Proof

The reconstruction-error formula is monotone in retained power. Minimizing the error is equivalent to maximizing

$$
\sum_{j\in S}Q_j
$$

over all admissible sets $S$ containing $K$ units.

If a proposed set contains $Q_a$ and omits a larger $Q_b$, replacing $Q_a$ by $Q_b$ increases retained power. Repeating this exchange produces the set $\{Q_1,\ldots,Q_K\}$. No other $K$-unit set retains more power.

If the first $K$ ranked units fail the error target, every other $K$-unit set retains no more power and must also fail. The first successful ranked value of $K$ is minimal.

## Real-valued fields and conjugate pairs

If $f$ is real-valued,

$$
\rho_{-m}(r)=\rho_m(r)^*.
$$

Consequently,

$$
P_{-m}=P_m.
$$

### Proof

Taking the complex conjugate of the coefficient gives

$$
\rho_m(r)^*
=
\frac{1}{2\pi}
\int_0^{2\pi}
f(r,\theta)e^{im\theta}\,d\theta,
$$

because $f^*=f$. The right-hand side is $\rho_{-m}(r)$.

A conjugate pair reconstructs as

$$
\rho_me^{im\theta}
+
\rho_{-m}e^{-im\theta}
=
2\operatorname{Re}
\left[
\rho_me^{im\theta}
\right],
$$

which is real. This is why PETAL2D keeps conjugate pairs together for real input.

## Rotation covariance

For an active rotation

$$
f_\phi(r,\theta)=f(r,\theta-\phi),
$$

the angular coefficients transform as

$$
\rho_m^{(\phi)}(r)
=
e^{-im\phi}\rho_m(r).
$$

Therefore

$$
P_m^{(\phi)}=P_m.
$$

### Proof

Insert the rotated field into the coefficient definition:

$$
\rho_m^{(\phi)}(r)
=
\frac{1}{2\pi}
\int_0^{2\pi}
f(r,\theta-\phi)e^{-im\theta}\,d\theta.
$$

Set $u=\theta-\phi$. Periodicity restores the integration interval to $[0,2\pi)$ and

$$
e^{-im\theta}
=
e^{-im\phi}e^{-imu}.
$$

Hence

$$
\rho_m^{(\phi)}(r)
=
e^{-im\phi}\rho_m(r).
$$

Taking the modulus squared removes the phase factor.

## Rotational selection rules

Suppose a field transforms under a rotation by $2\pi/n$ as

$$
f(r,\theta+2\pi/n)
=
e^{i2\pi\ell/n}f(r,\theta),
$$

where $\ell$ is an integer modulo $n$. Then any nonzero angular coefficient satisfies

$$
m\equiv\ell\pmod n.
$$

The invariant scalar case corresponds to $\ell=0$.

### Proof

Insert the Fourier expansion:

$$
\sum_m
\rho_m(r)e^{im\theta}e^{i2\pi m/n}
=
e^{i2\pi\ell/n}
\sum_m\rho_m(r)e^{im\theta}.
$$

Uniqueness of Fourier coefficients gives

$$
\rho_m(r)
\left[
e^{i2\pi m/n}
-
e^{i2\pi\ell/n}
\right]
=
0.
$$

For $\rho_m(r)\ne0$,

$$
e^{i2\pi(m-\ell)/n}=1,
$$

so $m-\ell$ is divisible by $n$.

## Discrete angular aliasing

For

$$
\theta_j=\frac{2\pi j}{N_\theta},
\qquad
j=0,\ldots,N_\theta-1,
$$

harmonics whose indices differ by an integer multiple of $N_\theta$ are identical on the sampled angular grid.

### Proof

For any integer $k$,

$$
e^{i(m+kN_\theta)\theta_j}
=
e^{im\theta_j}
e^{i2\pi kj}
=
e^{im\theta_j}.
$$

The samples cannot distinguish $m$ from $m+kN_\theta$. This is the origin of angular aliasing and the Nyquist limit.

## Angular-momentum interpretation for a wavefunction

If the analyzed field is a normalized complex wavefunction,

$$
\psi(r,\theta)
=
\sum_m\rho_m(r)e^{im\theta},
$$

the fixed-$m$ subspace is the eigenspace of the planar operator

$$
L_z=-i\hbar\frac{\partial}{\partial\theta}
$$

with eigenvalue $m\hbar$.

The squared norm of the projection onto this subspace is

$$
\|\Pi_m\psi\|^2
=
2\pi P_m.
$$

Since normalization gives

$$
1=2\pi\sum_nP_n,
$$

the angular-momentum probability is

$$
\Pr(L_z=m\hbar)
=
\frac{P_m}{\sum_nP_n}.
$$

This interpretation applies to the wavefunction itself. A density such as $|\psi|^2$ does not retain the phase required for the same angular-momentum interpretation.

## Translation of the power centroid

Define

$$
\mathbf r_c[f]
=
\frac{
\int \mathbf r\,|f(\mathbf r)|^2\,d^2r
}{
\int |f(\mathbf r)|^2\,d^2r
}.
$$

For

$$
g(\mathbf r)=f(\mathbf r-\mathbf a),
$$

and an integration domain translated with the field,

$$
\mathbf r_c[g]
=
\mathbf r_c[f]+\mathbf a.
$$

### Proof

Set $\mathbf u=\mathbf r-\mathbf a$. The denominator is unchanged. The numerator becomes

$$
\int
(\mathbf u+\mathbf a)
|f(\mathbf u)|^2\,d^2u
=
\int
\mathbf u|f(\mathbf u)|^2\,d^2u
+
\mathbf a
\int
|f(\mathbf u)|^2\,d^2u.
$$

Dividing by the common denominator gives the result.
