# Applications and interpretation patterns

PETAL2D is a general angular-analysis tool, but its interpretation depends on what the input field represents. The same mathematical spectrum can carry different physical meaning for a wavefunction, a density, an order parameter, or an image.

## Complex wavefunctions and angular momentum

For a normalized complex wavefunction on the analyzed disk, with angular momentum defined about the chosen PETAL2D origin, the power fraction

$$
w_m=\frac{P_m}{\sum_nP_n}
$$

is exactly the probability weight of the $L_z=m\hbar$ angular sector within that disk, with radial degrees of freedom left unresolved. This follows because $e^{im\theta}$ diagonalizes $L_z$ and the angular sectors are orthogonal.

This is useful for planar model wavefunctions, defect states, excitonic relative-coordinate states, and localized continuum eigenmodes.

## Atomic-orbital and Wannier character

For a localized orbital centered on an atomic, molecular, or Wannier site, the dominant $m$ channels provide an angular fingerprint. Familiar real orbitals appear as conjugate combinations: $p_x$ contains $m=\pm1$, while $d_{x^2-y^2}$ contains $m=\pm2$.

In a crystal, however, $m$ is usually a basis label, not an exact isolated-atom quantum number. Crystal fields can mix angular channels allowed by the local point-group symmetry. PETAL2D is especially useful for quantifying that mixing rather than forcing one orbital label.

## Local crystal symmetry

A converged spectrum can expose rotational selection rules. An exactly $C_3$-invariant scalar field about the correct center contains only $m=0,\pm3,\pm6,\ldots$. A field transforming under a nontrivial one-dimensional $C_n$ character occupies one congruence class $m\equiv\ell\pmod n$.

This makes the method useful for localized crystal states, high-symmetry regions, and local orbital textures. Weak symmetry-incompatible modes should be convergence-tested before being interpreted as symmetry breaking because Cartesian-to-polar interpolation can generate small leakage.

## Densities and scalar observables

For charge density, probability density, a scalar spin-density component, or another real observable, PETAL2D describes shape anisotropy. It does not recover phase information that is absent from the observable. For example, $|e^{im\theta}|^2=1$, so an angular-momentum eigenstate can have an isotropic probability density.

## Complex order parameters and vortices

A complex scalar order parameter with phase winding $\Delta(r,\theta)\sim g(r)e^{im\theta}$ appears naturally as a dominant complex $m$ channel. PETAL2D can separate mixed windings while retaining each radial envelope. This can be useful for superfluid/superconducting order parameters, scalar wave vortices, and related 2D fields.

## Photonic, phononic, acoustic, and PDE modes

Any localized scalar eigenmode of a two-dimensional wave or differential equation can be analyzed in the same way. The angular spectrum can compare mode families, symmetry breaking, defect-induced mixing, or changes in radial localization across parameters.

## Scientific imaging and morphology

For a localized scalar intensity map, PETAL2D can serve as an interpretable angular descriptor. Unlike a generic image feature vector, the output retains a direct meaning in terms of angular harmonics and radial profiles. This is most appropriate when a physically meaningful center exists.

## When PETAL2D is the wrong tool

Use a different representation when the problem is intrinsically three-dimensional, lacks a meaningful center, is dominated by translation/plane-wave structure, or requires coupled vector/tensor symmetry analysis. A Cartesian FFT, spherical harmonics, or problem-specific symmetry decomposition may then be more natural.
