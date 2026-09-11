from __future__ import annotations

from collections.abc import Mapping
from enum import Enum

import numpy as np
from matplotlib import pyplot as plt
from numpy.fft import fft, fftfreq, ifft
from scipy.integrate import trapezoid
from scipy.ndimage import map_coordinates


class _FunctionType(Enum):
    CALLABLE = 0
    DISCRETE = 1


class PolarDecomposition(Mapping):
    r"""Polar angular-harmonic decomposition of a localized 2D scalar field.

    ``PolarDecomposition`` converts a scalar field from Cartesian coordinates
    to a polar representation about a chosen origin and expands the angular
    dependence at every radius as

    .. math::

        f(r,\theta) = \sum_m \rho_m(r)e^{im\theta}.

    PETAL2D defines the radial power of angular harmonic ``m`` as

    .. math::

        P_m = \int_0^{r_{\max}} |\rho_m(r)|^2 r\,dr.

    The angular FFT is computed first for the complete discrete spectrum.
    Adaptive truncation then ranks admissible angular channels by power and
    retains the smallest set that satisfies ``recon_err_tol``. For real input,
    nonzero conjugate harmonics ``(+m, -m)`` are selected together so the
    automatic reconstruction remains real-valued.

    Parameters
    ----------
    f_xy : callable or ndarray of shape (Nx, Ny)
        Input field. A callable must accept array-valued ``x`` and ``y`` and
        return either a scalar or an array with the supplied mesh shape.
        Callable inputs are evaluated directly on the final polar grid after
        centering. Array inputs are interpolated from the Cartesian grid.
    x, y : array_like
        One-dimensional, finite, strictly increasing, uniformly spaced
        Cartesian coordinates. For array input, ``f_xy.shape`` must equal
        ``(len(x), len(y))``.
    Nr : int or None, optional
        Number of radial samples. If ``None``, use
        ``min(len(x), len(y))``. The radial grid includes both zero and
        ``rmax``.
    Ntheta : int, default=512
        Number of uniformly spaced angular samples on ``[0, 2*pi)``.
    rmax : float or None, optional
        Maximum analyzed polar radius. If ``None``, use the largest disk
        centered on ``origin`` that is fully contained in the Cartesian
        rectangle. The resulting value is available as ``safe_rmax``.
    m_abs_max : int or None, optional
        Optional bound on admissible ``|m|`` during adaptive selection. It may
        not exceed ``m_nyquist``. This does not change the complete spectrum
        stored in ``m_all`` and ``rho_all``.
    recon_err_tol : float, default=1.0
        Target relative polar-area-weighted L2 reconstruction error, in
        percent. The selector uses a Parseval power identity; ``target_reached``
        reports whether the admissible spectrum can meet the target.
    radial_power_tail_fraction : float, default=1e-6
        Maximum fraction of a retained mode's integrated radial power allowed
        outside ``radial_power_support_radius[m]``. Must satisfy
        ``0 <= radial_power_tail_fraction < 1``.
    radial_relative_amplitude_threshold : float, default=1e-3
        Fraction of ``max_r |rho_m(r)|`` used to define
        ``radial_amplitude_support_radius[m]``. The outermost threshold
        crossing is used so significant outer radial lobes are retained.
        Must satisfy ``0 <= value < 1``.
    origin : {"centroid"} or array_like of length 2, default="centroid"
        Expansion origin. ``"centroid"`` uses the Cartesian L2 power centroid
        of ``|f|^2``. A finite length-two input is interpreted as the explicit
        Cartesian coordinate ``(x0, y0)``.
    normalize : bool, default=False
        If ``True``, divide the Cartesian field by the square root of its
        Cartesian L2 power before polar analysis.
    interp_method : {"linear", "cubic"}, default="cubic"
        Cartesian-to-polar interpolation used only for sampled array input.

    Attributes
    ----------
    origin : tuple of float
        Expansion origin ``(x0, y0)`` actually used.
    centroid : tuple of float
        Numerically evaluated Cartesian L2 power centroid, regardless of
        whether an explicit ``origin`` was requested.
    r : ndarray of shape (Nr,)
        Radial sample coordinates.
    theta : ndarray of shape (Ntheta,)
        Angular sample coordinates in radians.
    safe_rmax : float
        Radius of the largest centered disk fully contained in the Cartesian
        rectangle.
    m_nyquist : int
        Largest representable absolute angular-frequency class implied by
        ``Ntheta``.
    nyquist_mode : int or None
        For even ``Ntheta``, the single NumPy FFT Nyquist index
        ``-Ntheta/2``; otherwise ``None``.
    f_polar : ndarray of shape (Nr, Ntheta)
        Input field evaluated or interpolated on the polar grid.
    f_recon : ndarray of shape (Nr, Ntheta)
        Adaptive retained-mode reconstruction.
    cartesian_power : float
        Cartesian physical L2 power after optional normalization.
    polar_power : float
        Physical L2 power in the analyzed polar domain.
    domain_consistency : float
        Ratio ``polar_power / cartesian_power``. It quantifies captured field
        weight, not pointwise interpolation accuracy.
    is_real_input : bool
        Whether the analyzed input is real-valued to PETAL2D's numerical
        tolerance.
    m_all : ndarray of shape (Ntheta,)
        Complete discrete FFT harmonic indices in NumPy FFT order.
    rho_all : ndarray of shape (Nr, Ntheta)
        Complete radial harmonic coefficients.
    powers_all : ndarray of shape (Ntheta,)
        Complete radial mode powers ``P_m``.
    power_fracs_all : ndarray of shape (Ntheta,)
        ``powers_all`` normalized by their sum.
    m_sorted : list of int
        Adaptively retained harmonics, ordered by selected power units.
    rho : dict[int, ndarray]
        Retained radial harmonic coefficients keyed by ``m``.
    power_fracs : dict[int, float]
        Retained mode power fractions keyed by ``m``.
    selected_pairs : list of tuple[int, ...]
        Adaptive selection units actually retained. Real fields use conjugate
        ``(+m, -m)`` pairs except for ``m=0`` and the even-grid Nyquist
        singleton. Complex fields use singleton channels.
    selected_pair_powers : ndarray
        Power of each retained adaptive selection unit.
    selected_pair_power_fracs : ndarray
        Fractional power of each retained adaptive selection unit.
    retained_power_fraction : float
        Fraction of the complete discrete polar spectral power retained.
    recon_error : float
        Parseval-predicted adaptive reconstruction error in percent.
    recon_error_measured : float
        Directly evaluated polar weighted L2 reconstruction error in percent.
    target_reached : bool
        Whether the admissible angular spectrum satisfies ``recon_err_tol``.
    parseval_relative_error : float
        Relative mismatch between direct polar power and the sum of discrete
        harmonic powers; useful as a numerical consistency diagnostic.
    radial_power_support_radius : dict[int, float]
        Per-retained-mode smallest radius satisfying the integrated radial
        power-tail criterion.
    radial_amplitude_support_radius : dict[int, float]
        Per-retained-mode outermost radius satisfying the relative-amplitude
        criterion.
    cutoff_radius : dict[int, float]
        Per-retained-mode maximum of the power- and amplitude-support radii.
        This is a diagnostic only and never truncates reconstruction.

    Notes
    -----
    PETAL2D uses the Fourier convention

    .. math::

        \rho_m(r) = \frac{1}{2\pi}\int_0^{2\pi}
        f(r,\theta)e^{-im\theta}\,d\theta.

    Therefore

    .. math::

        \|f\|_{L^2(D)}^2 = 2\pi\sum_m P_m.

    The factor ``2*pi`` cancels from all power fractions and reconstruction
    error ratios. The package documentation contains full proofs of Parseval
    error certification, real-field conjugacy, rotation covariance, centering,
    discrete aliasing, and rotational selection rules.

    Examples
    --------
    Analyze a real ``p_x``-like orbital. The automatic selector retains the
    conjugate pair ``(+1, -1)``::

        >>> import numpy as np
        >>> from petal2d import PolarDecomposition
        >>> x = np.linspace(-5.0, 5.0, 101)
        >>> y = np.linspace(-5.0, 5.0, 101)
        >>> def px(x, y):
        ...     return x * np.exp(-0.5 * (x*x + y*y))
        >>> dec = PolarDecomposition(px, x, y, Ntheta=128)
        >>> dec.selected_pairs
        [(1, -1)]
        >>> dec.reconstruction_error() < 1e-10
        True

    See Also
    --------
    numpy.fft.fft : Discrete Fourier transform used for the angular spectrum.
    scipy.ndimage.map_coordinates : Interpolation used for sampled input.
    """

    _REAL_INPUT_TOL = 1e-12

    def __init__(
        self,
        f_xy,
        x: np.ndarray,
        y: np.ndarray,
        Nr: int | None = None,
        Ntheta: int = 512,
        rmax: float | None = None,
        m_abs_max: int | None = None,
        recon_err_tol: float = 1.0,
        radial_power_tail_fraction: float = 1e-6,
        radial_relative_amplitude_threshold: float = 1e-3,
        origin="centroid",
        normalize: bool = False,
        interp_method: str = "cubic",
    ):
        self._classification = self._classify_function_type(f_xy)
        self._f_xy = f_xy

        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self._Nx = self.x.size
        self._Ny = self.y.size

        self._Nr = Nr
        self._Ntheta = Ntheta
        self._rmax = rmax
        self._m_abs_max = m_abs_max
        self._recon_err_tol = recon_err_tol
        self.radial_power_tail_fraction = float(radial_power_tail_fraction)
        self.radial_relative_amplitude_threshold = float(
            radial_relative_amplitude_threshold
        )
        self._origin_spec = origin
        self._normalize = bool(normalize)
        self._interp_method = interp_method

        self._validate_inputs()
        self._prepare_field()
        self._decompose()

    # ------------------------------------------------------------------
    # Mapping interface: retained spectrum
    # ------------------------------------------------------------------
    def __getitem__(self, key):
        if isinstance(key, tuple):
            if len(key) == 0:
                raise KeyError("Empty index.")
            m = int(key[0])
            radial_index = key[1:]
        else:
            m = int(key)
            radial_index = ()

        if m not in self.rho:
            raise KeyError(
                f"Harmonic m={m} is not retained. Available modes: "
                f"{list(self.rho.keys())}"
            )
        rho_m = self.rho[m]
        return rho_m[radial_index] if radial_index else rho_m

    def __iter__(self):
        return iter(self.rho)

    def __len__(self):
        return len(self.rho)

    def __contains__(self, m):
        try:
            return int(m) in self.rho
        except (TypeError, ValueError):
            return False

    # ------------------------------------------------------------------
    # Validation and field preparation
    # ------------------------------------------------------------------
    @staticmethod
    def _classify_function_type(input_func) -> _FunctionType:
        if callable(input_func):
            return _FunctionType.CALLABLE
        try:
            arr = np.asarray(input_func)
        except Exception as exc:  # pragma: no cover
            raise TypeError(f"Could not convert f_xy to an array: {exc}") from exc
        if arr.size == 0:
            raise TypeError("Input field is empty.")
        if arr.ndim != 2:
            raise TypeError(f"Discrete input field is {arr.ndim}D; expected 2D.")
        return _FunctionType.DISCRETE

    @staticmethod
    def _validate_coordinate_grid(values: np.ndarray, name: str) -> None:
        if values.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional.")
        if values.size < 2:
            raise ValueError(f"{name} must contain at least two points.")
        if not np.all(np.isfinite(values)):
            raise ValueError(f"{name} contains non-finite values.")
        delta = np.diff(values)
        if np.any(delta <= 0):
            raise ValueError(f"{name} must be strictly increasing.")
        if not np.allclose(delta, delta[0], rtol=1e-9, atol=1e-12):
            raise ValueError(
                f"{name} must be uniformly spaced for map_coordinates interpolation."
            )

    @staticmethod
    def _validate_positive_int(value, name: str, minimum: int = 1) -> int:
        if not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer.")
        value = int(value)
        if value < minimum:
            raise ValueError(f"{name} must be >= {minimum}.")
        return value

    @staticmethod
    def _coerce_callable_output(values, shape: tuple[int, int], name: str) -> np.ndarray:
        arr = np.asarray(values)
        if arr.ndim == 0:
            arr = np.full(shape, arr, dtype=arr.dtype)
        if arr.shape != shape:
            raise ValueError(
                f"Callable {name} returned shape {arr.shape}; expected {shape}."
            )
        return arr

    @classmethod
    def _is_effectively_real(cls, values: np.ndarray) -> bool:
        arr = np.asarray(values)
        if not np.iscomplexobj(arr):
            return True
        norm = float(np.linalg.norm(arr.ravel()))
        if norm == 0.0:
            return True
        imag_norm = float(np.linalg.norm(np.imag(arr).ravel()))
        return imag_norm / norm <= cls._REAL_INPUT_TOL

    def _validate_inputs(self) -> None:
        self._validate_coordinate_grid(self.x, "x")
        self._validate_coordinate_grid(self.y, "y")

        if self._classification is _FunctionType.DISCRETE:
            arr = np.asarray(self._f_xy)
            expected = (self._Nx, self._Ny)
            if arr.shape != expected:
                raise ValueError(
                    f"Discrete f_xy has shape {arr.shape}; expected {expected} "
                    "to match x and y."
                )

        if self._Nr is None:
            self._Nr = min(self._Nx, self._Ny)
        self._Nr = self._validate_positive_int(self._Nr, "Nr", minimum=2)
        self._Ntheta = self._validate_positive_int(self._Ntheta, "Ntheta", minimum=2)

        self.m_nyquist = self._Ntheta // 2
        self.nyquist_mode = -self.m_nyquist if self._Ntheta % 2 == 0 else None

        if self._m_abs_max is not None:
            self._m_abs_max = self._validate_positive_int(
                self._m_abs_max, "m_abs_max", minimum=0
            )
            if self._m_abs_max > self.m_nyquist:
                raise ValueError(
                    f"m_abs_max={self._m_abs_max} exceeds the angular Nyquist "
                    f"limit m_nyquist={self.m_nyquist} for Ntheta={self._Ntheta}."
                )

        if not np.isfinite(self._recon_err_tol) or self._recon_err_tol < 0:
            raise ValueError("recon_err_tol must be a finite non-negative percentage.")
        if (
            not np.isfinite(self.radial_power_tail_fraction)
            or self.radial_power_tail_fraction < 0
            or self.radial_power_tail_fraction >= 1
        ):
            raise ValueError(
                "radial_power_tail_fraction must satisfy 0 <= value < 1."
            )
        if (
            not np.isfinite(self.radial_relative_amplitude_threshold)
            or self.radial_relative_amplitude_threshold < 0
            or self.radial_relative_amplitude_threshold >= 1
        ):
            raise ValueError(
                "radial_relative_amplitude_threshold must satisfy "
                "0 <= value < 1."
            )

        if self._rmax is not None:
            if not np.isfinite(self._rmax) or self._rmax <= 0:
                raise ValueError("rmax must be a finite positive number.")
            self._rmax = float(self._rmax)

        if self._interp_method not in {"linear", "cubic"}:
            raise ValueError("interp_method must be either 'linear' or 'cubic'.")

        if isinstance(self._origin_spec, str):
            if self._origin_spec != "centroid":
                raise ValueError("origin string must be 'centroid'.")
        else:
            try:
                origin = np.asarray(self._origin_spec, dtype=float)
            except Exception as exc:
                raise TypeError(
                    "origin must be 'centroid' or a finite length-two coordinate."
                ) from exc
            if origin.shape != (2,) or not np.all(np.isfinite(origin)):
                raise ValueError(
                    "origin must be 'centroid' or a finite length-two coordinate."
                )
            self._origin_spec = (float(origin[0]), float(origin[1]))

    def _integrate_cartesian(self, values: np.ndarray):
        return trapezoid(trapezoid(values, self.y, axis=1), self.x, axis=0)

    def _polar_mean_power(self, values: np.ndarray) -> float:
        radial = np.mean(np.abs(values) ** 2, axis=1) * self.r
        return float(np.real(trapezoid(radial, self.r)))

    def _prepare_field(self) -> None:
        """Prepare Cartesian and polar representations without retaining work arrays.

        The Cartesian and polar coordinate meshes used during preparation can be
        large. They are intentionally local temporaries rather than object
        attributes: after ``f_polar`` has been constructed, PETAL2D does not need
        those meshes for decomposition, reconstruction, plotting, or diagnostics.
        This keeps the persistent and peak memory footprint tied more closely to
        the actual field/spectrum arrays rather than coordinate bookkeeping.
        """
        if self._classification is _FunctionType.CALLABLE:
            # Preserve the original callable contract (full Cartesian mesh input)
            # while keeping these meshes local to preparation.
            X_cart, Y_cart = np.meshgrid(self.x, self.y, indexing="ij")
            sampled = self._coerce_callable_output(
                self._f_xy(X_cart, Y_cart),
                (self._Nx, self._Ny),
                "f_xy(x, y)",
            )
        else:
            sampled = np.asarray(self._f_xy)

        if not np.all(np.isfinite(sampled)):
            raise ValueError("f_xy contains non-finite values.")

        density = np.abs(sampled) ** 2
        raw_cartesian_power = float(np.real(self._integrate_cartesian(density)))
        if not np.isfinite(raw_cartesian_power) or raw_cartesian_power <= 0:
            raise ValueError("f_xy must have non-zero finite L2 power.")

        self._normalization_factor = (
            np.sqrt(raw_cartesian_power) if self._normalize else 1.0
        )
        normalized_sampled = (
            sampled
            if self._normalization_factor == 1.0
            else sampled / self._normalization_factor
        )
        density = np.abs(normalized_sampled) ** 2
        self.cartesian_power = float(np.real(self._integrate_cartesian(density)))

        # The Cartesian first moments can be evaluated by broadcasting 1D
        # coordinate vectors; full persistent X/Y meshes are unnecessary.
        self._x_centroid = float(
            np.real(self._integrate_cartesian(self.x[:, None] * density))
            / self.cartesian_power
        )
        self._y_centroid = float(
            np.real(self._integrate_cartesian(self.y[None, :] * density))
            / self.cartesian_power
        )
        self.centroid = (self._x_centroid, self._y_centroid)

        if self._origin_spec == "centroid":
            self.origin = self.centroid
            self.origin_mode = "centroid"
        else:
            self.origin = self._origin_spec
            self.origin_mode = "explicit"

        x0, y0 = self.origin
        if not (self.x[0] <= x0 <= self.x[-1] and self.y[0] <= y0 <= self.y[-1]):
            raise ValueError("origin must lie inside the supplied Cartesian domain.")

        self.x_c = self.x - x0
        self.y_c = self.y - y0
        self.safe_rmax = float(
            min(
                x0 - self.x[0],
                self.x[-1] - x0,
                y0 - self.y[0],
                self.y[-1] - y0,
            )
        )

        if self._rmax is None:
            if self.safe_rmax <= 0:
                raise ValueError(
                    "The selected origin leaves no non-zero complete circle inside "
                    "the Cartesian domain; provide an interior origin."
                )
            self._rmax = self.safe_rmax

        self.r = np.linspace(0.0, self._rmax, self._Nr)
        self.theta = np.linspace(0.0, 2.0 * np.pi, self._Ntheta, endpoint=False)
        cos_theta = np.cos(self.theta)
        sin_theta = np.sin(self.theta)
        radial_column = self.r[:, None]

        cartesian_is_real = self._is_effectively_real(normalized_sampled)

        if self._classification is _FunctionType.CALLABLE:
            # Only two polar coordinate arrays are required, and both are local
            # temporaries released when this method returns.
            X_phys = x0 + radial_column * cos_theta[None, :]
            Y_phys = y0 + radial_column * sin_theta[None, :]
            polar_values = self._coerce_callable_output(
                self._f_xy(X_phys, Y_phys),
                (self._Nr, self._Ntheta),
                "f_xy on the polar grid",
            )
            self.f_polar = (
                polar_values
                if self._normalization_factor == 1.0
                else polar_values / self._normalization_factor
            )
        else:
            # map_coordinates accepts coordinates with shape
            # (ndim, output_shape...). Build that array directly so we do not
            # retain or duplicate R/TH/Xp/Yp meshes.
            dx = self.x[1] - self.x[0]
            dy = self.y[1] - self.y[0]
            coords = np.empty((2, self._Nr, self._Ntheta), dtype=float)
            coords[0] = (
                x0 - self.x[0] + radial_column * cos_theta[None, :]
            ) / dx
            coords[1] = (
                y0 - self.y[0] + radial_column * sin_theta[None, :]
            ) / dy
            interp_order = 1 if self._interp_method == "linear" else 3
            self.f_polar = map_coordinates(
                normalized_sampled,
                coords,
                order=interp_order,
                mode="constant",
                cval=0.0,
            )

        if not np.all(np.isfinite(self.f_polar)):
            raise ValueError("Polar sampling produced non-finite values.")

        self.is_real_input = cartesian_is_real and self._is_effectively_real(
            self.f_polar
        )
        if self.is_real_input:
            self.f_polar = np.real(self.f_polar)

        # Physical polar power includes the angular integral. The internal FFT
        # power convention uses the angular mean, so the factor 2*pi is applied
        # only for comparison with Cartesian-domain power.
        self.polar_power = 2.0 * np.pi * self._polar_mean_power(self.f_polar)
        self.domain_consistency = self.polar_power / self.cartesian_power

    # ------------------------------------------------------------------
    # Spectral decomposition and adaptive selection
    # ------------------------------------------------------------------
    def _weighted_relative_error_percent(
        self, f_orig: np.ndarray, f_recon: np.ndarray
    ) -> float:
        num = self._polar_mean_power(f_recon - f_orig)
        den = self._polar_mean_power(f_orig)
        if den <= 0:
            return np.nan
        return 100.0 * np.sqrt(max(num, 0.0) / den)

    def _compute_cutoff_radii(
        self, rho_m_r: np.ndarray, powers: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute the two support criteria and the final cutoff radius for retained modes.

        Three radii are returned for every retained harmonic ``m``:

        ``radial_power_support_radius``
            Smallest radius enclosing at least
            ``1 - radial_power_tail_fraction`` of the integrated radial mode
            power ``int |rho_m(r)|^2 r dr``.

        ``radial_amplitude_support_radius``
            Outermost radius where ``|rho_m(r)|`` is at least
            ``radial_relative_amplitude_threshold`` times that mode's peak
            radial amplitude.

        ``cutoff_radius``
            The larger of the power-support and amplitude-support radii. Beyond
            this radius, the retained mode satisfies both configured tail
            criteria. It is the dashed radius drawn in the radial-harmonic plot.
            It is diagnostic only and never truncates the reconstruction.
        """
        amplitudes = np.abs(rho_m_r)
        radial_power_density = amplitudes**2 * self.r[:, None]
        dr = self.r[1] - self.r[0]
        cumulative_power = np.zeros_like(radial_power_density, dtype=float)
        cumulative_power[1:] = dr * np.cumsum(
            0.5 * (radial_power_density[:-1] + radial_power_density[1:]), axis=0
        )

        power_targets = (1.0 - self.radial_power_tail_fraction) * powers
        power_support = np.full(powers.size, np.nan, dtype=float)
        amplitude_support = np.full(powers.size, np.nan, dtype=float)

        for j, (total_power, target_power) in enumerate(zip(powers, power_targets)):
            if total_power <= 0:
                continue

            # Power-support radius: linearly interpolate the first cumulative
            # crossing of the requested enclosed-power target.
            idx = int(np.searchsorted(cumulative_power[:, j], target_power))
            if idx <= 0:
                power_support[j] = self.r[0]
            elif idx >= self.r.size:
                power_support[j] = self.r[-1]
            else:
                r0, r1 = self.r[idx - 1], self.r[idx]
                c0, c1 = cumulative_power[idx - 1, j], cumulative_power[idx, j]
                power_support[j] = (
                    r0
                    if c1 == c0
                    else r0 + (target_power - c0) * (r1 - r0) / (c1 - c0)
                )

            # Amplitude-support radius: use the outermost threshold crossing,
            # so oscillatory radial profiles retain every significant outer lobe.
            amp = amplitudes[:, j]
            peak = float(np.max(amp))
            if peak <= 0:
                continue
            threshold = self.radial_relative_amplitude_threshold * peak
            above = np.flatnonzero(amp >= threshold)
            if above.size == 0:
                amplitude_support[j] = self.r[0]
                continue
            k = int(above[-1])
            if k >= self.r.size - 1:
                amplitude_support[j] = self.r[-1]
            else:
                a0, a1 = float(amp[k]), float(amp[k + 1])
                r0, r1 = float(self.r[k]), float(self.r[k + 1])
                if a1 == a0:
                    amplitude_support[j] = r0
                else:
                    crossing = r0 + (threshold - a0) * (r1 - r0) / (a1 - a0)
                    amplitude_support[j] = float(np.clip(crossing, r0, r1))

        combined_support = np.fmax(power_support, amplitude_support)
        return power_support, amplitude_support, combined_support

    def _build_selection_pairs(self, allowed_mask: np.ndarray):
        """Return mode/index pairs used by the adaptive selector."""
        mode_to_index = {
            int(m): int(i)
            for i, m in enumerate(self.m_all)
            if allowed_mask[i]
        }

        pairs: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
        if not self.is_real_input:
            for i in np.flatnonzero(allowed_mask):
                pairs.append(((int(self.m_all[i]),), (int(i),)))
            return pairs

        if 0 in mode_to_index:
            pairs.append(((0,), (mode_to_index[0],)))

        for m_abs in range(1, self.m_nyquist + 1):
            if self._Ntheta % 2 == 0 and m_abs == self.m_nyquist:
                m = -m_abs
                if m in mode_to_index:
                    pairs.append(((m,), (mode_to_index[m],)))
                continue

            plus = mode_to_index.get(m_abs)
            minus = mode_to_index.get(-m_abs)
            if plus is not None and minus is not None:
                pairs.append(((m_abs, -m_abs), (plus, minus)))
            elif plus is not None:  # defensive; symmetric m_abs_max normally prevents this
                pairs.append(((m_abs,), (plus,)))
            elif minus is not None:
                pairs.append(((-m_abs,), (minus,)))
        return pairs

    def _decompose(self) -> None:
        self.m_all = np.rint(fftfreq(self._Ntheta, d=1.0 / self._Ntheta)).astype(int)
        self.rho_all = fft(self.f_polar, axis=1) / self._Ntheta

        radial_density_all = np.abs(self.rho_all) ** 2 * self.r[:, None]
        self.powers_all = trapezoid(radial_density_all, self.r, axis=0)
        self.total_power = float(np.sum(self.powers_all))
        if self.total_power <= 0:
            raise ValueError("The polar field has zero L2 power over the selected radius.")
        self.power_fracs_all = self.powers_all / self.total_power

        polar_mean_power = self._polar_mean_power(self.f_polar)
        self.parseval_relative_error = abs(polar_mean_power - self.total_power) / polar_mean_power

        if self._m_abs_max is None:
            allowed_mask = np.ones(self._Ntheta, dtype=bool)
        else:
            allowed_mask = np.abs(self.m_all) <= self._m_abs_max

        allowed_indices = np.flatnonzero(allowed_mask)
        self.m_vals = self.m_all[allowed_indices].copy()
        self.allowed_power_fraction = float(
            np.sum(self.powers_all[allowed_indices]) / self.total_power
        )

        raw_pairs = self._build_selection_pairs(allowed_mask)
        pair_powers = np.array(
            [float(np.sum(self.powers_all[list(indices)])) for _, indices in raw_pairs]
        )
        order = np.argsort(-pair_powers, kind="stable")
        ranked_pairs = [raw_pairs[i] for i in order]
        ranked_powers = pair_powers[order]

        cumulative_power = np.cumsum(ranked_powers)
        residual_fraction = np.maximum(0.0, 1.0 - cumulative_power / self.total_power)
        pair_errors = 100.0 * np.sqrt(residual_fraction)

        reached = np.flatnonzero(pair_errors <= self._recon_err_tol)
        if reached.size:
            n_pairs_keep = int(reached[0] + 1)
            self.target_reached = True
        else:
            n_pairs_keep = len(ranked_pairs)
            self.target_reached = bool(
                pair_errors[-1]
                <= self._recon_err_tol + 100.0 * np.finfo(float).eps
            )

        selected_pair_records = ranked_pairs[:n_pairs_keep]
        self.selected_pairs = [modes for modes, _ in selected_pair_records]
        self.selected_pair_powers = ranked_powers[:n_pairs_keep].copy()
        self.selected_pair_power_fracs = self.selected_pair_powers / self.total_power
        selected_indices = np.array(
            [i for _, indices in selected_pair_records for i in indices], dtype=int
        )
        self.m_sorted = [
            int(m) for modes, _ in selected_pair_records for m in modes
        ]
        self.Nm = len(self.m_sorted)

        self.powers = self.powers_all[selected_indices].copy()
        self._power_fracs = self.powers / self.total_power
        self._rho_m_r = self.rho_all[:, selected_indices].copy()
        self.retained_power_fraction = float(np.sum(self.powers) / self.total_power)
        self.recon_error = 100.0 * np.sqrt(
            max(0.0, 1.0 - self.retained_power_fraction)
        )

        self.f_recon = self.reconstruct()
        self.recon_error_measured = self._weighted_relative_error_percent(
            self.f_polar, self.f_recon
        )
        if self.is_real_input:
            denom = max(float(np.linalg.norm(np.real(self.f_recon))), np.finfo(float).tiny)
            self.real_reconstruction_residual = float(
                np.linalg.norm(np.imag(self.f_recon)) / denom
            )
        else:
            self.real_reconstruction_residual = np.nan

        (
            self._radial_power_support_radius,
            self._radial_amplitude_support_radius,
            self._cutoff_radius,
        ) = self._compute_cutoff_radii(self._rho_m_r, self.powers)
        self.rho = {
            int(m): self._rho_m_r[:, i] for i, m in enumerate(self.m_sorted)
        }
        self.power_fracs = {
            int(m): float(self._power_fracs[i]) for i, m in enumerate(self.m_sorted)
        }
        self.radial_power_support_radius = {
            int(m): float(self._radial_power_support_radius[i])
            for i, m in enumerate(self.m_sorted)
        }
        self.radial_amplitude_support_radius = {
            int(m): float(self._radial_amplitude_support_radius[i])
            for i, m in enumerate(self.m_sorted)
        }
        self.cutoff_radius = {
            int(m): float(self._cutoff_radius[i])
            for i, m in enumerate(self.m_sorted)
        }

    # ------------------------------------------------------------------
    # Public reconstruction API
    # ------------------------------------------------------------------
    def reconstruct(self, modes=None) -> np.ndarray:
        """Reconstruct the polar field from selected angular harmonics.

        Parameters
        ----------
        modes : None, {"all"}, or iterable of int, optional
            Harmonics to include. ``None`` uses the adaptively retained
            spectrum. ``"all"`` uses the complete discrete FFT spectrum. An
            iterable uses exactly the requested harmonics. Explicit mode lists
            are never silently conjugate-paired.

        Returns
        -------
        ndarray of shape (Nr, Ntheta)
            Field reconstructed on the stored polar grid. For real input,
            adaptive, complete, or explicitly conjugate-closed requests return
            a real array. An intentionally asymmetric explicit subset can
            return a complex array.

        Notes
        -----
        Reconstruction uses the full analyzed radial grid. ``cutoff_radius``
        is diagnostic and never crops the result.
        """
        if modes is None:
            requested = list(self.m_sorted)
        elif isinstance(modes, str):
            if modes != "all":
                raise ValueError("modes string must be 'all'.")
            requested = [int(m) for m in self.m_all]
        else:
            try:
                requested = [int(m) for m in modes]
            except TypeError as exc:
                raise TypeError("modes must be None, 'all', or an iterable of integers.") from exc

        if len(set(requested)) != len(requested):
            raise ValueError("modes contains duplicate harmonic indices.")

        mode_to_index = {int(m): i for i, m in enumerate(self.m_all)}
        missing = [m for m in requested if m not in mode_to_index]
        if missing:
            raise ValueError(
                f"Requested mode(s) {missing} are outside the discrete angular spectrum."
            )

        spectrum = np.zeros_like(self.rho_all)
        for m in requested:
            i = mode_to_index[m]
            spectrum[:, i] = self.rho_all[:, i]
        reconstructed = ifft(spectrum * self._Ntheta, axis=1)

        if self.is_real_input and self._modes_preserve_real_valuedness(requested):
            return np.real(reconstructed)
        return reconstructed

    def _modes_preserve_real_valuedness(self, modes) -> bool:
        """Return whether a requested mode set is closed under conjugation."""
        requested = set(int(m) for m in modes)
        for m in requested:
            if m == 0:
                continue
            if self.nyquist_mode is not None and m == self.nyquist_mode:
                continue
            if -m not in requested:
                return False
        return True

    def reconstruction_error(self, modes=None) -> float:
        """Return polar-area-weighted relative L2 reconstruction error.

        Parameters
        ----------
        modes : None, {"all"}, or iterable of int, optional
            Uses the same mode-selection semantics as :meth:`reconstruct`.

        Returns
        -------
        float
            Relative weighted L2 error in **percent**.

        Notes
        -----
        The norm uses the polar area measure ``r dr dtheta``. For the adaptive
        subset this directly measured value should agree with ``recon_error``,
        which is predicted from omitted spectral power by Parseval's identity.
        """
        return self._weighted_relative_error_percent(
            self.f_polar, self.reconstruct(modes=modes)
        )

    def spectrum_table(self, retained_only: bool = True) -> list[dict]:
        """Return spectrum rows as dependency-free Python dictionaries.

        Parameters
        ----------
        retained_only : bool, default=True
            If ``True``, return only adaptively retained harmonics. If ``False``,
            return every discrete FFT harmonic.

        Returns
        -------
        list of dict
            Rows sorted by decreasing individual harmonic power. Every row has
            keys ``m``, ``power``, ``power_fraction``,
            ``radial_power_support_radius``,
            ``radial_amplitude_support_radius``, ``cutoff_radius``, and
            ``retained``. Radius entries are ``None`` for non-retained modes.

        Notes
        -----
        ``power_fraction`` refers to an individual harmonic, whereas
        ``selected_pair_power_fracs`` refers to adaptive selection units when
        real-field conjugate pairing is active. No pandas dependency is used.
        """
        retained = set(self.m_sorted)
        rows = []
        for m, power, fraction in zip(
            self.m_all, self.powers_all, self.power_fracs_all
        ):
            m_int = int(m)
            is_retained = m_int in retained
            if retained_only and not is_retained:
                continue
            rows.append(
                {
                    "m": m_int,
                    "power": float(power),
                    "power_fraction": float(fraction),
                    "radial_power_support_radius": (
                        float(self.radial_power_support_radius[m_int])
                        if is_retained else None
                    ),
                    "radial_amplitude_support_radius": (
                        float(self.radial_amplitude_support_radius[m_int])
                        if is_retained else None
                    ),
                    "cutoff_radius": (
                        float(self.cutoff_radius[m_int])
                        if is_retained else None
                    ),
                    "retained": is_retained,
                }
            )
        rows.sort(key=lambda row: (-row["power"], row["m"]))
        return rows

    def print_spectrum(self, retained_only: bool = True, file=None) -> None:
        """Print a compact human-readable table of the angular spectrum.

        Parameters
        ----------
        retained_only : bool, default=True
            Print only retained modes by default. Set ``False`` to print the
            complete discrete FFT spectrum, including roundoff/interpolation
            leakage.
        file : file-like object or None, optional
            Destination passed to :func:`print`. ``None`` writes to standard
            output.

        Returns
        -------
        None
        """
        rows = self.spectrum_table(retained_only=retained_only)
        header = (
            f"{'m':>6}  {'power':>14}  {'fraction':>12}  "
            f"{'r_power_support':>16}  {'r_amplitude_support':>20}  "
            f"{'r_cutoff':>12}  {'retained':>8}"
        )
        print(header, file=file)
        print("-" * len(header), file=file)
        for row in rows:
            def fmt(value):
                return "-" if value is None else f"{value:.6g}"

            print(
                f"{row['m']:6d}  {row['power']:14.6e}  "
                f"{row['power_fraction']:12.6e}  "
                f"{fmt(row['radial_power_support_radius']):>16}  "
                f"{fmt(row['radial_amplitude_support_radius']):>20}  "
                f"{fmt(row['cutoff_radius']):>12}  "
                f"{str(row['retained']):>8}",
                file=file,
            )

    # ------------------------------------------------------------------
    # Plotting helpers
    # ------------------------------------------------------------------
    # ------------------------------------------------------------------
    # Plotting helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _plot_rc_params():
        """Return local Matplotlib settings used by PETAL2D plotting helpers.

        The settings are applied with ``matplotlib.pyplot.rc_context`` so the
        plotting helpers do not modify a user's global Matplotlib style.
        """
        return {
            "font.size": 10,
            "axes.titlesize": 10.5,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "legend.fontsize": 8.5,
            "legend.title_fontsize": 8.5,
            "lines.linewidth": 1.6,
            "figure.titlesize": 12,
            "savefig.bbox": "tight",
        }

    @staticmethod
    def _clean_axes(ax, *, grid: bool = False):
        """Apply restrained, publication-style axis formatting in-place."""
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if grid:
            ax.grid(True, which="major", linewidth=0.6, alpha=0.20)
            ax.set_axisbelow(True)

    def plot_harmonics_with_hist(self, title: str = ""):
        """Plot retained radial profiles and ranked angular power fractions.

        Parameters
        ----------
        title : str, optional
            Figure-level title.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Created figure.
        axes : ndarray of matplotlib.axes.Axes
            Two axes: radial profiles on the left and retained power-ranked
            harmonics on the right.

        Notes
        -----
        The left panel shows ``|rho_m(r)|`` for each retained harmonic. A
        color-matched dashed line marks ``cutoff_radius[m]``. The cutoff radius
        is diagnostic only and never truncates reconstruction.

        The right panel is categorical rather than a linear ``m`` axis. Only
        retained harmonics are shown, ordered from greatest to least fractional
        power, so sparse spectra with large angular indices do not create large empty ranges.
        """
        mode_to_index = {int(m): i for i, m in enumerate(self.m_sorted)}
        ranked_modes = sorted(
            (int(m) for m in self.m_sorted),
            key=lambda m: (-float(self.power_fracs[m]), mode_to_index[m]),
        )

        with plt.rc_context(self._plot_rc_params()):
            fig, axs = plt.subplots(
                1,
                2,
                figsize=(9.0, 3.6),
                gridspec_kw={"width_ratios": (1.35, 1.0)},
                layout="constrained",
            )
            ax1, ax2 = axs

            color_by_mode = {}
            for m in ranked_modes:
                i = mode_to_index[m]
                abs_rho = np.abs(self._rho_m_r[:, i])
                line, = ax1.plot(
                    self.r,
                    abs_rho,
                    label=rf"$m={m}$",
                )
                color_by_mode[m] = line.get_color()
                if np.isfinite(self._cutoff_radius[i]):
                    ax1.axvline(
                        self._cutoff_radius[i],
                        color=line.get_color(),
                        linestyle=(0, (4, 2.5)),
                        linewidth=1.15,
                        alpha=0.85,
                    )

            ax1.set_xlabel(r"$r$")
            ax1.set_ylabel(r"$|\rho_m(r)|$")
            ax1.set_ylim(bottom=0.0)
            ax1.set_title("Radial profiles", pad=8)
            self._clean_axes(ax1, grid=True)

            mode_legend = ax1.legend(
                title="Retained harmonics",
                frameon=True,
                framealpha=0.92,
                facecolor="white",
                edgecolor="none",
                loc="upper right",
                handlelength=2.3,
                borderaxespad=0.4,
            )
            ax1.add_artist(mode_legend)

            cutoff_handle = plt.Line2D(
                [0],
                [0],
                color="0.25",
                linestyle=(0, (4, 2.5)),
                linewidth=1.15,
                label="Cutoff radius",
            )
            ax1.legend(
                handles=[cutoff_handle],
                frameon=False,
                loc="lower right",
                handlelength=2.3,
                borderaxespad=0.4,
            )

            ranked_fractions = np.array(
                [self.power_fracs[m] for m in ranked_modes], dtype=float
            )
            positions = np.arange(len(ranked_modes), dtype=float)
            bars = ax2.bar(
                positions,
                ranked_fractions,
                width=0.68,
                linewidth=0.7,
                edgecolor="0.25",
                color=[color_by_mode[m] for m in ranked_modes],
            )

            max_fraction = float(np.max(ranked_fractions))
            ymax = min(1.08, max(0.12, 1.18 * max_fraction))
            ax2.set_ylim(0.0, ymax)
            ax2.set_xlim(-0.65, max(0.65, len(ranked_modes) - 0.35))
            ax2.set_xticks(positions)
            ax2.set_xticklabels([rf"${m}$" for m in ranked_modes])
            ax2.set_xlabel(r"Angular harmonic $m$")
            ax2.set_ylabel("Power fraction")
            ax2.set_title("Retained harmonics (power ranked)", pad=8)
            self._clean_axes(ax2, grid=True)

            if len(ranked_modes) <= 12:
                label_offset = 0.025 * ymax
                for bar, frac in zip(bars, ranked_fractions):
                    ax2.text(
                        bar.get_x() + bar.get_width() / 2.0,
                        min(float(frac) + label_offset, 0.985 * ymax),
                        f"{100.0 * float(frac):.1f}%",
                        ha="center",
                        va="bottom",
                        fontsize=8,
                    )

            if title:
                fig.suptitle(title)
            return fig, axs

    def plot_original_vs_reconstructions(
        self,
        title: str = "Reconstruction comparison",
        display_radius: float | None = None,
    ):
        """Compare the analyzed polar field with the adaptive reconstruction.

        Parameters
        ----------
        title : str, default="Reconstruction comparison"
            Figure-level title.
        display_radius : float or None, optional
            Half-width of the displayed Cartesian window in the same length
            units as ``x`` and ``y``. If ``None``, use 115% of the largest
            retained ``cutoff_radius``, capped at the analyzed radial-domain
            radius. This changes only the viewport.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Created figure.
        axes : ndarray of matplotlib.axes.Axes
            Two axes containing the analyzed field and adaptive reconstruction.

        Notes
        -----
        Real sign-changing fields are displayed with sign preserved and a
        symmetric zero-centered color scale. Real nonnegative fields use a
        sequential map. Complex fields are displayed as magnitudes. Both panels
        always share the same color normalization.
        """
        original = np.asarray(self.f_polar)
        reconstruction = np.asarray(self.f_recon)

        if self.is_real_input:
            original_plot = np.real(original)
            reconstruction_plot = np.real(reconstruction)
            combined_min = float(min(np.min(original_plot), np.min(reconstruction_plot)))
            combined_max = float(max(np.max(original_plot), np.max(reconstruction_plot)))
            scale = max(abs(combined_min), abs(combined_max), np.finfo(float).tiny)
            sign_changing = combined_min < -1e-12 * scale and combined_max > 1e-12 * scale

            if sign_changing:
                cmap = "RdBu_r"
                vmin, vmax = -scale, scale
            else:
                cmap = "magma"
                vmin = min(0.0, combined_min)
                vmax = combined_max
            colorbar_label = r"$f(x,y)$"
        else:
            original_plot = np.abs(original)
            reconstruction_plot = np.abs(reconstruction)
            cmap = "magma"
            vmin = 0.0
            vmax = float(max(np.max(original_plot), np.max(reconstruction_plot)))
            colorbar_label = r"$|f(x,y)|$"

        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin + 1.0

        labels = (
            "Original",
            f"Reconstruction\n$L^2$ error = {self.recon_error_measured:.3g}%",
        )
        fields = (original_plot, reconstruction_plot)

        dr = self.r[1] - self.r[0]
        dtheta = self.theta[1] - self.theta[0]
        r_edges = np.concatenate(
            ([max(0.0, self.r[0] - dr / 2.0)],
             (self.r[:-1] + self.r[1:]) / 2.0,
             [self.r[-1] + dr / 2.0])
        )
        theta_edges = np.concatenate(
            ([self.theta[0] - dtheta / 2.0],
             (self.theta[:-1] + self.theta[1:]) / 2.0,
             [self.theta[-1] + dtheta / 2.0])
        )
        R_edges, TH_edges = np.meshgrid(r_edges, theta_edges, indexing="ij")
        X_edges = R_edges * np.cos(TH_edges)
        Y_edges = R_edges * np.sin(TH_edges)

        analyzed_radius = float(self.r[-1])
        if display_radius is None:
            finite_support = [
                float(value)
                for value in self.cutoff_radius.values()
                if np.isfinite(value)
            ]
            if finite_support:
                view_radius = min(analyzed_radius, 1.15 * max(finite_support))
            else:
                view_radius = analyzed_radius
        else:
            view_radius = float(display_radius)
            if not np.isfinite(view_radius) or view_radius <= 0.0:
                raise ValueError("display_radius must be a finite positive number or None")
            view_radius = min(view_radius, analyzed_radius)

        with plt.rc_context(self._plot_rc_params()):
            fig, axs = plt.subplots(
                1,
                2,
                figsize=(8.2, 3.75),
                sharex=True,
                sharey=True,
                layout="constrained",
            )
            pcm = None
            for ax, label, field in zip(axs, labels, fields):
                pcm = ax.pcolormesh(
                    X_edges,
                    Y_edges,
                    field,
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    shading="auto",
                    rasterized=True,
                )
                ax.set_aspect("equal")
                ax.set_xlim(-view_radius, view_radius)
                ax.set_ylim(-view_radius, view_radius)
                ax.set_title(label, pad=8)
                ax.set_xlabel(r"$x-x_0$")
                ax.tick_params(direction="out", length=3.5)
                for spine in ax.spines.values():
                    spine.set_linewidth(0.8)

            axs[0].set_ylabel(r"$y-y_0$")
            axs[1].tick_params(labelleft=False)

            fig.colorbar(
                pcm,
                ax=axs,
                shrink=0.88,
                pad=0.025,
                label=colorbar_label,
            )
            if title:
                fig.suptitle(title)
            return fig, axs
