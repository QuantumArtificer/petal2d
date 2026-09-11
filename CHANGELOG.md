# Changelog

All notable changes to PETAL2D will be documented in this file.

## [0.1.0] - 2026-09-11
- Add a full Sphinx/PyData documentation site with pedagogical getting-started,
  narrative user guide, theory, formal proofs, guided physics examples,
  validation/performance documentation, glossary, limitations, and API reference.
- Expand NumPy-style public API docstrings and add documentation CI.
- Add sampled-field and adaptive-truncation executable examples.
- Add a safe cleanup utility that removes generated local artifacts only when
  they are not Git-tracked.
- Rewrite the README around installation, quick start, applications, diagnostics,
  documentation, validation, citation, and contribution guidance.
- Change the automatic radial-grid default from `Nr=min(Nx,Ny)//2` to
  `Nr=min(Nx,Ny)`, based on the validated accuracy/cost tradeoff.
- Add runtime and peak-memory cost studies for the three radial-resolution
  strategies `Nr=Nxy/2`, `Nr=Nxy`, and `Nr=2*Nxy`.
- Strengthen sampled-input memory benchmarking with 7 fresh-process repeats in
  full mode, larger Cartesian grids through `2048^2`, median/IQR reporting, and
  no forced affine slope fit for non-linear allocator/interpolation behavior.
- Add optional Linux CPU-affinity control and record CPU governor/affinity in
  benchmark metadata.
- Add a clean full-run script that removes stale generated results before
  regenerating all validation and benchmark outputs.
- Reduce PETAL2D memory use by keeping Cartesian/polar coordinate meshes as
  preparation-local temporaries instead of persistent object attributes.
- Correct sampled-input memory benchmarking so the reported incremental peak
  excludes the pre-existing user field and benchmark-only X/Y construction meshes.
- Replace the forced Ntheta*log(Ntheta) runtime reference with an empirical affine
  fit to measured complete-runtime scaling; FFT asymptotics remain a theoretical
  component-level statement only.
- Add an end-to-end validation study comparing Nr=floor(Nxy/2), Nr=Nxy, and
  Nr=2*Nxy before changing the package default.
- Separate dimensionless radial power error from dimensional cutoff-radius error
  in validation figures.
- Present the retained angular-power panel with x-axis label `Angular harmonic m`
  and title `Retained harmonics (power ranked)`; ranking remains categorical so
  sparse/high-|m| spectra never create empty linear harmonic ranges.
- Split paper-oriented numerical validation into accuracy/interpolation, angular
  resolution/Nyquist, and adaptive-truncation studies with JSON-first outputs.
- Replace legacy performance scripts with reproducible runtime and fresh-process
  peak-RSS benchmarks, each with separate JSON plotting scripts.
- Rename the final combined radial diagnostic from `radial_support_radius` to `cutoff_radius`; the two underlying power- and amplitude-support radii remain explicit.

### Added
- Polar angular-harmonic decomposition for analytic and sampled 2D fields.
- Power-centroid or explicit-origin expansion.
- Polar/Cartesian `domain_consistency` diagnostic.
- Angular Nyquist validation and complete-spectrum access.
- Adaptive truncation with conjugate-pair preservation for real-valued fields.
- Public reconstruction, reconstruction-error, and spectrum-table APIs.
- Mathematical test suite, JSON validation workflows, physics-motivated examples,
  and reproducible benchmark utilities.
- Consistent physics-example reporting and plotting; clarified radial-cutoff diagnostics.
- Replaced the ambiguous radial-cutoff API with explicit power-tail and
  relative-amplitude support diagnostics; the plotted cutoff radius is
  the stricter of the two criteria.
- Replaced the singular C3 example with a smooth Cartesian polynomial
  realization of the m=+/-3 harmonic.
