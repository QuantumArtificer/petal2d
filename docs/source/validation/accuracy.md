# Numerical accuracy

The validation scripts construct fields with known harmonic content and save machine-readable JSON before plotting. This keeps calculations separate from presentation and makes paper figures reproducible.

## What is tested

- Cartesian interpolation leakage for known $m=0,\pm3$ fields
- mode-power errors versus Cartesian resolution
- radial quadrature convergence
- default radial-resolution strategies
- angular Nyquist behavior and exact FFT aliasing
- adaptive truncation against analytically known spectral powers
- Parseval-predicted versus directly measured reconstruction errors.

## Main conclusions from the current validation campaign

The current benchmark set shows approximately second-order radial quadrature convergence, exact recovery of the discrete FFT alias in the tested analytic angular cases, and near-machine-precision agreement between Parseval-predicted and directly measured truncation errors. Cubic Cartesian-to-polar interpolation strongly suppresses spurious harmonic leakage relative to linear interpolation on smooth fields.

The automatic choice `Nr=min(len(x), len(y))` was selected because it provides a substantially better radial-accuracy/cost compromise than the previous half-resolution rule.

```{figure} ../_static/validation/accuracy_default_nr_strategy.png
:alt: Accuracy comparison of radial resolution strategies
:width: 90%

Accuracy/cost motivation for the default radial resolution.
```

```{figure} ../_static/validation/accuracy_radial_convergence.png
:alt: Radial convergence
:width: 90%

Independent radial-convergence study.
```

```{figure} ../_static/validation/truncation_error_certification.png
:alt: Predicted versus measured adaptive truncation error
:width: 75%

Parseval error certification for adaptive truncation.
```

## Reproduce

```bash
python validation/run_accuracy_validation.py
python validation/plot_accuracy_validation.py
python validation/run_angular_resolution_validation.py
python validation/plot_angular_resolution_validation.py
python validation/run_truncation_validation.py
python validation/plot_truncation_validation.py
```

The exact quantitative values in a publication should always be taken from the archived JSON corresponding to the cited software version and benchmark environment.
