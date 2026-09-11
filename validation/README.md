# PETAL2D numerical validation

The validation suite is deliberately separate from the fast unit/mathematical
tests. Compute scripts write **JSON only**; plotting scripts consume the saved
JSON and never recompute the numerical study.

## Accuracy and interpolation

```bash
python3 validation/run_accuracy_validation.py --quick
python3 validation/run_accuracy_validation.py
python3 validation/plot_accuracy_validation.py
```

This studies Cartesian-to-polar interpolation leakage, recovery of known
harmonic powers, radial quadrature convergence, and convergence of the two
components that define `cutoff_radius`. It also compares three end-to-end radial
resolution strategies, `Nr=floor(Nxy/2)`, `Nr=Nxy`, and `Nr=2*Nxy`, so the
default radial-grid rule is chosen from measured accuracy rather than assumption.

## Angular resolution and Nyquist behavior

```bash
python3 validation/run_angular_resolution_validation.py --quick
python3 validation/run_angular_resolution_validation.py
python3 validation/plot_angular_resolution_validation.py
```

This verifies discrete FFT aliasing, pure-harmonic recovery below the Nyquist
limit, and real-valued conjugate-pair behavior at the even-grid Nyquist bin.

## Adaptive truncation

```bash
python3 validation/run_truncation_validation.py --quick
python3 validation/run_truncation_validation.py
python3 validation/plot_truncation_validation.py
```

The synthetic fields have independently normalized radial basis functions and
known angular power fractions, so theoretical truncation errors and selected
modes/pairs are known before PETAL2D is run.

Full-result JSON files are written to `validation/results/`. Figures are written
to `validation/results/figures/`.
