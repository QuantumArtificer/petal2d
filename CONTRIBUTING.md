# Contributing to PETAL2D

PETAL2D is a small scientific package, so contributions should keep the public API narrow and the mathematics explicit.

## Development setup

```bash
python -m pip install -e ".[test,docs,benchmark,release]"
```

## Before submitting a change

```bash
python -m pytest -v
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

If a change affects numerical accuracy or performance, rerun the relevant JSON validation/benchmark study and explain what changed.

Public behavior should include tests and documentation. Prefer analytic reference fields for mathematical tests. Do not commit local caches, build products, benchmark-analysis ZIPs, extracted bundles, or `.petal2d_backup_*` directories.

## Scope

Core additions should support PETAL2D's central purpose: polar angular-harmonic analysis of localized 2D scalar fields. Large scope expansions such as 3D spherical harmonics, GUIs, unrelated file readers, or general crystallographic classification should be discussed before implementation.
