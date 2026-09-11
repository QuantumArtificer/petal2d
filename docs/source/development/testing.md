# Tests and documentation checks

Install development dependencies:

```bash
python -m pip install -e ".[test,docs]"
```

Run the mathematical tests:

```bash
python -m pytest -v
```

Build documentation with warnings treated as errors:

```bash
python -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

Before release, also build and check both distribution artifacts:

```bash
python -m pip install -e ".[release]"
python -m build
python -m twine check dist/*
```
