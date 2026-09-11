# Building the PETAL2D documentation

From the repository root:

```bash
python3 -m pip install -e ".[docs]"
python3 -m sphinx -W --keep-going -b html docs/source docs/_build/html
```

Open `docs/_build/html/index.html` in a browser.

The documentation is written as narrative MyST Markdown plus an autodoc API reference.  The build uses the PyData Sphinx Theme so the local documentation has the same broad information architecture as many Scientific Python projects: getting started, narrative user guide, theory, examples, validation, API reference, and development notes.
