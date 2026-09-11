# PETAL2D performance benchmarks

Performance benchmarks are separate from numerical validation. Compute scripts
write **JSON only**; plotting scripts consume those JSON files.

## Runtime

```bash
python3 benchmarks/benchmark_runtime.py --quick
python3 benchmarks/benchmark_runtime.py
python3 benchmarks/plot_runtime.py
```

Runtime measurements use warm-up runs, `time.perf_counter`, and report the
median and interquartile range. Callable and sampled inputs are benchmarked
separately. The benchmark-only subclass records two coarse stages:

- `prepare_field`: input evaluation/sampling, Cartesian preparation, and
  Cartesian-to-polar interpolation or direct polar evaluation;
- `decompose`: FFT, harmonic powers, adaptive selection, reconstruction, and
  cutoff diagnostics.

No benchmark timing is part of PETAL2D's public API.

## Peak memory

Install the optional benchmark dependency first:

```bash
python3 -m pip install -e ".[benchmark]"
```

Then run:

```bash
python3 benchmarks/benchmark_memory.py --quick
python3 benchmarks/benchmark_memory.py
python3 benchmarks/plot_memory.py
```

Each repetition runs in a fresh Python process. The user inputs are constructed
first; for sampled input, benchmark-only Cartesian X/Y construction meshes are
deleted and garbage-collected. Peak resident-set size (RSS) is then measured
relative to that **input-resident baseline**, so the reported incremental peak
counts PETAL2D work rather than the user's pre-existing field array or benchmark
construction temporaries. This is a memory-footprint measurement; PETAL2D does
**not** infer or claim hardware memory bandwidth from RSS measurements.

Results are written to `benchmarks/results/`; figures go to
`benchmarks/results/figures/`.
