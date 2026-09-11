# Performance and memory

Performance is measured separately from mathematical validation. Runtime results use repeated measurements and report medians and interquartile ranges. Peak-memory measurements use fresh processes and establish baselines after user inputs are resident.

## Computational structure

For a polar grid with $N_rN_\theta$ points, the angular FFT contribution scales as

$$
O(N_rN_\theta\log N_\theta),
$$

while field evaluation, interpolation, radial integration, reconstruction, and several array operations require $O(N_rN_\theta)$ work. Over the measured range, the complete implementation is close to linear even though the FFT term retains its $N\log N$ asymptotic complexity.

For sampled Cartesian input, interpolation can dominate runtime and temporary memory at large Cartesian grids.

```{figure} ../_static/validation/runtime_default_nr_cost.png
:alt: Runtime cost of radial resolution strategies
:width: 80%

Runtime tradeoff among the tested radial-resolution strategies.
```

```{figure} ../_static/validation/memory_polar_grid.png
:alt: Peak memory versus polar grid size
:width: 80%

Incremental peak memory for callable input versus polar-grid size.
```

## Reproduce

```bash
python benchmarks/benchmark_runtime.py
python benchmarks/plot_runtime.py

python -m pip install -e ".[benchmark]"
python benchmarks/benchmark_memory.py
python benchmarks/plot_memory.py
```

Absolute timings are hardware- and power-state-specific. The benchmark JSON records software versions, CPU information, affinity, and relevant threading environment variables so comparisons are interpretable.
