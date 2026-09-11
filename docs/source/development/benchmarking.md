# Benchmarking policy

PETAL2D treats performance claims as reproducible measurements rather than properties inferred from one profiler run.

## Rules

1. Numerical validation and performance benchmarking remain separate.
2. Compute scripts write JSON first. Plot scripts consume saved JSON.
3. Runtime is reported from repeated runs using robust summaries such as medians and interquartile ranges.
4. Peak memory is measured in fresh processes relative to an explicitly defined baseline.
5. Hardware, software versions, CPU affinity, power governor, and relevant threading variables are recorded.
6. Asymptotic algorithmic complexity is distinguished from an empirical fit over a finite benchmark range.
7. No hardware-memory-bandwidth claim is inferred from RSS measurements.

Use

```bash
./tools/run_full_validation_and_benchmarks.sh
```

to delete stale generated results and reproduce the full suite.
