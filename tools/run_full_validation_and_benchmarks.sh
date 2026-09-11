#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$ROOT"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

CPU_ARGS=()
if [ -n "${PETAL2D_BENCH_CPU:-}" ]; then
    CPU_ARGS=(--cpu "$PETAL2D_BENCH_CPU")
fi

echo "============================================================"
echo "PETAL2D full validation + benchmark run"
echo "============================================================"
echo "Repository: $ROOT"
echo "Threads: OMP=$OMP_NUM_THREADS MKL=$MKL_NUM_THREADS OPENBLAS=$OPENBLAS_NUM_THREADS NUMEXPR=$NUMEXPR_NUM_THREADS"
if [ ${#CPU_ARGS[@]} -gt 0 ]; then
    echo "Benchmark CPU affinity: ${PETAL2D_BENCH_CPU}"
else
    echo "Benchmark CPU affinity: not pinned (set PETAL2D_BENCH_CPU=<logical CPU> to pin)"
fi

GOVERNORS="$(for f in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do [ -r "$f" ] && cat "$f"; done 2>/dev/null | sort -u | tr '\n' ' ')"
if [ -n "$GOVERNORS" ]; then
    echo "CPU governor(s): $GOVERNORS"
    if ! echo "$GOVERNORS" | grep -qw performance; then
        echo "WARNING: CPU governor is not 'performance'. Scaling conclusions remain useful, but absolute timing values may be noisier/hardware-state dependent."
    fi
fi

echo
echo "Removing stale generated results..."
rm -rf validation/results benchmarks/results
mkdir -p validation/results benchmarks/results

echo
echo "Running tests..."
python3 -m pytest -q

echo
echo "Running full numerical validation..."
python3 validation/run_accuracy_validation.py
python3 validation/plot_accuracy_validation.py
python3 validation/run_angular_resolution_validation.py
python3 validation/plot_angular_resolution_validation.py
python3 validation/run_truncation_validation.py
python3 validation/plot_truncation_validation.py

echo
echo "Running full runtime benchmarks..."
python3 benchmarks/benchmark_runtime.py "${CPU_ARGS[@]}"
python3 benchmarks/plot_runtime.py

echo
echo "Running full memory benchmarks..."
python3 benchmarks/benchmark_memory.py "${CPU_ARGS[@]}"
python3 benchmarks/plot_memory.py

echo
echo "============================================================"
echo "COMPLETE"
echo "Validation JSON: validation/results/"
echo "Benchmark JSON:  benchmarks/results/"
echo "============================================================"
