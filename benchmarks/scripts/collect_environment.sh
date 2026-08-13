#!/usr/bin/env bash
set -eu

out_dir="${1:-benchmarks/environments}"
mkdir -p "$out_dir"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out="$out_dir/environment_$stamp.txt"

{
  echo "timestamp_utc=$stamp"
  echo "repository=$(pwd)"
  git rev-parse --show-toplevel 2>/dev/null || true
  git rev-parse HEAD 2>/dev/null || true
  git status --short 2>/dev/null || true
  uname -a || true
  command -v mpicxx >/dev/null 2>&1 && mpicxx --version || true
  command -v mpirun >/dev/null 2>&1 && mpirun --version || true
  command -v cmake >/dev/null 2>&1 && cmake --version || true
  command -v make >/dev/null 2>&1 && make --version | sed -n '1,2p' || true
  command -v lscpu >/dev/null 2>&1 && lscpu || true
  command -v free >/dev/null 2>&1 && free -h || true
} > "$out"

echo "$out"
