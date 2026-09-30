#!/usr/bin/env python3
"""Run standalone GPU_DEM scaffold tests and collect benchmark metadata."""

from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_EXECUTABLES = (
    "gpu_particle_data_test",
    "gpu_integrator_test",
    "gpu_neighbor_builder_test",
    "gpu_contact_pipeline_test",
    "gpu_hooke_collision_test",
    "gpu_hooke_damped_collision_test",
    "gpu_wall_contact_test",
    "gpu_timestep_driver_test",
    "gpu_timestep_profiling_test",
    "gpu_synthetic_benchmark",
)

TIMESTEP_RE = re.compile(r"total=([0-9]+(?:\.[0-9]+)?) ms")
DEVICE_RE = re.compile(r"device\s+(\d+)\s+\(([^)]+)\)")
PARTICLES_RE = re.compile(r"particles=([0-9]+)")
CONTACTS_RE = re.compile(r"contacts=([0-9]+)")
PARTICLE_RATE_RE = re.compile(r"particles_per_second=([0-9.eE+-]+)")
CONTACT_RATE_RE = re.compile(r"contacts_per_second=([0-9.eE+-]+)")


@dataclass
class BenchmarkResult:
    executable: str
    iteration: int
    returncode: int
    passed: bool
    wall_time_ms: float
    cuda_timestep_ms: float | None
    particles: int | None
    contacts: int | None
    particles_per_second: float | None
    contacts_per_second: float | None
    device_id: int | None
    device_name: str | None
    stdout: str
    stderr: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run standalone GPU_DEM benchmark/test executables."
    )
    parser.add_argument(
        "--build-dir",
        default="build/gpu_dem_scaffold_cuda",
        type=Path,
        help="Directory containing GPU_DEM standalone executables.",
    )
    parser.add_argument(
        "--out-dir",
        default="build/gpu_dem_scaffold_cuda/benchmark_results",
        type=Path,
        help="Directory where JSON and CSV result files are written.",
    )
    parser.add_argument(
        "--repeat",
        default=3,
        type=int,
        help="Number of times to run each executable.",
    )
    parser.add_argument(
        "--timeout",
        default=60.0,
        type=float,
        help="Per-executable timeout in seconds.",
    )
    parser.add_argument(
        "--exe",
        action="append",
        default=[],
        help="Executable name to run. May be repeated. Defaults to the full GPU_DEM set.",
    )
    return parser.parse_args()


def extract_cuda_timestep_ms(output: str) -> float | None:
    match = TIMESTEP_RE.search(output)
    if not match:
      return None
    return float(match.group(1))


def extract_device(output: str) -> tuple[int | None, str | None]:
    match = DEVICE_RE.search(output)
    if not match:
        return None, None
    return int(match.group(1)), match.group(2)


def extract_int(pattern: re.Pattern[str], output: str) -> int | None:
    match = pattern.search(output)
    if not match:
        return None
    return int(match.group(1))


def extract_float(pattern: re.Pattern[str], output: str) -> float | None:
    match = pattern.search(output)
    if not match:
        return None
    return float(match.group(1))


def run_one(executable: Path, iteration: int, timeout: float) -> BenchmarkResult:
    start = time.perf_counter()
    completed = subprocess.run(
        [str(executable)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=False,
    )
    wall_time_ms = (time.perf_counter() - start) * 1000.0
    combined_output = completed.stdout + "\n" + completed.stderr
    device_id, device_name = extract_device(combined_output)
    return BenchmarkResult(
        executable=executable.name,
        iteration=iteration,
        returncode=completed.returncode,
        passed=completed.returncode == 0,
        wall_time_ms=wall_time_ms,
        cuda_timestep_ms=extract_cuda_timestep_ms(combined_output),
        particles=extract_int(PARTICLES_RE, combined_output),
        contacts=extract_int(CONTACTS_RE, combined_output),
        particles_per_second=extract_float(PARTICLE_RATE_RE, combined_output),
        contacts_per_second=extract_float(CONTACT_RATE_RE, combined_output),
        device_id=device_id,
        device_name=device_name,
        stdout=completed.stdout.strip(),
        stderr=completed.stderr.strip(),
    )


def write_json(results: Iterable[BenchmarkResult], path: Path) -> None:
    data = [asdict(result) for result in results]
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def write_csv(results: Iterable[BenchmarkResult], path: Path) -> None:
    rows = [asdict(result) for result in results]
    fields = (
        "executable",
        "iteration",
        "returncode",
        "passed",
        "wall_time_ms",
        "cuda_timestep_ms",
        "particles",
        "contacts",
        "particles_per_second",
        "contacts_per_second",
        "device_id",
        "device_name",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def main() -> int:
    args = parse_args()
    if args.repeat <= 0:
        print("--repeat must be positive", file=sys.stderr)
        return 2

    executable_names = tuple(args.exe) if args.exe else DEFAULT_EXECUTABLES
    missing = [name for name in executable_names if not (args.build_dir / name).is_file()]
    if missing:
        print(
            "Missing GPU_DEM executable(s): " + ", ".join(missing),
            file=sys.stderr,
        )
        print(f"Build directory checked: {args.build_dir}", file=sys.stderr)
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    results: list[BenchmarkResult] = []
    for name in executable_names:
        executable = args.build_dir / name
        for iteration in range(1, args.repeat + 1):
            result = run_one(executable, iteration, args.timeout)
            results.append(result)
            status = "PASS" if result.passed else "FAIL"
            cuda_ms = (
                f", cuda_timestep={result.cuda_timestep_ms:.6f} ms"
                if result.cuda_timestep_ms is not None
                else ""
            )
            throughput = (
                f", particles/s={result.particles_per_second:.6e}, "
                f"contacts/s={result.contacts_per_second:.6e}"
                if result.particles_per_second is not None
                and result.contacts_per_second is not None
                else ""
            )
            print(
                f"{status} {name} iter={iteration} "
                f"wall={result.wall_time_ms:.3f} ms{cuda_ms}{throughput}"
            )

    json_path = args.out_dir / "gpu_dem_benchmark_results.json"
    csv_path = args.out_dir / "gpu_dem_benchmark_results.csv"
    write_json(results, json_path)
    write_csv(results, csv_path)

    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    return 0 if all(result.passed for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
