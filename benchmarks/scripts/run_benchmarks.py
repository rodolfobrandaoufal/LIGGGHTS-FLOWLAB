#!/usr/bin/env python3
"""Run benchmark manifest cases and preserve raw logs plus parsed CSV output."""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from parse_liggghts_log import parse_log, write_rows  # noqa: E402


def repo_root():
    return SCRIPT_DIR.parents[1]


def build_command(executable, input_file, ranks, mpirun, variables):
    cmd = []
    if ranks > 1:
        cmd.extend([mpirun, "-np", str(ranks)])
    cmd.extend([executable, "-in", input_file])
    for name, value in variables.items():
        cmd.extend(["-var", str(name), str(value)])
    return cmd


def run_case(root, executable, mpirun, case, repeat, ranks, raw_dir):
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = f"{case['id']}_r{ranks}_rep{repeat}_{timestamp}"
    log_path = raw_dir / f"{run_id}.log"
    input_path = root / case["input"]
    cmd = build_command(executable, str(input_path), ranks, mpirun, case.get("variables", {}))

    with log_path.open("w") as log:
        log.write("# command: " + " ".join(cmd) + "\n")
        log.write("# cwd: " + str(root) + "\n")
        log.write("# run_id: " + run_id + "\n")
        log.flush()
        completed = subprocess.run(cmd, cwd=str(root), stdout=log, stderr=subprocess.STDOUT)

    expected_exit_code = int(case.get("expected_exit_code", 0))
    status = "pass" if completed.returncode == expected_exit_code else "fail"
    expected_text = case.get("expected_stderr_contains")
    if expected_text and expected_text not in log_path.read_text(errors="replace"):
        status = "fail"

    row = parse_log(log_path)
    row.update({
        "case": case["id"],
        "run_id": run_id,
        "input": case["input"],
        "ranks": str(ranks),
        "exit_code": str(completed.returncode),
        "status": status,
    })
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="benchmarks/manifests/baseline_smoke.json")
    parser.add_argument("--executable", required=True, help="Path to lmp executable")
    parser.add_argument("--mpirun", default=os.environ.get("MPIRUN", "mpirun"))
    parser.add_argument("--case", action="append", help="Run only named case; may be repeated")
    parser.add_argument("--repeat", type=int, help="Override manifest repeat count")
    parser.add_argument("--ranks", type=int, action="append", help="Override ranks; may be repeated")
    args = parser.parse_args()

    root = repo_root()
    manifest_path = root / args.manifest
    manifest = json.loads(manifest_path.read_text())
    repeats = args.repeat if args.repeat is not None else int(manifest.get("repeats", 1))
    raw_dir = root / "benchmarks" / "raw_logs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    csv_path = root / manifest.get("output", {}).get("parsed_csv", "benchmarks/parsed_csv/benchmark_results.csv")

    selected = set(args.case or [])
    rows = []
    for case in manifest["cases"]:
        if selected and case["id"] not in selected:
            continue
        rank_list = args.ranks or case.get("ranks") or manifest.get("default_ranks", [1])
        for ranks in rank_list:
            for repeat in range(1, repeats + 1):
                rows.append(run_case(root, args.executable, args.mpirun, case, repeat, int(ranks), raw_dir))

    write_rows(rows, csv_path, append=False)
    failures = [row for row in rows if row["status"] != "pass"]
    print(f"wrote {len(rows)} rows to {csv_path}")
    if failures:
        for row in failures:
            print(f"FAIL {row['case']} {row['run_id']} exit={row['exit_code']} log={row['log_file']}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
