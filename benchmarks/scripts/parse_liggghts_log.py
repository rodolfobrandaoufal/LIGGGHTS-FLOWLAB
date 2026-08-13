#!/usr/bin/env python3
"""Parse LIGGGHTS/LAMMPS-style log files into benchmark CSV rows."""

import argparse
import csv
import re
from pathlib import Path


FIELDNAMES = [
    "case", "run_id", "input", "ranks", "exit_code", "status",
    "loop_time", "procs", "steps", "atoms", "time_per_step",
    "pair_time", "neigh_time", "comm_time", "output_time", "modify_time",
    "nlocal_avg", "nlocal_max", "nlocal_min",
    "nghost_avg", "nghost_max", "nghost_min",
    "total_neighbors", "neighbor_builds", "dangerous_builds",
    "log_file", "error_message",
]


LOOP_RE = re.compile(
    r"Loop time of\s+(?P<loop>[0-9.eE+-]+)\s+on\s+(?P<procs>\d+)\s+procs\s+for\s+"
    r"(?P<steps>\d+)\s+steps\s+with\s+(?P<atoms>\d+)\s+atoms"
)
TIMING_RE = re.compile(
    r"^(?P<label>Pair|Neigh|Comm|Output|Modify)\s+time(?:\s*\(%\))?\s*=\s+(?P<value>[0-9.eE+-]+)",
    re.IGNORECASE,
)
RANGE_RE = re.compile(
    r"^(?P<label>Nlocal|Nghost)\s*:\s+(?P<avg>[0-9.eE+-]+)\s+ave\s+"
    r"(?P<max>[0-9.eE+-]+)\s+max\s+(?P<min>[0-9.eE+-]+)\s+min",
    re.IGNORECASE,
)
TOTAL_NEIGH_RE = re.compile(r"^Total\s+#\s+of\s+neighbors\s*=\s+(?P<value>[0-9.eE+-]+)", re.IGNORECASE)
NEIGH_BUILDS_RE = re.compile(r"^Neighbor\s+list\s+builds\s*=\s+(?P<value>\d+)", re.IGNORECASE)
DANGEROUS_RE = re.compile(r"^Dangerous\s+(?:builds|reneighborings)\s*=\s+(?P<value>\d+)", re.IGNORECASE)
ERROR_RE = re.compile(r"\b(ERROR|Error)\b[: ]+(?P<message>.*)")


def parse_log(path):
    row = {name: "" for name in FIELDNAMES}
    row["log_file"] = str(path)
    text = Path(path).read_text(errors="replace")

    for raw_line in text.splitlines():
        line = raw_line.strip()

        match = LOOP_RE.search(line)
        if match:
            row["loop_time"] = match.group("loop")
            row["procs"] = match.group("procs")
            row["steps"] = match.group("steps")
            row["atoms"] = match.group("atoms")
            try:
                steps = float(match.group("steps"))
                loop_time = float(match.group("loop"))
                if steps > 0.0:
                    row["time_per_step"] = "%.12g" % (loop_time / steps)
            except ValueError:
                pass
            continue

        match = TIMING_RE.search(line)
        if match:
            label = match.group("label").lower()
            row["neigh_time" if label == "neigh" else f"{label}_time"] = match.group("value")
            continue

        match = RANGE_RE.search(line)
        if match:
            label = match.group("label").lower()
            row[f"{label}_avg"] = match.group("avg")
            row[f"{label}_max"] = match.group("max")
            row[f"{label}_min"] = match.group("min")
            continue

        match = TOTAL_NEIGH_RE.search(line)
        if match:
            row["total_neighbors"] = match.group("value")
            continue

        match = NEIGH_BUILDS_RE.search(line)
        if match:
            row["neighbor_builds"] = match.group("value")
            continue

        match = DANGEROUS_RE.search(line)
        if match:
            row["dangerous_builds"] = match.group("value")
            continue

        match = ERROR_RE.search(line)
        if match and not row["error_message"]:
            row["error_message"] = match.group("message").strip()

    return row


def write_rows(rows, output, append=False):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    exists = output.exists() and output.stat().st_size > 0
    with output.open("a" if append else "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        if not append or not exists:
            writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", help="Log files to parse")
    parser.add_argument("--output", "-o", default="benchmarks/parsed_csv/parsed_logs.csv")
    parser.add_argument("--append", action="store_true")
    parser.add_argument("--case", default="")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--input", default="")
    parser.add_argument("--ranks", default="")
    parser.add_argument("--exit-code", default="")
    parser.add_argument("--status", default="")
    args = parser.parse_args()

    rows = []
    for log in args.logs:
        row = parse_log(log)
        row.update({
            "case": args.case,
            "run_id": args.run_id,
            "input": args.input,
            "ranks": args.ranks,
            "exit_code": args.exit_code,
            "status": args.status,
        })
        rows.append(row)

    write_rows(rows, args.output, append=args.append)


if __name__ == "__main__":
    main()
