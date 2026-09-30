#!/usr/bin/env python3
"""Merge audit/findings/*.csv into audit/01_findings.csv.

Every row is kept. A `duplicate_of` column points a duplicate at its canonical
finding; `source` names the originating fragment. Rows are sorted by severity,
then non-legacy before legacy, then id.
"""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.path.dirname(HERE)
ORDER = ["phase0", "contact", "fixes", "vv", "perf", "sota", "quality", "phaseC"]
SEV = {"P0-correctness": 0, "P1-performance": 1, "P1-scalability": 1,
       "P1-build": 1, "P2-physics": 2, "P2-io": 2, "P3-usability": 3}

# duplicate -> canonical
DUP = {
    "F-01": "C-01", "S-01": "C-01",
    "S-02": "F-02",
    "C-02": "P0-01", "S-10": "P0-01",
    "C-03": "P0-03",
    "C-05": "P0-04",
    "S-14": "F-15", "P0-11": "F-15",
    "S-20": "C-13",
    "S-12": "C-06",
    "S-03": "F-07",
    "S-18": "Q-05", "F-19": "Q-05",
    "P0-09": "F-26",
    # verification rows that measured an earlier code-inspection finding
    "V-01": "C-18", "V-02": "C-23", "V-03": "C-16", "V-06": "C-06",
    "V-07": "C-19", "V-12": "C-14", "V-13": "F-04", "V-14": "F-05",
    "V-15": "F-06", "V-16": "F-22", "V-17": "F-23",
    "PF-01": "F-15", "PF-02": "F-25", "PF-06": "S-07", "PF-11": "S-05",
}

FIELDS = ["id", "severity", "category", "file_line", "summary",
          "failure_scenario", "evidence_label", "evidence_ref",
          "recommendation", "effort", "legacy"]

rows = []
for name in ORDER:
    path = os.path.join(AUDIT, "findings", name + ".csv")
    if not os.path.exists(path):
        continue
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            r = {k: (r.get(k) or "").strip() for k in FIELDS}
            r["source"] = name
            r["duplicate_of"] = DUP.get(r["id"], "")
            rows.append(r)

# A measured duplicate upgrades the canonical row's evidence.
byid = {r["id"]: r for r in rows}
for r in rows:
    canon = byid.get(r["duplicate_of"])
    if canon and r["evidence_label"].lower().startswith("measured") \
            and not canon["evidence_label"].lower().startswith("measured"):
        canon["evidence_label"] = "measured (via %s)" % r["id"]
        canon["evidence_ref"] += "; " + r["evidence_ref"]

rows.sort(key=lambda r: (SEV.get(r["severity"], 9),
                         r["legacy"].lower().startswith("y"),
                         r["id"]))

out = os.path.join(AUDIT, "01_findings.csv")
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=FIELDS + ["duplicate_of", "source"],
                       quoting=csv.QUOTE_MINIMAL)
    w.writeheader()
    w.writerows(rows)

uniq = [r for r in rows if not r["duplicate_of"]]
by = {}
for r in uniq:
    key = (r["severity"], "legacy" if r["legacy"].lower().startswith("y") else "branch")
    by[key] = by.get(key, 0) + 1
print(f"{len(rows)} rows, {len(uniq)} unique -> {out}")
for k in sorted(by, key=lambda k: (SEV.get(k[0], 9), k[1])):
    print(f"  {k[0]:16s} {k[1]:7s} {by[k]}")
