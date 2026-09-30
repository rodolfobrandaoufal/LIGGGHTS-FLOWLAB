#!/usr/bin/env python3
"""
Parse every input deck under examples/ (and tests/, benchmarks/) for
'pair_style gran ...' and 'fix ... wall/gran ...' contact-model selections and
check them against:
  (1) the local Make.sh whitelist src/style_contact_model.whitelist (untracked, git-ignored),
  (2) the whitelist CMake's WRITE_WHITELIST generates with default ENABLE_MODEL_* options
      (audit/logs/build_release.cmake_generated_style_contact_model.h).
Any selection absent from a list would have run through the (now removed)
fallback path and now hard-errors.  Output: audit/logs/contact_deck_whitelist.txt
"""
import os, re, glob, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SRC = os.path.join(ROOT, "src")

def ids(prefix):
    m = {}
    for f in glob.glob(os.path.join(SRC, prefix + "_model_*.h")):
        for line in open(f, errors="ignore"):
            mm = re.match(r"^%s_MODEL\((\w+),\s*([\w/]+),\s*(\d+)\)" % prefix.upper(), line)
            if mm: m[mm.group(2)] = mm.group(1)
    return m
NM, TM, CM, RM, SM = ids("normal"), ids("tangential"), ids("cohesion"), ids("rolling"), ids("surface")
NM["off"] = "NORMAL_OFF"; TM["off"] = "TANGENTIAL_OFF"; CM["off"] = "COHESION_OFF"; RM["off"] = "ROLLING_OFF"

def load(path):
    s = set()
    for line in open(path):
        mm = re.match(r"GRAN_MODEL\((.*)\)", line.strip())
        if mm: s.add(tuple(x.strip() for x in mm.group(1).split(",")))
    return s
WL_LOCAL = load(os.path.join(SRC, "style_contact_model.whitelist"))
WL_CMAKE = load(os.path.join(ROOT, "audit/logs/build_release.cmake_generated_style_contact_model.h"))

def parse(tokens):
    sel = {"model": "NORMAL_OFF", "tangential": "TANGENTIAL_OFF", "cohesion": "COHESION_OFF",
           "rolling_friction": "ROLLING_OFF", "surface": "SURFACE_DEFAULT"}
    maps = {"model": NM, "tangential": TM, "cohesion": CM, "rolling_friction": RM, "surface": SM}
    i = 0
    for key in ["model", "tangential", "cohesion", "rolling_friction", "surface"]:
        if i + 1 < len(tokens) and tokens[i] == key:
            sel[key] = maps[key].get(tokens[i + 1], "UNKNOWN(" + tokens[i + 1] + ")")
            i += 2
    return (sel["model"], sel["tangential"], sel["cohesion"], sel["rolling_friction"], sel["surface"])

rows = []
decks = [f for pat in ("examples/**/in.*", "tests/**/in.*", "benchmarks/**/in.*")
         for f in glob.glob(os.path.join(ROOT, pat), recursive=True)]
for f in sorted(set(decks)):
    try: lines = open(f, errors="ignore").read().replace("&\n", " ").splitlines()
    except IsADirectoryError: continue
    for ln in lines:
        ln = ln.split("#")[0].strip()
        t = ln.split()
        if len(t) >= 3 and t[0] == "pair_style" and t[1] == "gran":
            combo = parse(t[2:])
        elif len(t) >= 4 and t[0] == "fix" and t[3] == "wall/gran":
            combo = parse(t[4:])
        else:
            continue
        if "NORMAL_OFF" == combo[0] or any(c.startswith("UNKNOWN") for c in combo):
            continue  # variables / custom / SPH walls: not resolvable statically
        rows.append((os.path.relpath(f, ROOT), combo, combo in WL_LOCAL, combo in WL_CMAKE))

rows = sorted(set(rows))
out = ["# deck | combination | in local Make.sh whitelist (124) | in CMake-default whitelist (%d)" % len(WL_CMAKE)]
for r in rows:
    out.append("%s | %s | %s | %s" % (r[0], ",".join(r[1]), "yes" if r[2] else "NO", "yes" if r[3] else "NO"))
nl = sum(1 for r in rows if not r[2]); nc = sum(1 for r in rows if not r[3])
out.append("")
out.append("selections parsed: %d ; missing from local whitelist: %d ; missing from CMake-default whitelist: %d" % (len(rows), nl, nc))
txt = "\n".join(out); print(txt)
open(os.path.join(ROOT, "audit/logs/contact_deck_whitelist.txt"), "w").write(txt + "\n")
