#!/usr/bin/env python3
"""
Coverage check for the tracked contact-model whitelist (dispatch fix A4, P0-02/P0-03/P0-05).
Every 'pair_style gran' / 'fix ... wall/gran' selection in examples/, tests/ and benchmarks/ decks
must be listed in src/contact_model_whitelist.txt, so that shipped decks use the fast static path
(any other combination still runs through the runtime fallback, with a warning).
Exit status 0 = all covered. Adapted from audit/scripts/contact/check_deck_whitelist.py.
(LIGGGHTS modernization branch)
"""
import os, re, glob, sys
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
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
WL = load(os.path.join(SRC, "contact_model_whitelist.txt"))
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
# scratch output of test runners is not part of the deck set
decks = [f for f in decks if "/work/" not in f]
for f in sorted(set(decks)):
    try: lines = open(f, errors="ignore").read().replace("&\n", " ").splitlines()
    except IsADirectoryError: continue
    if "whitelist-coverage: skip" in "\n".join(lines): continue
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
        rows.append((os.path.relpath(f, ROOT), combo, combo in WL))


rows = sorted(set(rows))
missing = [r for r in rows if not r[2]]
for r in missing:
    print("NOT WHITELISTED: %s | GRAN_MODEL(%s)" % (r[0], ", ".join(r[1])))
print("whitelist entries: %d ; deck selections: %d ; missing: %d" % (len(WL), len(rows), len(missing)))
print("WHITELIST COVERAGE: " + ("PASS" if not missing else "FAIL"))
sys.exit(1 if missing else 0)
