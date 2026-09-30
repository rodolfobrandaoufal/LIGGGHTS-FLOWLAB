#!/usr/bin/env python3
"""Rewrite *_component(i,c) accessor calls back to direct array indexing and
diff against HEAD, to prove the accessor refactor is a pure textual rename.
Usage: normalize_accessors.py src/dump_custom.cpp"""
import re, subprocess, sys, difflib
path = sys.argv[1]
new = open(path).read()
old = subprocess.run(["git","show","HEAD:"+path],capture_output=True,text=True).stdout
for name,arr in [("x","x"),("v","v"),("f","f"),("omega","omega"),("torque","torque")]:
    new = re.sub(r"\b%s_component\(\s*([^,()]+?)\s*,\s*([^()]+?)\s*\)"%name, r"%s[\1][\2]"%arr, new)
# drop accessor definitions block and local array declarations for the comparison
def strip(s):
    out=[]
    for l in s.splitlines():
        if re.match(r"\s*double\s*\*\*\s*(x|v|f|omega|torque)\s*=\s*atom->(x|v|f|omega|torque)\s*;\s*$", l): continue
        out.append(l.rstrip())
    return out
d = list(difflib.unified_diff(strip(old), strip(new), "HEAD", "normalized", n=0, lineterm=""))
print("\n".join(d))
print("DIFF_LINES", sum(1 for l in d if l[:1] in "+-" and not l.startswith(("+++","---"))))
