"""Line/branch coverage of src/ from the --coverage build after the CTest run
(gcovr unavailable: its lxml wheel did not load). Runs gcov -n -b per object
and aggregates per source file; headers instantiated in several units take the
best unit (approximation: no union of line sets)."""
import os, re, subprocess, glob, json, collections
W = "/tmp/claude-1000/-media-storage-LIGGGHTS-PUBLIC-v6/0010f84b-565a-43d6-8523-c6e56b72c1bf/scratchpad/audit2_work"
SRC = W + "/clone/src/"
OUT = "/media/storage/LIGGGHTS-PUBLIC-v6/audit2/data/coverage"
files = {}
for gcda in glob.glob(W + "/builds/cov/**/*.gcda", recursive=True):
    obj_dir = os.path.dirname(gcda)
    src = SRC + os.path.basename(gcda)[:-len(".cpp.gcda")] + ".cpp"
    if not os.path.exists(src): continue
    p = subprocess.run(["gcov", "-n", "-b", "-o", gcda[:-len(".gcda")] + ".o", src], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=obj_dir)
    cur = None
    for line in p.stdout.splitlines():
        m = re.match(r"File '(.+)'", line)
        if m: cur = os.path.realpath(m.group(1)); continue
        if cur and cur.startswith(os.path.realpath(SRC)):
            name = cur[len(os.path.realpath(SRC))+1:]
            m = re.match(r"Lines executed:([\d.]+)% of (\d+)", line)
            if m:
                pct, n = float(m.group(1)), int(m.group(2)); d = files.setdefault(name, {})
                if pct*n >= d.get("lcov", -1)*d.get("lines", 0) or "lines" not in d: d["lcov"], d["lines"] = pct, n
            m = re.match(r"Taken at least once:([\d.]+)% of (\d+)", line)
            if m:
                pct, n = float(m.group(1)), int(m.group(2)); d = files.setdefault(name, {})
                if "branches" not in d or pct*n >= d["bcov"]*d["branches"]: d["bcov"], d["branches"] = pct, n
tl = sum(d.get("lines", 0) for d in files.values()); tc = sum(d.get("lines", 0)*d.get("lcov", 0)/100 for d in files.values())
tb = sum(d.get("branches", 0) for d in files.values()); tbc = sum(d.get("branches", 0)*d.get("bcov", 0)/100 for d in files.values())
summary = {"files": len(files), "lines": tl, "lines_covered": round(tc), "line_pct": 100*tc/tl, "branches": tb, "branch_pct": 100*tbc/tb if tb else None}
# branch-changed files
changed = subprocess.run(["git", "-C", "/media/storage/LIGGGHTS-PUBLIC-v6", "diff", "--name-only", "3d5c00f2", "HEAD", "--", "src"], stdout=subprocess.PIPE, text=True).stdout.split()
ch = {f[4:] for f in changed}
cl = sum(d.get("lines", 0) for f, d in files.items() if f in ch); cc = sum(d.get("lines", 0)*d.get("lcov", 0)/100 for f, d in files.items() if f in ch)
summary["changed_files_line_pct"] = 100*cc/cl if cl else None; summary["changed_files_lines"] = cl
json.dump({"summary": summary, "files": files}, open(OUT + "/coverage.json", "w"), indent=1)
with open(OUT + "/coverage_by_file.txt", "w") as f:
    for name, d in sorted(files.items(), key=lambda x: -x[1].get("lines", 0)):
        f.write(f"{name:50s} lines {d.get('lines',0):6d} {d.get('lcov',0):6.1f}%  branches {d.get('branches',0):6d} {d.get('bcov',0):6.1f}%{'  [changed]' if name in ch else ''}\n")
print(json.dumps(summary, indent=1))
