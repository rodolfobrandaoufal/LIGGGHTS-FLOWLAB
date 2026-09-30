#!/usr/bin/env python3
"""Summaries for Phase-3 campaigns: mean ± std (n), min, paired ratio vs reference with 95% CI, Welch p."""
import json, sys, statistics as st
from scipy import stats
L = '/media/storage/LIGGGHTS-PUBLIC-v6/audit/logs/perf/'
def load(c):
    return [json.loads(l) for l in open(L + c + '/results.jsonl')]
def table(c, ref, keys=('loop','pair','neigh','comm','modfy','outpt','other'), extra=(), rows=None):
    rows = rows or load(c)
    labs = list(dict.fromkeys(r['label'] for r in rows))
    ok = lambda r: r.get('loop') is not None  # 3 runs have rc=1: deck edited after their 'run' finished (timings valid)
    out = [f'### {c}', '', '| variant | n | ' + ' | '.join(keys + tuple(extra)) + ' |', '|' + '---|' * (2 + len(keys) + len(extra))]
    for lab in labs:
        rs = [r for r in rows if r['label'] == lab and ok(r)]
        cells = []
        for k in keys + tuple(extra):
            v = [r[k] for r in rs if r.get(k) is not None]
            cells.append('%.4g ± %.2g (min %.4g)' % (st.mean(v), st.stdev(v) if len(v) > 1 else 0, min(v)) if v else '-')
        out.append(f'| {lab} | {len(rs)} | ' + ' | '.join(cells) + ' |')
    out.append('')
    if ref:
        out += ['| variant / ' + ref + ' | metric | paired ratio mean ± std | 95% CI | Welch p | ratio of mins |', '|---|---|---|---|---|---|']
        for lab in labs:
            if lab == ref: continue
            for k in keys[:5]:
                a = {r['rep']: r[k] for r in rows if r['label'] == lab and ok(r) and r.get(k)}
                b = {r['rep']: r[k] for r in rows if r['label'] == ref and ok(r) and r.get(k)}
                common = sorted(set(a) & set(b))
                if len(common) < 2: continue
                q = [a[i] / b[i] for i in common]
                m, s = st.mean(q), st.stdev(q); h = stats.t.ppf(0.975, len(q) - 1) * s / len(q) ** 0.5
                p = stats.ttest_ind(list(a.values()), list(b.values()), equal_var=False).pvalue
                out.append(f'| {lab} | {k} | {m:.4f} ± {s:.4f} (n={len(q)}) | [{m-h:.4f}, {m+h:.4f}] | {p:.3g} | {min(a.values())/min(b.values()):.4f} |')
        out.append('')
    return '\n'.join(out)
if __name__ == '__main__':
    print(table(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
