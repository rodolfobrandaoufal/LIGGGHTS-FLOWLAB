import json, statistics as st, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
L='/media/storage/LIGGGHTS-PUBLIC-v6/audit/logs/perf/'
def rows(c): return [json.loads(l) for l in open(L+c+'/results.jsonl')]
def agg(rs, k): v=[r[k] for r in rs if r.get(k) is not None]; return (st.mean(v), st.stdev(v) if len(v)>1 else 0, min(v))
out=[]
for camp in ('t5_strong_200k','t5_weak_12k'):
    R=rows(camp); nps=sorted(set(r['np'] for r in R)); t1=agg([r for r in R if r['np']==1],'loop')
    out.append(f'\n#### {camp}\n\n| ranks | atoms | loop (mean ± std, min) | speedup / eff. (mean) | eff. (min-based) | Pair % | Comm % | Neigh % | Modify % | Nlocal max/avg | Neighs max/avg | mem/proc MB |\n|---|---|---|---|---|---|---|---|---|---|---|---|')
    for n in nps:
        rs=[r for r in R if r['np']==n and r.get('loop')]; l=agg(rs,'loop')
        if camp.startswith('t5_strong'): sp=t1[0]/l[0]; eff=sp/n; effm=t1[2]/l[2]/n; spt=f'{sp:.2f} / {eff:.2f}'
        else: eff=t1[0]/l[0]; effm=t1[2]/l[2]; spt=f'– / {eff:.2f}'
        fr=lambda k: 100*st.mean([r[k]/r['loop'] for r in rs])
        out.append(f"| {n} | {rs[0]['atoms']} | {l[0]:.2f} ± {l[1]:.2f} ({l[2]:.2f}) | {spt} | {effm:.2f} | {fr('pair'):.1f} | {fr('comm'):.1f} | {fr('neigh'):.1f} | {fr('modfy'):.1f} | {st.mean([r['nlocal_max']/r['nlocal_avg'] for r in rs]):.2f} | {st.mean([r['neighs_max']/r['neighs_avg'] for r in rs]):.2f} | {rs[0].get('mem_mb')} |")
R=rows('t5_imbalance_200k')
out.append('\n#### t5_imbalance_200k\n\n| ranks | variant | loop (mean ± std, min) | Pair % | Comm % | Nlocal max/avg (min) | Neighs max/avg | mem/proc MB |\n|---|---|---|---|---|---|---|---|')
for n in sorted(set(r['np'] for r in R)):
    for lab in dict.fromkeys(r['label'] for r in R):
        rs=[r for r in R if r['np']==n and r['label']==lab and r.get('loop')]
        if not rs: continue
        l=agg(rs,'loop'); fr=lambda k: 100*st.mean([r[k]/r['loop'] for r in rs])
        out.append(f"| {n} | {lab} | {l[0]:.2f} ± {l[1]:.2f} ({l[2]:.2f}) | {fr('pair'):.1f} | {fr('comm'):.1f} | {rs[0]['nlocal_max']/rs[0]['nlocal_avg']:.2f} ({rs[0]['nlocal_min']:.0f}) | {rs[0]['neighs_max']/rs[0]['neighs_avg']:.2f} | {rs[0].get('mem_mb')} |")
print('\n'.join(out))
# plot: parallel efficiency
fig,ax=plt.subplots(figsize=(7,4),dpi=130)
for camp,col,name in (('t5_strong_200k','#2a78d6','strong, 200k spheres'),('t5_weak_12k','#eb6834','weak, 12.5k spheres/rank')):
    R=rows(camp); nps=sorted(set(r['np'] for r in R)); t1=agg([r for r in R if r['np']==1],'loop')
    e=[]; em=[]
    for n in nps:
        l=agg([r for r in R if r['np']==n and r.get('loop')],'loop')
        e.append((t1[0]/l[0])/(n if 'strong' in camp else 1)); em.append((t1[2]/l[2])/(n if 'strong' in camp else 1))
    ax.plot(nps,e,'-o',color=col,lw=2,ms=6,label=name+' (mean)'); ax.plot(nps,em,'--',color=col,lw=1.5,alpha=0.6,label=name+' (min-based)')
ax.set_xscale('log',base=2); ax.set_xticks([1,2,4,8,16]); ax.set_xticklabels(['1','2','4','8','16'])
ax.set_ylim(0,1.15); ax.set_xlabel('MPI ranks (cpus 16..16+n-1, SMT siblings of busy cores)'); ax.set_ylabel('parallel efficiency')
ax.set_title('lmp_release scaling on the settled bed (shared node)',loc='left',fontsize=11)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.grid(axis='y',color='#e6e6e3',lw=0.8); ax.legend(frameon=False,fontsize=8)
fig.tight_layout(); fig.savefig('/media/storage/LIGGGHTS-PUBLIC-v6/audit/plots/perf/scaling_efficiency.png')
