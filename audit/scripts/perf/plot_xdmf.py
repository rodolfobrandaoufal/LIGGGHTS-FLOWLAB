import json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
L='/media/storage/LIGGGHTS-PUBLIC-v6/audit/logs/perf/'
rows=[json.loads(l) for l in open(L+'t4c_io.jsonl') if 'xdmf_growth' in l]
C={'hdf5':'#2a78d6','custombin':'#eb6834'}; N={'hdf5':'dump hdf5 (+ XDMF rewrite)','custombin':'dump custom .bin'}
fig,ax=plt.subplots(figsize=(7,4),dpi=130)
for fmt in ('hdf5','custombin'):
    for k,r in enumerate([r for r in rows if r['fmt']==fmt]):
        c=r['cpu']; x=[c[i][0] for i in range(1,len(c))]; y=[1e3*(c[i][1]-c[i-1][1])/100 for i in range(1,len(c))]
        ax.plot(x,y,color=C[fmt],lw=2,alpha=0.9 if k==0 else 0.35,label=N[fmt] if k==0 else None)
ax.set_xlabel('dump index (one dump per step, 64 atoms)'); ax.set_ylabel('wall time per dump [ms]')
ax.set_title('Per-dump cost grows linearly with dump count (F-25)',loc='left',fontsize=11)
for s in ('top','right'): ax.spines[s].set_visible(False)
ax.grid(axis='y',color='#e6e6e3',lw=0.8); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig('/media/storage/LIGGGHTS-PUBLIC-v6/audit/plots/perf/xdmf_growth.png')
