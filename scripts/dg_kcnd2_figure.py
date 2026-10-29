import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

E = pd.read_csv('/home/claude/kv4/hip/markers/allen_exon_markers.csv', index_col=0)
C = pd.read_csv('/home/claude/kv4/hip/markers/allen_markers.csv', index_col=0)
dg = pd.read_csv('/home/claude/kv4/hip/markers/dg_cells_markers.csv', index_col=0)
dep = dg['num_processed']
def cpm(df, g): return df.loc[dg.index, g] / dep * 1e6
ex = {g: cpm(E, g) for g in ['Kcnd2', 'Kcnip4', 'Dpp6']}
inn = {g: ((C.loc[dg.index, g] - E.loc[dg.index, g]).clip(lower=0) / dep * 1e6) for g in ['Kcnd2', 'Kcnip4', 'Dpp6']}
GC = {'Kv4-null': '#D62728', 'Kv4-low': '#F0B429', 'Kv4-high': '#2060B0'}

fig, axes = plt.subplots(1, 4, figsize=(22, 5.8))
for ax, g in zip(axes[:3], ['Kcnd2', 'Kcnip4', 'Dpp6']):
    x = ex[g].clip(lower=0.3); y = inn[g].clip(lower=0.3)
    for grp in ['Kv4-high', 'Kv4-low', 'Kv4-null']:
        m = dg.grp == grp
        ax.scatter(x[m], y[m], s=7, color=GC[grp], alpha=0.5, label=f'{grp} (n={m.sum()})', edgecolor='none')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel(f'{g} exonic (spliced mRNA), CPM — Allen exon-only counts', fontsize=10)
    ax.set_ylabel(f'{g} intronic (nascent pre-mRNA), CPM', fontsize=10)
    rho = spearmanr(ex[g], inn[g])[0]
    ax.set_title(f'{g}: exonic vs intronic per DG cell  (ρ = {rho:.2f})', fontsize=11, fontweight='bold')
    ax.axvline(0.42, color='#999', linewidth=0.6, linestyle=':')
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
axes[0].legend(fontsize=9, loc='lower right', markerscale=3)
axes[0].text(0.33, 2000, 'exonic = 0\n(these cells)', fontsize=8, color='#D62728', ha='left', va='top')

# panel 4: maturity / identity markers by group
mk = ['Dcx', 'Calb2', 'Sox11', 'Neurod1', 'Calb1', 'Prox1', 'C1ql2']
ax = axes[3]
w = 0.26
for k, grp in enumerate(['Kv4-null', 'Kv4-low', 'Kv4-high']):
    sub = dg[dg.grp == grp]
    pct = [(sub[m] > 0).mean() * 100 for m in mk]
    ax.bar(np.arange(len(mk)) + (k - 1) * w, pct, width=w, color=GC[grp], label=grp)
ax.set_xticks(np.arange(len(mk))); ax.set_xticklabels(mk, rotation=30, ha='right', fontsize=10)
ax.set_ylabel('% of DG cells with the marker detected', fontsize=10)
ax.set_title('maturity / identity markers by Kv4 group', fontsize=11, fontweight='bold')
ax.axvspan(-0.5, 3.5, color='#F5E6C8', alpha=0.4, zorder=0); ax.text(1.5, 103, 'immature / adult-born', ha='center', fontsize=9, color='#8a6d1f')
ax.axvspan(3.5, 6.5, color='#D6E4F5', alpha=0.4, zorder=0); ax.text(5, 103, 'mature granule cell', ha='center', fontsize=9, color='#2060B0')
ax.set_ylim(0, 110); ax.legend(fontsize=9, loc='center left')
for sp in ['top', 'right']: ax.spines[sp].set_visible(False)

fig.suptitle('DG granule cells (n = 2,458): Kv4-null cells are mature cells that transcribe Kcnd2 but lack spliced Kcnd2 mRNA in the soma\n'
             'Kv4 group from our salmon quantification (null = 0 reads on Kcnd1/2/3; low/high = below/above median Kcnd2); '
             'exonic and intronic counts from the Allen RSEM matrices (exon-only vs combined)', fontsize=11.5)
fig.tight_layout(rect=[0, 0, 1, 0.9])
fig.savefig('/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip/dg_kcnd2_exon_vs_intron.png', dpi=170)
print('done')
