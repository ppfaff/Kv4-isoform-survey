import pandas as pd, numpy as np
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = pd.read_pickle('/home/claude/kv4/iso/d.pkl')
counts = pd.read_pickle('/home/claude/kv4/counts.pkl')

FAMILY = {
    'Kv4':   ['Kcnd1', 'Kcnd2', 'Kcnd3'],
    'KChIP': ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'],
    'DPP':   ['Dpp6', 'Dpp10'],
}
MARK = {'Kv4': ('o', '#222222'), 'KChIP': ('D', '#1B9E77'), 'DPP': ('^', '#D62728')}
OFFSET = {'Kv4': -0.22, 'KChIP': 0.0, 'DPP': 0.22}

# per (sample_name, gene) total reads — already de-duplicated across that gene's transcripts
gene_reads = d[['sample_name', 'gene', 'gene_reads']].drop_duplicates()
gr = gene_reads.pivot(index='sample_name', columns='gene', values='gene_reads').fillna(0.0)

meta = d[['sample_name', 'type', 'primary_type_id', 'cls']].drop_duplicates().set_index('sample_name')

depth = counts.sum(axis=0)
depth.index.name = 'sample_name'

wide = meta.copy()
for fam, genes in FAMILY.items():
    wide[fam + '_raw'] = gr.reindex(index=wide.index, columns=genes).fillna(0.0).sum(axis=1)
wide['depth'] = depth.reindex(wide.index)
wide = wide.dropna(subset=['depth'])
for fam in FAMILY:
    wide[fam + '_cpm'] = wide[fam + '_raw'] / wide['depth'] * 1e6

CLS_ROWS = [('GABA', 'Inhibitory (GABA)'), ('Glut', 'Excitatory (Glut)'), ('NonN', 'Non-neuronal')]

fig, axes = plt.subplots(3, 1, figsize=(18, 11), sharey=False)
corr_report = {}

for ax, (cls_key, cls_label) in zip(axes, CLS_ROWS):
    sub = wide[wide['cls'] == cls_key].copy()
    if sub.empty:
        ax.set_visible(False)
        continue
    sub = sub.sort_values('Kv4_cpm', ascending=False)
    sub['x'] = np.arange(len(sub), dtype=float)

    for fam in ['Kv4', 'KChIP', 'DPP']:
        marker, color = MARK[fam]
        vals = sub[fam + '_cpm']
        mask = vals > 0
        ax.scatter(sub.loc[mask, 'x'] + OFFSET[fam], vals[mask], s=60, marker=marker,
                   facecolor=color, edgecolor='black', linewidth=0.3, alpha=0.75, label=fam, zorder=3)

    ax.set_yscale('log')
    ax.set_ylim(0.3, None)
    ax.set_xlim(-1, len(sub))
    ax.set_xticks([])
    ax.set_title(f'{cls_label}  (n = {len(sub)} cells)', fontsize=11, fontweight='bold', loc='left')
    ax.set_ylabel('reads per million\ntotal reads (CPM)', fontsize=9)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    rho_kk, _ = spearmanr(sub['Kv4_cpm'], sub['KChIP_cpm'])
    rho_kd, _ = spearmanr(sub['Kv4_cpm'], sub['DPP_cpm'])
    rho_cd, _ = spearmanr(sub['KChIP_cpm'], sub['DPP_cpm'])
    corr_report[cls_key] = (rho_kk, rho_kd, rho_cd, len(sub))

axes[-1].set_xlabel('individual cells, ordered by total Kv4 signal (highest to lowest), left to right', fontsize=10)
axes[0].legend(loc='upper right', fontsize=10, framealpha=0.9, ncol=3)

fig.suptitle('Total Kv4, KChIP, and DPP family expression per cell (normalized CPM)\n'
             'cells ordered by Kv4 signal — do KChIP and DPP decline in the same cells?\n'
             '460-cell preview, Tasic 2016 (GSE71585)', fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig('/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/kv4_kchip_dpp_total.png', dpi=200)

print('Spearman rank correlations (per class, based on CPM):')
print(f'{"class":8s} {"n":>5s} {"Kv4-KChIP":>10s} {"Kv4-DPP":>10s} {"KChIP-DPP":>10s}')
for cls_key, (rkk, rkd, rcd, n) in corr_report.items():
    print(f'{cls_key:8s} {n:5d} {rkk:10.3f} {rkd:10.3f} {rcd:10.3f}')

# zero-overlap check: of the cells with zero Kv4, how many also have zero KChIP / zero DPP?
for cls_key, _ in CLS_ROWS:
    sub = wide[wide['cls'] == cls_key]
    zero_kv4 = sub[sub['Kv4_cpm'] == 0]
    if len(zero_kv4) == 0:
        continue
    print(f'\n{cls_key}: {len(zero_kv4)} cells with zero Kv4 —',
          f'also zero KChIP: {(zero_kv4["KChIP_cpm"]==0).sum()},',
          f'also zero DPP: {(zero_kv4["DPP_cpm"]==0).sum()}')
