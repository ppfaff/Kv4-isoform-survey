"""(1) Lamp5: per-cell Kv4 / KChIP / DPP family totals sorted by Kv4, with the Lamp5 cluster identity strip.
(2) A-current complex index: each family's CPM divided by its mean over all neurons, then the three multiplied
    per cell; per-subclass fraction of 'depleted' cells and which family is limiting."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip'
d = pd.read_pickle(HIP_PKL).rename(columns={'subclass_label': 'subclass', 'cluster_label': 'cluster'})
d['cls'] = d['subclass'].map(cls_of)
meta = d[['sample_name', 'subclass', 'cluster', 'cls', 'num_processed']].drop_duplicates().set_index('sample_name')
gr = d[['sample_name', 'gene', 'gene_reads']].drop_duplicates().pivot(index='sample_name', columns='gene', values='gene_reads').reindex(meta.index).fillna(0.0)
cpm = gr.div(meta['num_processed'], axis=0) * 1e6
FAMILY = {'Kv4': ['Kcnd1', 'Kcnd2', 'Kcnd3'], 'KChIP': ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'], 'DPP': ['Dpp6', 'Dpp10']}
fam = meta.copy()
for f, genes in FAMILY.items():
    fam[f] = cpm[genes].sum(axis=1)
fam = fam[fam['cls'] != 'NonN']
MARK = {'Kv4': ('o', '#222222'), 'KChIP': ('D', '#1B9E77'), 'DPP': ('^', '#D62728')}

# ------------------------------------------------------------------ 1. Lamp5
L = fam[fam['subclass'] == 'Lamp5'].sort_values('Kv4', ascending=False).copy()
L['x'] = np.arange(len(L), dtype=float)
clusters = L['cluster'].value_counts().index.tolist()
CCOL = dict(zip(clusters, [cm.tab10(i % 10) for i in range(len(clusters))]))

import matplotlib.patches as mpatches
def lamp5_panel(ax, sax, L, title, bounds=None):
    for f, off in [('Kv4', -0.25), ('KChIP', 0.0), ('DPP', 0.25)]:
        mk, c = MARK[f]; m = L[f] > 0
        ax.scatter(L.loc[m, 'x'] + off, L.loc[m, f], s=26, marker=mk, facecolor=c, edgecolor='black', linewidth=0.2, alpha=0.75, label=f, zorder=3)
        z = L[L[f] == 0]
        if len(z):
            yz = {'Kv4': 0.42, 'KChIP': 0.33, 'DPP': 0.26}[f]
            ax.scatter(z['x'], np.full(len(z), yz), marker='x', color=c, s=18, linewidth=0.8, zorder=3, label=f'no {f} ({len(z)})')
    if bounds is None:
        nK = int((L['Kv4'] == 0).sum())
        ax.axvline(len(L) - nK - 0.5, color='red', linestyle='--', linewidth=1.2)
    else:
        for c, x0, x1 in bounds:
            ax.axvline(x0 - 0.5, color='#BBBBBB', linewidth=0.6)
            for f in ['KChIP', 'Kv4', 'DPP']:
                med = L.loc[L['cluster'] == c, f].median()
                ax.plot([x0 - 0.5, x1 + 0.5], [max(med, 0.2)] * 2, color=MARK[f][1], linewidth=2.2, alpha=0.9, zorder=4)
    ax.set_yscale('log'); ax.set_ylim(0.2, None); ax.set_xlim(-1, len(L)); ax.set_xticks([])
    ax.set_ylabel('reads per million (CPM)', fontsize=11)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    ax.legend(loc='lower left', fontsize=9, ncol=2, framealpha=0.9)
    ax.set_title(title, fontsize=12, loc='left')
    sax.bar(L['x'], 1.0, width=1.0, color=[CCOL[c] for c in L['cluster']], edgecolor='none')
    sax.set_xlim(-1, len(L)); sax.set_ylim(0, 1); sax.set_xticks([]); sax.set_yticks([])
    for sp in sax.spines.values(): sp.set_visible(False)
    if bounds is not None:
        for c, x0, x1 in bounds:
            sax.text((x0 + x1) / 2, -0.3, f'{c} ({int(x1 - x0 + 1)})', rotation=90, ha='center', va='top', fontsize=8.5, color=CCOL[c], fontweight='bold')

fig = plt.figure(figsize=(20, 17))
outer = gridspec.GridSpec(2, 1, hspace=0.35, top=0.93, bottom=0.12, left=0.06, right=0.98)
g1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0], height_ratios=[6, 0.45], hspace=0.04)
g2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1], height_ratios=[6, 0.45], hspace=0.04)
lamp5_panel(fig.add_subplot(g1[0]), fig.add_subplot(g1[1]), L,
            f'A. Lamp5 (n = {len(L)}): cells ordered by total Kv4 (high to low); red dashed line = start of the zero-Kv4 cells')
# B: grouped by cluster (clusters ordered by median KChIP), within cluster ordered by Kv4
cl_order = L.groupby('cluster')['KChIP'].median().sort_values(ascending=False).index.tolist()
x, pieces, bounds = 0.0, [], []
for c in cl_order:
    b = L[L['cluster'] == c].sort_values('Kv4', ascending=False).copy()
    b['x'] = np.arange(len(b)) + x; pieces.append(b); bounds.append((c, x, x + len(b) - 1)); x += len(b)
LB = pd.concat(pieces)
lamp5_panel(fig.add_subplot(g2[0]), fig.add_subplot(g2[1]), LB,
            'B. same cells grouped by Lamp5 cluster (clusters ordered by median KChIP, high to low); within cluster ordered by Kv4; '
            'thick bars = cluster medians (green KChIP, black Kv4, red DPP)', bounds=bounds)
fig.suptitle('Lamp5 interneurons: total Kv4, KChIP and DPP per cell — Allen SSv4 hippocampal formation\n'
             'crosses along the bottom mark cells with zero reads for that family', fontsize=13)
fig.savefig(f'{OUT}/hip_lamp5_family_totals.png', dpi=170); plt.close(fig)

# per-cluster breakdown within Lamp5
rows = []
for c, g in L.groupby('cluster'):
    rows.append(dict(cluster=c, n=len(g), **{f'zero_{f}': int((g[f] == 0).sum()) for f in FAMILY},
                     **{f'median_{f}': g[f].median() for f in FAMILY}))
LC = pd.DataFrame(rows).sort_values('n', ascending=False)
LC.to_csv(f'{OUT}/hip_lamp5_by_cluster.csv', index=False)
print('Lamp5 by cluster:'); print(LC.round(1).to_string(index=False))
# co-occurrence of zeros within Lamp5
zK, zC, zD = L['Kv4'] == 0, L['KChIP'] == 0, L['DPP'] == 0
print(f'\nLamp5 zeros: Kv4 {zK.sum()}, KChIP {zC.sum()}, DPP {zD.sum()}; Kv4&KChIP {(zK & zC).sum()}, Kv4&DPP {(zK & zD).sum()}, KChIP&DPP {(zC & zD).sum()}, all three {(zK & zC & zD).sum()}')
print('expected Kv4&KChIP if independent:', round(zK.mean() * zC.mean() * len(L), 1))
print('median KChIP in zero-Kv4 vs Kv4+ Lamp5:', round(L.loc[zK, 'KChIP'].median(), 1), round(L.loc[~zK, 'KChIP'].median(), 1))
print('median DPP   in zero-Kv4 vs Kv4+ Lamp5:', round(L.loc[zK, 'DPP'].median(), 1), round(L.loc[~zK, 'DPP'].median(), 1))

# ------------------------------------------------------------------ 2. A-current complex index
means = {f: fam[f].mean() for f in FAMILY}
print('\nfamily means over all neurons (CPM):', {k: round(v, 1) for k, v in means.items()})
for f in FAMILY:
    fam[f + '_n'] = fam[f] / means[f]
fam['index'] = fam['Kv4_n'] * fam['KChIP_n'] * fam['DPP_n']
fam['limiting'] = fam[['Kv4_n', 'KChIP_n', 'DPP_n']].idxmin(axis=1).str.replace('_n', '')
THR = 0.1
ncells = fam.groupby('subclass').size()
subs = [s for s in SUBCLASS_ORDER if ncells.get(s, 0) >= 15]
rows = []
for s in subs:
    g = fam[fam['subclass'] == s]
    dep = g[g['index'] < THR]
    r = dict(subclass=s, cls=cls_of(s), n=len(g), median_index=g['index'].median(),
             pct_index_zero=100 * (g['index'] == 0).mean(), pct_index_lt_0_01=100 * (g['index'] < 0.01).mean(),
             pct_depleted=100 * len(dep) / len(g))
    for f in FAMILY:
        r[f'pct_depleted_{f}_limiting'] = 100 * (dep['limiting'] == f).sum() / len(g)
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(f'{OUT}/hip_acurrent_index_by_subclass.csv', index=False)
fam[['subclass', 'cluster', 'cls', 'Kv4', 'KChIP', 'DPP', 'Kv4_n', 'KChIP_n', 'DPP_n', 'index', 'limiting']].to_csv(f'{OUT}/hip_acurrent_index_per_cell.csv.gz', index=True)
print(); print(R.round(1).to_string(index=False))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(17, 9), gridspec_kw={'width_ratios': [1.4, 1]}, sharey=True)
y = np.arange(len(subs))[::-1]
# left: distribution of log10 index per subclass (10-90 and 25-75 bands + median)
for yi, s in zip(y, subs):
    v = fam.loc[fam['subclass'] == s, 'index']
    v = np.log10(v.clip(lower=1e-4))
    q = v.quantile([0.1, 0.25, 0.5, 0.75, 0.9])
    c = CLS_COLOR[cls_of(s)]
    a1.plot([q[0.1], q[0.9]], [yi, yi], color=c, linewidth=2, alpha=0.4)
    a1.plot([q[0.25], q[0.75]], [yi, yi], color=c, linewidth=7, alpha=0.8, solid_capstyle='butt')
    a1.plot(q[0.5], yi, 'o', color='white', markeredgecolor='black', markersize=6, zorder=4)
    a1.text(-4.15, yi, f'{int((v <= -4).sum())}', va='center', ha='right', fontsize=8, color='#777')
a1.axvline(np.log10(THR), color='red', linestyle='--', linewidth=1.2)
a1.axvline(0, color='#999', linestyle=':', linewidth=1)
a1.set_xlim(-4.6, 2.2)
a1.set_xticks([-4, -3, -2, -1, 0, 1, 2]); a1.set_xticklabels(['0 / ≤10⁻⁴', '10⁻³', '10⁻²', '10⁻¹', '1', '10', '100'])
a1.set_xlabel('A-current complex index = (Kv4/mean) × (KChIP/mean) × (DPP/mean), log scale', fontsize=11)
a1.set_yticks(y); a1.set_yticklabels(subs, fontsize=10)
for lbl, s in zip(a1.get_yticklabels(), subs): lbl.set_color(CLS_COLOR[cls_of(s)])
a1.set_title('per-subclass distribution (10–90%, 25–75%, median); grey n = cells at zero', fontsize=11)
for sp in ['top', 'right']: a1.spines[sp].set_visible(False)
# right: % depleted stacked by limiting family
left = np.zeros(len(subs))
FCOL = {'Kv4': '#222222', 'KChIP': '#1B9E77', 'DPP': '#D62728'}
for f in FAMILY:
    vals = R.set_index('subclass').loc[subs, f'pct_depleted_{f}_limiting'].values
    a2.barh(y, vals, left=left, color=FCOL[f], height=0.7, label=f'{f} is the limiting family')
    left += vals
for yi, s, tot, n in zip(y, subs, left, R.set_index('subclass').loc[subs, 'n']):
    a2.text(tot + 0.5, yi, f'{tot:.1f}%  (n={n:,})', va='center', fontsize=8.5)
a2.set_xlim(0, max(left) * 1.35); a2.set_xlabel(f'% of cells with index < {THR}', fontsize=11)
a2.set_title(f'cells depleted for A-current subunits (index < {THR}), by limiting family', fontsize=11)
a2.legend(loc='lower right', fontsize=10, frameon=False)
for sp in ['top', 'right']: a2.spines[sp].set_visible(False)
fig.suptitle('Hippocampal formation: A-current complex index per cell\n'
             f'family CPM normalised to its mean over all {len(fam):,} neurons (Kv4 {means["Kv4"]:.0f}, KChIP {means["KChIP"]:.0f}, DPP {means["DPP"]:.0f} CPM), then multiplied; '
             'a cell at the neuron-average for all three scores 1', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig(f'{OUT}/hip_acurrent_index_by_subclass.png', dpi=170)
print('done')
