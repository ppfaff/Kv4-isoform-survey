"""Hippocampal per-cell Kv4 pore-subunit signal, clustered by subclass (the cortex v5 design), plus the
Kv4 / KChIP / DPP family-total comparison. Normalisation: reads per million reads processed by salmon
(num_processed) — the per-cell depth available for this dataset."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from scipy.stats import spearmanr, mannwhitneyu
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip'
d = load_hip()
MIN_CELLS = 15
ncells = d.groupby('subclass')['sample_name'].nunique()
keep = [s for s in SUBCLASS_ORDER if ncells.get(s, 0) >= MIN_CELLS]
d = d[d['subclass'].isin(keep)]

VARIANTS = ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L']
COLOR = {'Kv4.1': '#4FA8D8', 'Kv4.2': '#E8822D', 'Kv4.3S': '#B9A3D9', 'Kv4.3L': '#7B4FA3'}
OFFSET = {'Kv4.1': -0.27, 'Kv4.2': -0.09, 'Kv4.3S': 0.09, 'Kv4.3L': 0.27}

dd = d[d['vclass'].isin(VARIANTS)]
percell = dd.groupby(['sample_name', 'vclass'])['reads'].sum().reset_index()
wide = percell.pivot(index='sample_name', columns='vclass', values='reads').reindex(columns=VARIANTS).fillna(0.0)
meta = d[['sample_name', 'subclass', 'cls', 'num_processed', 'percent_mapped']].drop_duplicates().set_index('sample_name')
wide = wide.join(meta, how='right').fillna({v: 0.0 for v in VARIANTS})
for v in VARIANTS:
    wide[v + '_cpm'] = wide[v] / wide['num_processed'] * 1e6
VCPM = [v + '_cpm' for v in VARIANTS]
wide['total_cpm'] = wide[VCPM].sum(axis=1)
wide['total_raw'] = wide[VARIANTS].sum(axis=1)

palette = [cm.tab20(i % 20) if i < 20 else cm.tab20b((i - 20) % 20) for i in range(len(keep))]
TYPE_COLOR = dict(zip(keep, palette))
CLS_ROWS = [('GABA', 'Inhibitory (GABA)'), ('Glut', 'Excitatory / principal (Glut)'), ('NonN', 'Non-neuronal')]
FIG_ROWS = CLS_ROWS[:2]   # non-neuronal (37 cells, all but one with zero Kv4) reported in the table, not plotted

# ------------------------------------------------------------------ figure 1: per-cell, clustered by subclass
fig = plt.figure(figsize=(22, 14))
outer = gridspec.GridSpec(2, 1, hspace=0.45, top=0.91, bottom=0.17, left=0.05, right=0.98, figure=fig,
                          height_ratios=[1, 1.5])
axes = []
for block in range(2):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[block], height_ratios=[4, 0.35], hspace=0.03)
    axes.append(fig.add_subplot(inner[0])); axes.append(fig.add_subplot(inner[1]))

zero_stats = []
for i, (cls_key, cls_label) in enumerate(CLS_ROWS):
    sub = wide[wide['cls'] == cls_key].copy()
    if cls_key not in dict(FIG_ROWS):
        for t, g in sub.groupby('subclass'):
            nz = int((g['total_cpm'] == 0).sum())
            zero_stats.append(dict(cls=cls_key, subclass=t, n=len(g), n_zero_kv4=nz, pct_zero=100 * nz / len(g),
                                   mean_kv4_cpm=g['total_cpm'].mean(), median_kv4_cpm=g['total_cpm'].median()))
        continue
    ax, strip_ax = axes[2 * i], axes[2 * i + 1]
    type_mean = sub.groupby('subclass')['total_cpm'].mean().sort_values(ascending=False)
    order = type_mean.index.tolist()
    x, bounds, pieces = 0.0, [], []
    for t in order:
        block = sub[sub['subclass'] == t].sort_values('total_cpm', ascending=False).copy()
        block['x'] = np.arange(len(block)) + x
        pieces.append(block); bounds.append((t, x, x + len(block) - 1, type_mean[t])); x += len(block)
    sub = pd.concat(pieces); n_total = len(sub)
    ms = 14 if n_total > 2000 else (22 if n_total > 500 else 40)
    for v in VARIANTS:
        vals = sub[v + '_cpm']; m = vals > 0
        ax.scatter(sub.loc[m, 'x'] + OFFSET[v], vals[m], s=ms, color=COLOR[v], edgecolor='black',
                   linewidth=0.15, alpha=0.7, label=v, zorder=3)
    miss = sub[sub['total_cpm'] == 0]
    if len(miss):
        ax.scatter(miss['x'], np.full(len(miss), 0.42), marker='x', color='#888888', s=ms * 0.6,
                   linewidth=0.8, zorder=3, label='no Kv4 detected')
    for t, x0, x1, m in bounds:
        ax.plot([x0 - 0.5, x1 + 0.5], [m, m], color='black', linewidth=1.6, alpha=0.6, zorder=2)
        ax.axvline(x0 - 0.5, color='#CCCCCC', linewidth=0.5, zorder=0.5)
        n = x1 - x0 + 1; nz = int((sub[sub['subclass'] == t]['total_cpm'] == 0).sum())
        zero_stats.append(dict(cls=cls_key, subclass=t, n=n, n_zero_kv4=nz, pct_zero=100 * nz / n,
                               mean_kv4_cpm=m, median_kv4_cpm=sub[sub['subclass'] == t]['total_cpm'].median()))
    ax.set_yscale('log'); ax.set_ylim(0.3, None); ax.set_xlim(-1, n_total); ax.set_xticks([])
    ax.set_title(f'{cls_label}  (n = {n_total:,} cells, {len(order)} subclasses, blocks ordered by subclass-mean Kv4)',
                 fontsize=12, fontweight='bold', loc='left')
    ax.set_ylabel('Kv4 reads per\nmillion reads', fontsize=10)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    strip_ax.bar(sub['x'], height=1.0, width=1.0, color=[TYPE_COLOR[t] for t in sub['subclass']], edgecolor='none')
    strip_ax.set_xlim(-1, n_total); strip_ax.set_ylim(0, 1); strip_ax.set_xticks([]); strip_ax.set_yticks([])
    for sp in strip_ax.spines.values(): sp.set_visible(False)
    for t, x0, x1, m in bounds:
        strip_ax.text((x0 + x1) / 2, -1.6, f'{t} ({int(x1 - x0 + 1):,})', rotation=90 if (x1 - x0) < n_total * 0.06 else 0,
                      fontsize=9, ha='center', va='top', color=TYPE_COLOR[t], fontweight='bold')
axes[-2].set_xlabel('individual cells, grouped into subclass blocks (block = subclass-mean Kv4, high to low); '
                    'within block, high to low Kv4 left to right', fontsize=10)
h, l = axes[0].get_legend_handles_labels()
axes[0].legend(h, l, loc='upper center', fontsize=13, framealpha=0.9, ncol=5, markerscale=2.5)
fig.suptitle('Hippocampal formation: Kv4 pore-subunit signal per cell, clustered by subclass\n'
             'black step = per-subclass mean; gray X = cell with no detectable Kv4 read\n'
             f'Allen SSv4 (Yao 2021), {len(wide):,} cells; non-neuronal cells (n = 37, 36 with zero Kv4) omitted', fontsize=12)
fig.savefig(f'{OUT}/hip_percell_kv4_by_subclass.png', dpi=170)
Z = pd.DataFrame(zero_stats).sort_values(['cls', 'mean_kv4_cpm'], ascending=[True, False])
Z.to_csv(f'{OUT}/hip_kv4_zero_by_subclass.csv', index=False)
print(Z.round(2).to_string(index=False))

# ------------------------------------------------------------------ figure 2: family totals per cell
FAMILY = {'Kv4': ['Kcnd1', 'Kcnd2', 'Kcnd3'], 'KChIP': ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'], 'DPP': ['Dpp6', 'Dpp10']}
MARK = {'Kv4': ('o', '#222222'), 'KChIP': ('D', '#1B9E77'), 'DPP': ('^', '#D62728')}
FOFF = {'Kv4': -0.22, 'KChIP': 0.0, 'DPP': 0.22}
gr = d[['sample_name', 'gene', 'gene_reads']].drop_duplicates().pivot(index='sample_name', columns='gene', values='gene_reads').fillna(0.0)
fam = meta.copy()
for f, genes in FAMILY.items():
    fam[f + '_cpm'] = gr.reindex(index=fam.index, columns=genes).fillna(0).sum(axis=1) / fam['num_processed'] * 1e6

fig, axs = plt.subplots(3, 1, figsize=(22, 12), gridspec_kw={'height_ratios': [1, 1.6, 0.5]})
rep = []
for ax, (cls_key, cls_label) in zip(axs, CLS_ROWS):
    sub = fam[fam['cls'] == cls_key].sort_values('Kv4_cpm', ascending=False).copy()
    sub['x'] = np.arange(len(sub), dtype=float)
    ms = 10 if len(sub) > 2000 else 24
    for f in FAMILY:
        mk, c = MARK[f]; vals = sub[f + '_cpm']; m = vals > 0
        ax.scatter(sub.loc[m, 'x'] + FOFF[f], vals[m], s=ms, marker=mk, facecolor=c, edgecolor='black',
                   linewidth=0.15, alpha=0.65, label=f, zorder=3)
    ax.set_yscale('log'); ax.set_ylim(0.3, None); ax.set_xlim(-1, len(sub)); ax.set_xticks([])
    ax.set_title(f'{cls_label}  (n = {len(sub):,} cells)', fontsize=12, fontweight='bold', loc='left')
    ax.set_ylabel('reads per million\n(CPM)', fontsize=10)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    rkk = spearmanr(sub['Kv4_cpm'], sub['KChIP_cpm'])[0]; rkd = spearmanr(sub['Kv4_cpm'], sub['DPP_cpm'])[0]
    rcd = spearmanr(sub['KChIP_cpm'], sub['DPP_cpm'])[0]
    z = sub[sub['Kv4_cpm'] == 0]
    rep.append(dict(cls=cls_key, n=len(sub), rho_Kv4_KChIP=rkk, rho_Kv4_DPP=rkd, rho_KChIP_DPP=rcd,
                    n_zero_Kv4=len(z), zero_Kv4_also_zero_KChIP=int((z['KChIP_cpm'] == 0).sum()),
                    zero_Kv4_also_zero_DPP=int((z['DPP_cpm'] == 0).sum()),
                    n_zero_KChIP=int((sub['KChIP_cpm'] == 0).sum()), n_zero_DPP=int((sub['DPP_cpm'] == 0).sum())))
axs[-1].set_xlabel('individual cells, ordered by total Kv4 signal (highest to lowest), left to right', fontsize=10)
axs[0].legend(loc='upper right', fontsize=12, framealpha=0.9, ncol=3, markerscale=2.5)
fig.suptitle('Hippocampal formation: total Kv4, KChIP and DPP family expression per cell (CPM), cells ordered by Kv4\n'
             f'Allen SSv4 (Yao 2021), {len(fam):,} cells', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig(f'{OUT}/hip_kv4_kchip_dpp_total.png', dpi=170)
R = pd.DataFrame(rep); R.to_csv(f'{OUT}/hip_family_cooccurrence.csv', index=False)
print(); print(R.round(3).to_string(index=False))

# depth check: does zero-Kv4 track low depth?
for cls_key, _ in CLS_ROWS:
    sub = wide[wide['cls'] == cls_key]
    z, nz = sub[sub['total_cpm'] == 0]['num_processed'], sub[sub['total_cpm'] > 0]['num_processed']
    if len(z) >= 5:
        print(f'{cls_key}: depth median zero-Kv4 {z.median():,.0f} vs expressing {nz.median():,.0f}; MWU p={mannwhitneyu(z, nz).pvalue:.2g}')
