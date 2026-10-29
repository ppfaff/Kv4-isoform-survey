"""Per-cell, gene-level view of each Kv4-complex family (KChIP: Kcnip1-4; DPP: Dpp6/Dpp10), clustered by
subclass — same design as the per-cell Kv4 figure. Gene totals (all transcripts), CPM by num_processed."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip'
DS = sys.argv[1] if len(sys.argv) > 1 else 'hip'
d, ORDER, REGION, SOURCE, CLS = load_dataset(DS)
OUT = f'{OUT}/{DS}'
import os; os.makedirs(OUT, exist_ok=True)
ncells = d.groupby('subclass')['sample_name'].nunique()
keep = [s for s in ORDER if ncells.get(s, 0) >= (15 if DS == 'hip' else 5)]
d = d[d['subclass'].isin(keep)]

meta = d[['sample_name', 'subclass', 'cls', 'num_processed']].drop_duplicates().set_index('sample_name')
gr = d[['sample_name', 'gene', 'gene_reads']].drop_duplicates().pivot(index='sample_name', columns='gene', values='gene_reads').fillna(0.0)
gr = gr.reindex(meta.index).fillna(0.0)
cpm = gr.div(meta['num_processed'], axis=0) * 1e6

palette = [cm.tab20(i % 20) if i < 20 else cm.tab20b((i - 20) % 20) for i in range(len(keep))]
TYPE_COLOR = dict(zip(keep, palette))
CLS_ROWS = [('GABA', 'Inhibitory (GABA)'), ('Glut', 'Excitatory / principal (Glut)')]
NONN_N = int((meta['cls'] == 'NonN').sum())

FAMILIES = {
    'KChIP': dict(genes=['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'],
                  color={'Kcnip1': HEX['purple'], 'Kcnip2': HEX['sky_blue'], 'Kcnip3': HEX['bluish_green'], 'Kcnip4': HEX['orange']},
                  offset={'Kcnip1': -0.27, 'Kcnip2': -0.09, 'Kcnip3': 0.09, 'Kcnip4': 0.27},
                  label={'Kcnip1': 'KChIP1 (Kcnip1)', 'Kcnip2': 'KChIP2 (Kcnip2)', 'Kcnip3': 'KChIP3 (Kcnip3)', 'Kcnip4': 'KChIP4 (Kcnip4)'}),
    'DPP':   dict(genes=['Dpp6', 'Dpp10'],
                  color={'Dpp6': HEX['purple'], 'Dpp10': HEX['orange']},
                  offset={'Dpp6': -0.18, 'Dpp10': 0.18},
                  label={'Dpp6': 'DPP6 (Dpp6)', 'Dpp10': 'DPP10 (Dpp10)'}),
    'Kv4':   dict(genes=['Kcnd1', 'Kcnd2', 'Kcnd3'],
                  color={'Kcnd1': HEX['sky_blue'], 'Kcnd2': HEX['orange'], 'Kcnd3': HEX['purple']},
                  offset={'Kcnd1': -0.22, 'Kcnd2': 0.0, 'Kcnd3': 0.22},
                  label={'Kcnd1': 'Kv4.1 (Kcnd1)', 'Kcnd2': 'Kv4.2 (Kcnd2)', 'Kcnd3': 'Kv4.3 (Kcnd3)'}),
}

summary = []
for fam, spec in FAMILIES.items():
    genes = spec['genes']
    tot = cpm[genes].sum(axis=1)
    fig = plt.figure(figsize=(22, 14))
    outer = gridspec.GridSpec(2, 1, hspace=0.45, top=0.91, bottom=0.17, left=0.05, right=0.98, figure=fig, height_ratios=[1, 1.5])
    axes = []
    for block in range(2):
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[block], height_ratios=[4, 0.35], hspace=0.03)
        axes.append(fig.add_subplot(inner[0])); axes.append(fig.add_subplot(inner[1]))
    for i, (cls_key, cls_label) in enumerate(CLS_ROWS):
        ax, strip_ax = axes[2 * i], axes[2 * i + 1]
        idx = meta.index[meta['cls'] == cls_key]
        sub = pd.concat([meta.loc[idx], cpm.loc[idx, genes]], axis=1)
        sub['total'] = tot.loc[idx]
        type_mean = sub.groupby('subclass')['total'].mean().sort_values(ascending=False)
        order = type_mean.index.tolist()
        x, bounds, pieces = 0.0, [], []
        for t in order:
            block = sub[sub['subclass'] == t].sort_values('total', ascending=False).copy()
            block['x'] = np.arange(len(block)) + x
            pieces.append(block); bounds.append((t, x, x + len(block) - 1, type_mean[t])); x += len(block)
        sub = pd.concat(pieces); n_total = len(sub)
        ms = 14 if n_total > 2000 else (22 if n_total > 500 else 34)
        for g in genes:
            m = sub[g] > 0
            ax.scatter(sub.loc[m, 'x'] + spec['offset'][g], sub.loc[m, g], s=ms, color=spec['color'][g],
                       edgecolor='black', linewidth=0.15, alpha=0.7, label=spec['label'][g], zorder=3)
        miss = sub[sub['total'] == 0]
        if len(miss):
            ax.scatter(miss['x'], np.full(len(miss), 0.42), marker='x', color='#888888', s=ms * 0.6,
                       linewidth=0.8, zorder=3, label=f'no {fam} detected')
        for t, x0, x1, m in bounds:
            ax.plot([x0 - 0.5, x1 + 0.5], [m, m], color='black', linewidth=1.6, alpha=0.6, zorder=2)
            ax.axvline(x0 - 0.5, color='#CCCCCC', linewidth=0.5, zorder=0.5)
            blk = sub[sub['subclass'] == t]
            row = dict(family=fam, cls=cls_key, subclass=t, n=len(blk), n_zero=int((blk['total'] == 0).sum()),
                       pct_zero=100 * (blk['total'] == 0).mean(), mean_cpm=m, median_cpm=blk['total'].median())
            for g in genes:
                row[f'pct_zero_{g}'] = 100 * (blk[g] == 0).mean()
            summary.append(row)
        ax.set_yscale('log'); ax.set_ylim(0.3, None); ax.set_xlim(-1, n_total); ax.set_xticks([])
        ax.set_title(f'{cls_label}  (n = {n_total:,} cells, {len(order)} subclasses, blocks ordered by subclass-mean {fam})',
                     fontsize=12, fontweight='bold', loc='left')
        ax.set_ylabel(f'{fam} reads per\nmillion reads', fontsize=10)
        for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
        strip_ax.bar(sub['x'], height=1.0, width=1.0, color=[TYPE_COLOR[t] for t in sub['subclass']], edgecolor='none')
        strip_ax.set_xlim(-1, n_total); strip_ax.set_ylim(0, 1); strip_ax.set_xticks([]); strip_ax.set_yticks([])
        for sp in strip_ax.spines.values(): sp.set_visible(False)
        for t, x0, x1, m in bounds:
            strip_ax.text((x0 + x1) / 2, -1.6, f'{t} ({int(x1 - x0 + 1):,})', rotation=90 if (x1 - x0) < n_total * (0.06 if DS == 'hip' else 0.12) else 0,
                          fontsize=9, ha='center', va='top', color=TYPE_COLOR[t], fontweight='bold')
    axes[-2].set_xlabel(f'individual cells, grouped into subclass blocks (block = subclass-mean {fam}, high to low); '
                        f'within block, high to low {fam} left to right', fontsize=10)
    h, l = axes[0].get_legend_handles_labels()
    axes[0].legend(h, l, loc='upper center', fontsize=13, framealpha=0.9, ncol=len(genes) + 1, markerscale=2.5)
    fig.suptitle(f'{REGION}: {fam} family expression per cell (gene level, all transcripts), clustered by type\n'
                 f'black step = per-type mean; gray X = cell with no detectable {fam} read\n'
                 f'{SOURCE}, {len(meta) - NONN_N:,} neurons (non-neuronal cells omitted)', fontsize=12)
    fig.savefig(f'{OUT}/{DS}_percell_{fam}_genes_by_type.png', dpi=170)
    plt.close(fig)

S = pd.DataFrame(summary)
S.to_csv(f'{OUT}/{DS}_family_zero_by_type.csv', index=False)
pd.set_option('display.width', 250)
for fam in FAMILIES:
    x = S[S.family == fam].sort_values(['cls', 'mean_cpm'], ascending=[True, False])
    cols = ['cls', 'subclass', 'n', 'n_zero', 'pct_zero', 'mean_cpm', 'median_cpm'] + [c for c in x.columns if c.startswith('pct_zero_')]
    print(f'\n=== {fam} ==='); print(x[cols].round(1).to_string(index=False))
