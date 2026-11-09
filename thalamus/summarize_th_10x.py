#!/usr/bin/env python3
"""Summarize Kv4-complex expression across thalamic cell types and thalamocortical-loop elements.

Input: tables/th10x_cells.parquet from extract_th_10x.py (per-cell log2 CPM+1 for the gene panel).

Outputs (OUT/tables and OUT/figures):
  th10x_cluster_means.csv       per cluster: n, mean log2, fraction detected, nucleus, loop group
  th10x_loop_summary.csv        per loop element: cell-weighted mean, detection, Kv4 paralog shares,
                                KChIP and DPP partner composition
  th10x_kv4_loop_heatmap.png    A-current + burst genes x loop elements
  th10x_kv4_paralog_share.png   Kcnd1/2/3 share per loop element and per RT / relay subclass
  th10x_percell_kv4_scatter.png Kcnd2 vs Kcnd3 per cell, TRN vs first-order vs higher-order relay

Usage: python3 summarize_th_10x.py [TABLES_DIR] [OUT_DIR] [--min-cells 30]
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from th_common import (KV4, KCHIP, DPP, COMPLEX, BURST, MARKERS, LOOP_ORDER, cluster_to_nucleus,  # noqa: E402
                       loop_group_of_nucleus, subclass_fallback_group, RT_SUBCLASS, read_cells)

ap = argparse.ArgumentParser()
ap.add_argument('tables', nargs='?', default='tables')
ap.add_argument('out', nargs='?', default='.')
ap.add_argument('--min-cells', type=int, default=30)
ap.add_argument('--neurons-only', action='store_true', default=True)
args = ap.parse_args()
os.makedirs(f'{args.out}/tables', exist_ok=True); os.makedirs(f'{args.out}/figures', exist_ok=True)

df = read_cells(f'{args.tables}/th10x_cells')
genes = [g for g in COMPLEX + BURST + MARKERS if g in df.columns]
if args.neurons_only and 'neurotransmitter' in df.columns:
    df = df[df['neurotransmitter'].isin(['Glut', 'GABA', 'Glut-GABA', 'GABA-Glyc'])]
print(f'{len(df):,} neurons')

# ------------------------------------------------ cluster -> nucleus -> loop element
c2n, _ = cluster_to_nucleus()
df['nucleus'] = df['cluster'].map(c2n)
df['loop_group'] = df['nucleus'].map(loop_group_of_nucleus)
fallback = df['loop_group'].isna() | (df['loop_group'] == 'unassigned')
df.loc[fallback, 'loop_group'] = df.loc[fallback, 'subclass'].map(subclass_fallback_group)
df['loop_group'] = df['loop_group'].fillna('unassigned')

lin = (2.0 ** df[genes] - 1.0)                    # back to CPM for means that respect abundance
det = (df[genes] > 0).astype(float)

# ------------------------------------------------ per-cluster table
grp = df.groupby(['class', 'subclass', 'supertype', 'cluster'], observed=True)
cl = grp.size().rename('n_cells').to_frame()
cl = cl.join(lin.groupby([df['class'], df['subclass'], df['supertype'], df['cluster']], observed=True).mean()
             .add_prefix('cpm_'))
cl = cl.join(det.groupby([df['class'], df['subclass'], df['supertype'], df['cluster']], observed=True).mean()
             .add_prefix('det_'))
cl = cl.reset_index()
cl['nucleus'] = cl['cluster'].map(c2n)
cl['loop_group'] = cl['nucleus'].map(loop_group_of_nucleus)
m = cl['loop_group'].isin(['unassigned'])
cl.loc[m, 'loop_group'] = cl.loc[m, 'subclass'].map(subclass_fallback_group)
kv4tot = cl[[f'cpm_{g}' for g in KV4]].sum(axis=1)
for g in KV4:
    cl[f'share_{g}'] = cl[f'cpm_{g}'] / kv4tot
cl['KChIP_per_Kv4'] = cl[[f'cpm_{g}' for g in KCHIP]].sum(axis=1) / kv4tot
cl['DPP_per_Kv4'] = cl[[f'cpm_{g}' for g in DPP]].sum(axis=1) / kv4tot
cl = cl[cl.n_cells >= args.min_cells].sort_values(['loop_group', 'subclass', 'cluster'])
cl.to_csv(f'{args.out}/tables/th10x_cluster_means.csv', index=False)


# ------------------------------------------------ per loop element
def summarize(sub):
    L = 2.0 ** sub[genes] - 1.0
    r = {'n_cells': len(sub)}
    for g in genes:
        r[f'cpm_{g}'] = L[g].mean()
        r[f'det_{g}'] = (sub[g] > 0).mean()
    tot = sum(r[f'cpm_{g}'] for g in KV4)
    for g in KV4:
        r[f'share_{g}'] = r[f'cpm_{g}'] / tot if tot else np.nan
    ktot = sum(r[f'cpm_{g}'] for g in KCHIP); dtot = sum(r[f'cpm_{g}'] for g in DPP)
    for g in KCHIP:
        r[f'KChIP_comp_{g}'] = r[f'cpm_{g}'] / ktot if ktot else np.nan
    for g in DPP:
        r[f'DPP_comp_{g}'] = r[f'cpm_{g}'] / dtot if dtot else np.nan
    r['KChIP_per_Kv4'] = ktot / tot if tot else np.nan
    r['DPP_per_Kv4'] = dtot / tot if tot else np.nan
    # per-cell co-expression: fraction of Kcnd2+ cells, fraction Kcnd3-only, fraction Kv4-null
    k2, k3, k1 = sub['Kcnd2'] > 0, sub['Kcnd3'] > 0, sub['Kcnd1'] > 0
    r['frac_Kcnd2_pos'] = k2.mean(); r['frac_Kcnd3_only'] = (k3 & ~k2).mean(); r['frac_Kv4_null'] = (~(k1 | k2 | k3)).mean()
    return pd.Series(r)


loop = df.groupby('loop_group', observed=True).apply(summarize)
order = [g for g in LOOP_ORDER if g in loop.index] + [g for g in loop.index if g not in LOOP_ORDER]
loop = loop.loc[order]
loop.to_csv(f'{args.out}/tables/th10x_loop_summary.csv')
# also by nucleus (where annotated) and by subclass
by_nuc = df[df.nucleus.notna()].groupby('nucleus', observed=True).apply(summarize)
by_nuc = by_nuc[by_nuc.n_cells >= args.min_cells]
by_nuc.to_csv(f'{args.out}/tables/th10x_nucleus_summary.csv')
by_sub = df.groupby('subclass', observed=True).apply(summarize)
by_sub = by_sub[by_sub.n_cells >= args.min_cells]
by_sub.to_csv(f'{args.out}/tables/th10x_subclass_summary.csv')

# ------------------------------------------------ figures
show_genes = COMPLEX + [g for g in ['Cacna1g', 'Cacna1h', 'Cacna1i', 'Hcn1', 'Hcn2', 'Hcn4', 'Kcnc1', 'Kcnc2',
                                    'Kcna1', 'Scn8a', 'Pvalb', 'Sst', 'Slc17a6', 'Gad2'] if g in genes]
mat = np.log2(loop[[f'cpm_{g}' for g in show_genes]].values.astype(float) + 1.0)
fig, ax = plt.subplots(figsize=(0.55 * len(show_genes) + 3, 0.45 * len(loop) + 2))
im = ax.imshow(mat, aspect='auto', cmap='magma_r')
ax.set_xticks(range(len(show_genes))); ax.set_xticklabels(show_genes, rotation=60, ha='right', fontsize=8)
ax.set_yticks(range(len(loop))); ax.set_yticklabels([f'{i}  (n={int(n):,})' for i, n in zip(loop.index, loop.n_cells)], fontsize=8)
ax.axvline(len(COMPLEX) - 0.5, color='w', lw=2)
for i in range(len(loop)):
    for j in range(len(show_genes)):
        ax.text(j, i, f'{loop.iloc[i][f"det_{show_genes[j]}"]*100:.0f}', ha='center', va='center', fontsize=5.5,
                color='w' if mat[i, j] > mat.max() * 0.6 else 'k')
plt.colorbar(im, ax=ax, fraction=0.02, label='log2(mean CPM + 1)')
ax.set_title('Kv4 complex (left) and burst machinery (right) by thalamocortical-loop element\n(cell numbers = % of cells detecting the gene)', fontsize=10)
fig.savefig(f'{args.out}/figures/th10x_kv4_loop_heatmap.png', dpi=160, bbox_inches='tight'); plt.close(fig)

# paralog share + partner composition
fig, axes = plt.subplots(1, 3, figsize=(15, 0.4 * len(loop) + 2.5), sharey=True)
y = np.arange(len(loop))
cols = {'Kcnd1': '#4FA8D8', 'Kcnd2': '#E8822D', 'Kcnd3': '#B9A3D9'}
left = np.zeros(len(loop))
for g in KV4:
    v = loop[f'share_{g}'].values
    axes[0].barh(y, v, left=left, color=cols[g], label=g, edgecolor='k', lw=0.3); left += v
axes[0].set_yticks(y); axes[0].set_yticklabels(loop.index, fontsize=8); axes[0].invert_yaxis()
axes[0].set_xlabel('share of Kv4 (Kcnd) CPM'); axes[0].set_title('Kv4 paralog share'); axes[0].legend(fontsize=8, frameon=False)
kc = {'Kcnip1': '#7B4FA3', 'Kcnip2': '#4FA8D8', 'Kcnip3': '#4CAF50', 'Kcnip4': '#E8822D'}
left = np.zeros(len(loop))
for g in KCHIP:
    v = loop[f'KChIP_comp_{g}'].values
    axes[1].barh(y, v, left=left, color=kc[g], label=g, edgecolor='k', lw=0.3); left += v
axes[1].set_xlabel('share of KChIP CPM'); axes[1].set_title('KChIP partner composition'); axes[1].legend(fontsize=8, frameon=False)
left = np.zeros(len(loop))
for g, c in zip(DPP, ['#7B4FA3', '#E8822D']):
    v = loop[f'DPP_comp_{g}'].values
    axes[2].barh(y, v, left=left, color=c, label=g, edgecolor='k', lw=0.3); left += v
axes[2].set_xlabel('share of DPP CPM'); axes[2].set_title('DPP partner composition'); axes[2].legend(fontsize=8, frameon=False)
fig.suptitle('ABC Atlas 10x thalamus: Kv4 complex by loop element', y=1.02)
fig.savefig(f'{args.out}/figures/th10x_kv4_paralog_share.png', dpi=160, bbox_inches='tight'); plt.close(fig)

# per-cell Kcnd2 vs Kcnd3 for the three loop elements that matter most for spike-wave
fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True, sharey=True)
for ax, g in zip(axes, ['TRN (RT)', 'first-order sensory relay', 'higher-order relay']):
    sub = df[df.loop_group == g]
    if len(sub) == 0:
        ax.set_title(f'{g}: no cells'); continue
    ax.scatter(sub['Kcnd3'] + np.random.uniform(-0.08, 0.08, len(sub)), sub['Kcnd2'] + np.random.uniform(-0.08, 0.08, len(sub)),
               s=3, alpha=0.3, color='#3A3A3A')
    ax.set_title(f'{g} (n={len(sub):,})\nKcnd2+ {100*(sub.Kcnd2>0).mean():.0f}%, Kcnd3-only {100*((sub.Kcnd3>0)&(sub.Kcnd2==0)).mean():.0f}%, null {100*((sub[KV4]>0).sum(axis=1)==0).mean():.0f}%', fontsize=9)
    ax.set_xlabel('Kcnd3 log2(CPM+1)')
axes[0].set_ylabel('Kcnd2 log2(CPM+1)')
fig.savefig(f'{args.out}/figures/th10x_percell_kv4_scatter.png', dpi=160, bbox_inches='tight'); plt.close(fig)

pd.set_option('display.width', 220)
print(loop[['n_cells'] + [f'share_{g}' for g in KV4] + ['frac_Kcnd2_pos', 'frac_Kcnd3_only', 'frac_Kv4_null',
                                                     'KChIP_per_Kv4', 'DPP_per_Kv4']].round(3).to_string())
