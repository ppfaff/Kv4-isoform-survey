#!/usr/bin/env python3
"""Project cluster-level Kv4-complex expression onto the THALMANAC thalamus MERFISH subset.

The ABC Atlas MERFISH panel (500 genes) was built for cell typing and (as of this writing) does not
directly measure most of Kcnd1-3 / Kcnip1-4 / Dpp6 / Dpp10. We therefore paint each MERFISH cell with the
mean expression of its *cluster* taken from the 10x data (th10x_cluster_means.csv from summarize_th_10x.py)
- the same logic as the Allen 'imputed' MERFISH release but restricted to our panel, deterministic, and
without an 8,460-gene download. Whatever panel genes ARE measured directly are reported alongside so the
painting can be sanity-checked (measured vs painted correlation per nucleus).

Uses THALMANAC's ThalamusWrapper (pip install git+https://github.com/AllenNeuralDynamics/abc-merfish-analysis)
for the TH+ZI spatial subset and CCFv3 nucleus labels.

Usage:
  python3 merfish_nucleus_map.py --cache /path/to/abc_atlas_cache --cluster-means tables/th10x_cluster_means.csv --out .

Writes:
  OUT/tables/merfish_th_cells_painted.parquet   per MERFISH cell: section, xyz, nucleus, cluster, painted panel
  OUT/tables/merfish_nucleus_kv4.csv            per CCFv3 nucleus: composition-weighted Kv4 complex
  OUT/figures/merfish_th_sections_Kcnd2.png     painted Kcnd2 and Kcnd3 on the standard THALMANAC sections
  OUT/figures/merfish_nucleus_kv4_bars.png      per-nucleus Kv4 paralog share and KChIP/DPP composition
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from th_common import KV4, KCHIP, DPP, COMPLEX, PANEL, loop_group_of_nucleus, write_cells  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--cache', required=True)
ap.add_argument('--cluster-means', default='tables/th10x_cluster_means.csv')
ap.add_argument('--out', default='.')
ap.add_argument('--release', default='20250531')
args = ap.parse_args()
os.makedirs(f'{args.out}/tables', exist_ok=True); os.makedirs(f'{args.out}/figures', exist_ok=True)

from abc_merfish_analysis.abc_load_thalamus import ThalamusWrapper  # noqa: E402
tw = ThalamusWrapper(directory=args.cache, version=args.release)

# ---- thalamus MERFISH cells with CCF nucleus labels (THALMANAC standard subset, neurons of TH/ZI classes)
obs = tw.get_combined_metadata(drop_unused=True)
obs = tw.filter_by_thalamus_coords(obs, buffer=0, include_white_matter=True)
obs = tw.filter_by_class_thalamus(obs, include=tw.TH_ZI_CLASSES)   # TH Glut, MH-LH Glut, HY GABA (RT-ZI, ZI)
nuc_col = 'parcellation_substructure' if 'parcellation_substructure' in obs.columns else \
    [c for c in obs.columns if 'substructure' in c][0]
print(f'{len(obs):,} thalamic neurons in MERFISH subset; nucleus column = {nuc_col}')

# ---- paint with 10x cluster means
cm = pd.read_csv(args.cluster_means).set_index('cluster')
paint_genes = [g for g in PANEL if f'cpm_{g}' in cm.columns]
painted = obs[['cluster']].join(cm[[f'cpm_{g}' for g in paint_genes]], on='cluster')
painted.columns = ['cluster'] + paint_genes
n_unpainted = painted[paint_genes[0]].isna().sum()
print(f'{n_unpainted:,} cells belong to clusters without a 10x mean (left NaN)')
cells = obs[[c for c in ['brain_section_label', 'x_section', 'y_section', 'z_section', 'x_reconstructed',
                        'y_reconstructed', 'z_reconstructed', nuc_col, 'parcellation_structure',
                        'class', 'subclass', 'supertype', 'cluster'] if c in obs.columns]].copy()
cells = cells.join(painted.drop(columns='cluster'))
cells['loop_group'] = cells[nuc_col].map(loop_group_of_nucleus)

# ---- genes measured directly in the MERFISH panel: record them for validation
gene_md = tw.get_gene_metadata()
measured = [g for g in PANEL if g in set(gene_md['gene_symbol'])]
print('panel genes measured directly by MERFISH:', measured or 'none')
if measured:
    ad = tw.load_adata(transform='log2cpm', from_metadata=obs)
    for g in measured:
        cells[f'measured_{g}'] = np.asarray(ad[:, g].X).ravel()
write_cells(cells, f'{args.out}/tables/merfish_th_cells_painted')

# ---- per-nucleus summary = composition-weighted cluster means
rows = []
for nuc, sub in cells.groupby(nuc_col, observed=True):
    if len(sub) < 30 or nuc in ('unassigned', 'TH-unassigned'):
        continue
    r = {'nucleus': nuc, 'loop_group': loop_group_of_nucleus(nuc), 'n_cells': len(sub),
         'n_clusters': sub['cluster'].nunique(), 'top_cluster': sub['cluster'].value_counts().idxmax(),
         'top_cluster_frac': sub['cluster'].value_counts(normalize=True).iloc[0]}
    for g in paint_genes:
        r[f'cpm_{g}'] = sub[g].mean()
    tot = sum(r[f'cpm_{g}'] for g in KV4)
    for g in KV4:
        r[f'share_{g}'] = r[f'cpm_{g}'] / tot if tot else np.nan
    ktot = sum(r[f'cpm_{g}'] for g in KCHIP); dtot = sum(r[f'cpm_{g}'] for g in DPP)
    for g in KCHIP:
        r[f'KChIP_comp_{g}'] = r[f'cpm_{g}'] / ktot if ktot else np.nan
    for g in DPP:
        r[f'DPP_comp_{g}'] = r[f'cpm_{g}'] / dtot if dtot else np.nan
    for g in measured:
        r[f'measured_{g}'] = sub[f'measured_{g}'].mean()
    rows.append(r)
nuc = pd.DataFrame(rows).sort_values(['loop_group', 'nucleus'])
nuc.to_csv(f'{args.out}/tables/merfish_nucleus_kv4.csv', index=False)
if measured:
    for g in measured:
        if f'cpm_{g}' in nuc:
            print(f'measured vs painted {g}: Spearman r = {nuc[[f"measured_{g}", f"cpm_{g}"]].corr(method="spearman").iloc[0,1]:.2f}')

# ---- sections painted with Kcnd2 and Kcnd3
secs = tw.TH_EXAMPLE_SECTION_LABELS
fig, axes = plt.subplots(2, len(secs), figsize=(5 * len(secs), 9))
for j, s in enumerate(secs):
    sub = cells[cells['brain_section_label'] == s]
    for i, g in enumerate(['Kcnd2', 'Kcnd3']):
        ax = axes[i, j]
        v = np.log2(sub[g] + 1)
        sc = ax.scatter(sub['x_section'], sub['y_section'], c=v, s=2, cmap='magma_r', vmin=0, vmax=np.nanpercentile(np.log2(cells[g] + 1), 99))
        ax.set_aspect('equal'); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f'{s}  painted {g}', fontsize=9)
        if j == len(secs) - 1:
            plt.colorbar(sc, ax=ax, fraction=0.03, label='log2(CPM+1), cluster mean')
fig.suptitle('THALMANAC thalamus sections painted with 10x cluster-level Kcnd2 / Kcnd3', y=0.98)
fig.savefig(f'{args.out}/figures/merfish_th_sections_Kcnd2.png', dpi=150, bbox_inches='tight'); plt.close(fig)

# ---- per-nucleus bars
fig, axes = plt.subplots(1, 3, figsize=(15, 0.32 * len(nuc) + 2), sharey=True)
y = np.arange(len(nuc))
for ax, fam, cols, lab in [(axes[0], KV4, ['#4FA8D8', '#E8822D', '#B9A3D9'], 'share'),
                           (axes[1], KCHIP, ['#7B4FA3', '#4FA8D8', '#4CAF50', '#E8822D'], 'KChIP_comp'),
                           (axes[2], DPP, ['#7B4FA3', '#E8822D'], 'DPP_comp')]:
    left = np.zeros(len(nuc))
    for g, c in zip(fam, cols):
        v = nuc[f'{lab}_{g}'].values
        ax.barh(y, v, left=left, color=c, label=g, edgecolor='k', lw=0.3); left += v
    ax.legend(fontsize=8, frameon=False); ax.set_xlim(0, 1)
axes[0].set_yticks(y); axes[0].set_yticklabels([f'{a}  [{b}]  n={int(c):,}' for a, b, c in zip(nuc.nucleus, nuc.loop_group, nuc.n_cells)], fontsize=7)
axes[0].invert_yaxis()
axes[0].set_title('Kv4 paralog share'); axes[1].set_title('KChIP composition'); axes[2].set_title('DPP composition')
fig.suptitle('Per CCFv3 nucleus (MERFISH cells, 10x cluster means)', y=1.01)
fig.savefig(f'{args.out}/figures/merfish_nucleus_kv4_bars.png', dpi=160, bbox_inches='tight'); plt.close(fig)
print(nuc[['nucleus', 'loop_group', 'n_cells', 'n_clusters'] + [f'share_{g}' for g in KV4]].round(3).to_string(index=False))
