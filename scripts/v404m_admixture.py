#!/usr/bin/env python3
"""Per-library tissue-admixture check for the V404M bulk RNA-seq.

Motivation (2026-11-09): the 'cortex' libraries carry striatal MSN markers (Gpr88, Adora2a, Ppp1r1b, Penk,
Drd1/2) at 20-25x hippocampal levels, and Kcnip2 - at ~5% of Kcnip3 in Tasic 2016 cortical neurons - equals
Kcnip3 in the bulk. This script scores every library for striatal, cortical, hippocampal, thalamic and
cerebellar signatures, tests whether admixture differs by genotype (which would confound the per-tissue DE),
and correlates the striatal index with Kcnip2.

Runs where the per-sample counts matrix lives (tables/v404m_star_counts.csv.gz or v404m_gene_counts.csv.gz).
Usage: python3 scripts/v404m_admixture.py [TABLES_DIR] [OUT_DIR]
Writes OUT/tables/v404m_admixture_per_sample.csv and OUT/figures/v404m_admixture.png
"""
import os, sys
import numpy as np, pandas as pd
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:          # the kv4 env has no matplotlib; tables are the deliverable, figure is optional
    plt = None
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v404m_panels import load_counts, load_meta  # noqa: E402


def mannwhitney_p(a, b):
    """Two-sided Mann-Whitney U, normal approximation with tie correction (no scipy in the kv4 env)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return np.nan
    allv = np.concatenate([a, b]); ranks = pd.Series(allv).rank().values
    u1 = ranks[:n1].sum() - n1 * (n1 + 1) / 2
    mu = n1 * n2 / 2
    _, cnt = np.unique(allv, return_counts=True); n = n1 + n2
    sigma = np.sqrt(n1 * n2 / 12 * ((n + 1) - (cnt ** 3 - cnt).sum() / (n * (n - 1))))
    if sigma == 0:
        return 1.0
    z = max(abs(u1 - mu) - 0.5, 0.0) / sigma
    import math
    return float(min(1.0, 2 * (1 - 0.5 * (1 + math.erf(z / np.sqrt(2))))))

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)

SIG = {
    'striatum (MSN)': ['Gpr88', 'Adora2a', 'Ppp1r1b', 'Penk', 'Drd1', 'Drd2', 'Rgs9', 'Pde10a'],
    'cortex (glut)': ['Slc17a7', 'Satb2', 'Cux2', 'Rorb', 'Fezf2', 'Foxp2', 'Syt6'],
    'hippocampus': ['Prox1', 'C1ql2', 'Dsp', 'Zbtb20', 'Amigo2', 'Pcp4'],
    'thalamus (relay)': ['Tcf7l2', 'Prkcd', 'Gbx2', 'Slc17a6', 'Ntng1'],
    'reticular / ZI': ['Gnb3', 'Pvalb', 'Six3', 'Isl1'],
    'cerebellum': ['Pcp2', 'Gabra6', 'Cbln1', 'Zic1', 'Calb1'],
}
meta = load_meta(TABLES)
C, src = load_counts(TABLES)
common = [s for s in C.columns if s in meta.index]
C, meta = C[common], meta.loc[common]
cpm = C / C.sum() * 1e6
snap = cpm.loc['Snap25']

rows = []
for s in common:
    r = {'sample': s, **meta.loc[s].to_dict(), 'Snap25_cpm': snap[s]}
    for sig, genes in SIG.items():
        g = [x for x in genes if x in cpm.index]
        r[f'idx_{sig}'] = np.log2(cpm.loc[g, s] + 1).mean() - np.log2(snap[s] + 1)   # neuron-normalised log index
    for g in ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Kcnd2', 'Kcnd3']:
        r[f'{g}_cpm'] = cpm.loc[g, s]
    r['Kcnip2_over_Kcnip3'] = cpm.loc['Kcnip2', s] / cpm.loc['Kcnip3', s]
    rows.append(r)
A = pd.DataFrame(rows).set_index('sample')
A.to_csv(f'{OUT}/tables/v404m_admixture_per_sample.csv')

pd.set_option('display.width', 220)
idx_cols = [c for c in A.columns if c.startswith('idx_')]
print(f'counts source: {src}\n')
print('mean admixture index by tissue (log2 signature CPM relative to Snap25):')
print(A.groupby('tissue')[idx_cols].mean().round(2).to_string())
print('\nwithin-tissue SD of the striatal index (dissection variability):')
print(A.groupby('tissue')['idx_striatum (MSN)'].std().round(2).to_string())
print('\nstriatal index by genotype within tissue (Mann-Whitney p):')
for t, sub in A.groupby('tissue'):
    a, b = sub[sub.genotype.str.startswith('mut')]['idx_striatum (MSN)'], sub[~sub.genotype.str.startswith('mut')]['idx_striatum (MSN)']
    p = mannwhitney_p(a, b)
    print(f'  {t:12s} mut {a.mean():.2f}  wt {b.mean():.2f}  p={p:.3f}')
ctx = A[A.tissue == 'cortex']
if len(ctx) > 3:
    rho = ctx['idx_striatum (MSN)'].rank().corr(ctx['Kcnip2_over_Kcnip3'].rank())   # Spearman via ranks (no scipy)
    print(f'\ncortex: Spearman(striatal index, Kcnip2/Kcnip3) = {rho:.2f}  (n={len(ctx)})')

# per-gene genotype ratios for the striatal signature: admixture moves MSN-specific genes (Gpr88, Adora2a,
# Drd1/2) together with the shared ones; an activity-induced granule-cell response moves Penk/Ppp1r1b/Pde10a alone
print('\nstriatal-signature genes, log2(mutant/WT mean CPM) per tissue:')
rows = []
for t, sub in meta.groupby('tissue'):
    mut = sub.index[sub.genotype.str.startswith('mut')]; wt = sub.index[~sub.genotype.str.startswith('mut')]
    for g in SIG['striatum (MSN)'] + ['Kcnip2', 'Fosb', 'Npy', 'Bdnf', 'Arc', 'Fos', 'Egr1']:
        if g in cpm.index:
            rows.append(dict(tissue=t, gene=g, wt_cpm=cpm.loc[g, wt].mean(),
                             log2_mut_over_wt=np.log2((cpm.loc[g, mut].mean() + 0.5) / (cpm.loc[g, wt].mean() + 0.5)),
                             p=mannwhitney_p(cpm.loc[g, mut], cpm.loc[g, wt])))
G = pd.DataFrame(rows)
G.to_csv(f'{OUT}/tables/v404m_admixture_genes.csv', index=False)
print(G.pivot(index='gene', columns='tissue', values='log2_mut_over_wt').round(2).to_string())
print('\n(p-values in tables/v404m_admixture_genes.csv)')

if plt is None:
    print('\nmatplotlib not available: skipped the figure; plot from tables/v404m_admixture_per_sample.csv elsewhere')
    sys.exit(0)
fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
tis = ['cortex', 'hippocampus', 'thalamus', 'cerebellum']
cols = {'cortex': '#999999', 'hippocampus': '#666666', 'thalamus': '#D6372E', 'cerebellum': '#BBBBBB'}
for t in tis:
    sub = A[A.tissue == t]
    axes[0].scatter(sub['idx_striatum (MSN)'], sub['Kcnip2_over_Kcnip3'], color=cols[t], label=t,
                    marker=['o' if g.startswith('mut') else 's' for g in sub.genotype][0] if len(sub) else 'o', s=40, edgecolor='k')
    for _, r in sub.iterrows():
        axes[0].scatter(r['idx_striatum (MSN)'], r['Kcnip2_over_Kcnip3'], color=cols[t], s=40, edgecolor='k',
                        marker='o' if r.genotype.startswith('mut') else 's')
axes[0].set_xlabel('striatal MSN index (log2 vs Snap25)'); axes[0].set_ylabel('Kcnip2 / Kcnip3 CPM')
axes[0].set_title('Kcnip2 tracks striatal admixture?  (circle = mutant, square = WT)'); axes[0].legend(frameon=False, fontsize=8)
M = A.groupby('tissue')[idx_cols].mean().loc[tis]
im = axes[1].imshow(M.values, cmap='viridis', aspect='auto')
axes[1].set_xticks(range(len(idx_cols))); axes[1].set_xticklabels([c[4:] for c in idx_cols], rotation=40, ha='right', fontsize=8)
axes[1].set_yticks(range(len(tis))); axes[1].set_yticklabels(tis)
for i in range(len(tis)):
    for j in range(len(idx_cols)):
        axes[1].text(j, i, f'{M.values[i, j]:.1f}', ha='center', va='center', fontsize=7, color='w')
axes[1].set_title('signature index by tissue'); plt.colorbar(im, ax=axes[1], fraction=0.04)
fig.savefig(f'{OUT}/figures/v404m_admixture.png', dpi=150, bbox_inches='tight')
