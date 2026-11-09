#!/usr/bin/env python3
"""Figure for the dentate-gyrus activity program in V404M hippocampus (bulk RNA-seq).

Panel A: log2 fold change (mutant/WT) with standard error for the program genes, four tissues, from the
         original per-tissue DESeq2 (tables/v404m_de_<tissue>.csv); filled = padj<0.05.
Panel B: same genes in hippocampus, unadjusted vs composition-adjusted (tables/v404m_de_adjusted_hippocampus.csv).
Panel C: per-animal hippocampal signature index (tables/v404m_admixture_per_sample.csv): Penk/Drd1/Npy/Ppp1r1b
         driven, MSN genes flat.
Panel D (if tables/v404m_dg_program_per_sample.csv exists, written on Lambda by --extract): per-animal CPM for
         Penk, Pdyn, Drd1, Npy, Kcnip3.

Usage:  python3 scripts/v404m_dg_program.py [TABLES_DIR] [OUT_DIR]            # figure
        python3 scripts/v404m_dg_program.py TABLES OUT --extract                # on Lambda: per-sample CPM table
"""
import os, sys
import numpy as np, pandas as pd

TABLES = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith('--') else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith('--') else '.'
PROGRAM = ['Penk', 'Pdyn', 'Drd1', 'Npy', 'Ecel1', 'Igfbpl1', 'Rprm', 'Cck', 'Rasd1', 'Kcnip3', 'Fosb', 'Fos', 'Egr1', 'Npas4', 'Bdnf']
MSN = ['Gpr88', 'Adora2a', 'Drd2', 'Pde10a']

if '--extract' in sys.argv:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from v404m_panels import load_counts, load_meta
    meta = load_meta(TABLES); C, src = load_counts(TABLES)
    common = [s for s in C.columns if s in meta.index]
    cpm = C[common] / C[common].sum() * 1e6
    genes = [g for g in PROGRAM + MSN + ['Snap25', 'Prox1', 'C1ql2', 'Calb1'] if g in cpm.index]
    out = meta.loc[common].join(cpm.loc[genes, common].T)
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    out.to_csv(f'{OUT}/tables/v404m_dg_program_per_sample.csv')
    print(f'wrote {OUT}/tables/v404m_dg_program_per_sample.csv ({src})')
    sys.exit(0)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
os.makedirs(f'{OUT}/figures', exist_ok=True)
TIS = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']
tcol = {'hippocampus': '#D6372E', 'cortex': '#3A3A3A', 'thalamus': '#4FA8D8', 'cerebellum': '#BBBBBB'}
de = {t: pd.read_csv(f'{TABLES}/v404m_de_{t}.csv').set_index('gene') for t in TIS}
adj = pd.read_csv(f'{TABLES}/v404m_de_adjusted_hippocampus.csv', index_col=0)
per = pd.read_csv(f'{TABLES}/v404m_admixture_per_sample.csv', index_col=0)
has_d = os.path.exists(f'{TABLES}/v404m_dg_program_per_sample.csv')

fig = plt.figure(figsize=(16, 9.5 if has_d else 5.2))
gs = fig.add_gridspec(2 if has_d else 1, 3, width_ratios=[1.4, 1.0, 0.7], hspace=0.45, wspace=0.35)

# A: forest across tissues
ax = fig.add_subplot(gs[0, 0]); y = np.arange(len(PROGRAM))
for k, t in enumerate(TIS):
    d = de[t].reindex(PROGRAM)
    off = (k - 1.5) * 0.18
    sig = d.padj < 0.05
    ax.errorbar(d.log2FoldChange, y + off, xerr=d.lfcSE, fmt='none', ecolor=tcol[t], elinewidth=1, alpha=0.8)
    ax.scatter(d.log2FoldChange[sig], y[sig.values] + off, color=tcol[t], s=34, zorder=3, label=t)
    ax.scatter(d.log2FoldChange[~sig], y[~sig.values] + off, facecolor='w', edgecolor=tcol[t], s=34, zorder=3)
ax.axvline(0, color='k', lw=0.7); ax.set_yticks(y); ax.set_yticklabels(PROGRAM); ax.invert_yaxis()
ax.set_xlabel('log2 fold change, V404M/+ vs +/+ (± SE)'); ax.legend(fontsize=8, frameon=False, loc='lower right')
ax.set_title('A  Dentate activity program is hippocampus-specific\n(filled = padj<0.05, DESeq2 ~sex+genotype)', loc='left', fontsize=10)

# B: hippocampus unadjusted vs adjusted
ax = fig.add_subplot(gs[0, 1]); a = adj.reindex(PROGRAM)
ax.barh(y - 0.2, a.log2FoldChange_without, 0.38, color='#999999', label='~ sex + genotype')
ax.barh(y + 0.2, a.log2FoldChange_with, 0.38, color='#D6372E', label='+ choroid-plexus & cortex indices')
for i, g in enumerate(PROGRAM):
    for off, p in [(-0.2, a.padj_without.iloc[i]), (0.2, a.padj_with.iloc[i])]:
        if pd.notna(p) and p < 0.05:
            v = a.log2FoldChange_without.iloc[i] if off < 0 else a.log2FoldChange_with.iloc[i]
            ax.text(v + (0.03 if v >= 0 else -0.03), i + off, '*', ha='left' if v >= 0 else 'right', va='center', fontsize=9)
ax.axvline(0, color='k', lw=0.7); ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis()
ax.set_xlabel('log2 fold change (hippocampus)'); ax.legend(fontsize=7.5, frameon=False, loc='lower right')
ax.set_title('B  Robust to tissue-composition covariates\n(* padj<0.05)', loc='left', fontsize=10)

# C: per-animal signature index
ax = fig.add_subplot(gs[0, 2]); hp = per[per.tissue == 'hippocampus']; rng = np.random.default_rng(1)
for i, (g, mk, lab) in enumerate([('mutant', 'o', 'V404M/+'), ('wildtype', 's', '+/+')]):
    ss = hp[hp.genotype == g]
    ax.scatter(np.full(len(ss), i) + rng.uniform(-0.08, 0.08, len(ss)), ss['idx_striatum (MSN)'], marker=mk, s=60, color='#666666', edgecolor='k')
    ax.hlines(ss['idx_striatum (MSN)'].mean(), i - 0.25, i + 0.25, color='#D6372E', lw=2)
ax.set_xticks([0, 1]); ax.set_xticklabels(['V404M/+', '+/+']); ax.set_xlim(-0.6, 1.6)
ax.set_ylabel('Penk/Drd1/Ppp1r1b/Npy-containing index (log2)')
ax.set_title('C  Per animal, hippocampus\n(exact p = 0.001 for Penk, Npy, Drd1)', loc='left', fontsize=10)

if has_d:
    D = pd.read_csv(f'{TABLES}/v404m_dg_program_per_sample.csv', index_col=0); D = D[D.tissue == 'hippocampus']
    show = [g for g in ['Penk', 'Pdyn', 'Drd1', 'Npy', 'Ecel1', 'Kcnip3'] if g in D.columns]
    for j, g in enumerate(show):
        ax = fig.add_subplot(gs[1, :]) if j == 0 else ax
        if j == 0:
            ax.set_title('D  Per-animal CPM, hippocampus (circle = V404M/+, square = +/+)', loc='left', fontsize=10)
        for i, (geno, mk) in enumerate([('mutant', 'o'), ('wildtype', 's')]):
            v = D[D.genotype.str.startswith(geno[:3])][g]
            x = j * 3 + i + rng.uniform(-0.12, 0.12, len(v))
            ax.scatter(x, v / D[g].mean(), marker=mk, s=45, color='#D6372E' if i == 0 else '#999999', edgecolor='k')
            ax.hlines((v / D[g].mean()).mean(), j * 3 + i - 0.3, j * 3 + i + 0.3, color='k', lw=1.5)
    ax.set_xticks([j * 3 + 0.5 for j in range(len(show))]); ax.set_xticklabels(show); ax.set_ylabel('CPM / gene mean'); ax.axhline(1, color='k', lw=0.5, ls=':')

fig.suptitle('V404M/+ hippocampus: a seizure-type transcriptional program of dentate granule cells (bulk RNA-seq, 3-month littermates)', y=1.01)
fig.savefig(f'{OUT}/figures/v404m_dg_program.png', dpi=160, bbox_inches='tight')
print('wrote figures/v404m_dg_program.png', '(with per-animal panel)' if has_d else '(run --extract on Lambda for the per-animal panel)')
