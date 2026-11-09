#!/usr/bin/env python3
"""Thalamus view of the V404M bulk RNA-seq Kv4-complex tables.

Re-cuts the per-tissue tables written by scripts/v404m_kv4.py (tables/v404m_kv4_genes.csv,
tables/v404m_variant_fracs.csv) to answer, for the thalamus specifically:

  1. How much Kv4.2 does thalamus express relative to Kv4.1/Kv4.3, and relative to the other tissues?
  2. Which auxiliary-subunit N-terminal variants dominate in thalamus, and how does that composition
     differ from cortex/hippocampus (the tissues already surveyed at single-cell resolution)?
  3. Does the heterozygous V404M knock-in change any of this in thalamus (compensation test)?

Membrane-association classes of the KChIP N-termini follow the lab's 2015 review, extended with
KChIP1e as a transmembrane variant (AlphaFold). Variants not in the published scheme
(Kcnip1-new_short, Kcnip3-new-cys_rich, 'other/fragment') are left as 'unclassified'.

Usage: python3 scripts/v404m_thalamus_view.py [TABLES_DIR] [OUT_DIR]
Writes OUT/tables/thalamus_bulk_kv4_summary.csv, OUT/tables/thalamus_bulk_kchip_classes.csv,
       OUT/figures/thalamus_bulk_kv4_complex.png
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import VSPEC, HATCH  # noqa: E402

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)

TISSUES = ['cortex', 'hippocampus', 'thalamus', 'cerebellum']
KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']
AUX = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']
KCHIP_CLASS = {
    'cytoplasmic': ['KChIP2c', 'KChIP4b', 'KChIP4c'],
    'N-myristoyl': ['KChIP1a', 'KChIP1b'],
    'S-palmitoyl': ['KChIP1c', 'KChIP4d', 'KChIP2a', 'KChIP3a', 'KChIP2g', 'KChIP2b', 'KChIP2aL', 'KChIP4bL'],
    'transmembrane': ['KChIP2x', 'KChIP3x', 'KChIP4a', 'KChIP4e', 'KChIP1e'],
}
CLASS_OF = {v: k for k, vs in KCHIP_CLASS.items() for v in vs}
CLASS_COLOR = {'cytoplasmic': '#CFCFCF', 'N-myristoyl': '#4FA8D8', 'S-palmitoyl': '#E8822D',
               'transmembrane': '#7B4FA3', 'unclassified': '#FFFFFF'}

genes = pd.read_csv(f'{TABLES}/v404m_kv4_genes.csv')
fr = pd.read_csv(f'{TABLES}/v404m_variant_fracs.csv')
genes = genes[genes.tissue.isin(TISSUES)]
fr = fr[fr.tissue.isin(TISSUES)]

# ---------- 1. gene levels and Kv4 paralog share ----------
cpm_wt = genes.pivot(index='gene', columns='tissue', values='mean_logCPM_wt').rpow(2.0)[TISSUES]
cpm_mut = genes.pivot(index='gene', columns='tissue', values='mean_logCPM_mut').rpow(2.0)[TISSUES]
kv4_share = cpm_wt.loc[KV4] / cpm_wt.loc[KV4].sum()
kchip_total = cpm_wt.loc[['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']].sum()
dpp_total = cpm_wt.loc[['Dpp6', 'Dpp10']].sum()
kv4_total = cpm_wt.loc[KV4].sum()

summary = pd.DataFrame({
    'tissue': TISSUES,
    'Kcnd2_CPM_wt': cpm_wt.loc['Kcnd2'].values,
    'Kcnd2_share_of_Kv4': kv4_share.loc['Kcnd2'].values,
    'Kcnd3_share_of_Kv4': kv4_share.loc['Kcnd3'].values,
    'Kcnd1_share_of_Kv4': kv4_share.loc['Kcnd1'].values,
    'KChIP_total_CPM': kchip_total.values,
    'DPP_total_CPM': dpp_total.values,
    'KChIP_per_Kv4': (kchip_total / kv4_total).values,
    'DPP_per_Kv4': (dpp_total / kv4_total).values,
    'Dpp6_share_of_DPP': (cpm_wt.loc['Dpp6'] / dpp_total).values,
})
for g in KV4 + AUX:
    sub = genes[genes.gene == g].set_index('tissue').reindex(TISSUES)
    summary[f'{g}_log2FC_mut'] = sub['log2FC'].values
    summary[f'{g}_p_mut'] = sub['p'].values

# ---------- 2. KChIP membrane-association classes per tissue (WT), CPM-weighted across the 4 genes ----------
kchip = fr[fr.gene.str.startswith('Kcnip')].copy()
kchip['class'] = kchip.vclass.map(CLASS_OF).fillna('unclassified')
rows = []
for t in TISSUES:
    sub = kchip[kchip.tissue == t]
    w = sub.gene.map(cpm_wt[t])             # weight each gene's variant fractions by its CPM
    wm = sub.gene.map(cpm_mut[t])
    kg = sub.gene.unique()                  # denominator: total Kcnip CPM (fractions sum to 1 within a gene)
    for cls in list(KCHIP_CLASS) + ['unclassified']:
        m = sub['class'] == cls
        rows.append(dict(tissue=t, kchip_class=cls,
                         frac_wt=(sub.frac_wt[m] * w[m]).sum() / cpm_wt.loc[kg, t].sum(),
                         frac_mut=(sub.frac_mut[m] * wm[m]).sum() / cpm_mut.loc[kg, t].sum()))
classes = pd.DataFrame(rows)
classes['diff'] = classes.frac_mut - classes.frac_wt
classes.to_csv(f'{OUT}/tables/thalamus_bulk_kchip_classes.csv', index=False)
summary = summary.merge(classes.pivot(index='tissue', columns='kchip_class', values='frac_wt')
                        .add_prefix('KChIP_class_wt_').reset_index(), on='tissue')
summary.to_csv(f'{OUT}/tables/thalamus_bulk_kv4_summary.csv', index=False)

# ---------- 3. figure ----------
fig = plt.figure(figsize=(16, 11.5))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.35], hspace=0.5, wspace=0.4)
cols = {'cortex': '#999999', 'hippocampus': '#666666', 'thalamus': '#D6372E', 'cerebellum': '#BBBBBB'}

# A: WT CPM per gene per tissue
ax = fig.add_subplot(gs[0, 0:2])
x = np.arange(len(KV4 + AUX)); wbar = 0.2
for i, t in enumerate(TISSUES):
    ax.bar(x + (i - 1.5) * wbar, cpm_wt.loc[KV4 + AUX, t], wbar, color=cols[t], label=t,
           edgecolor='k' if t == 'thalamus' else 'none', lw=0.8)
ax.set_xticks(x); ax.set_xticklabels(KV4 + AUX, rotation=45, ha='right')
ax.set_ylabel('wild-type CPM (mean of log2 CPM)'); ax.set_yscale('log')
ax.set_title('A  Kv4-complex gene levels by tissue', loc='left'); ax.legend(fontsize=8, frameon=False)

# B: Kv4 paralog share
ax = fig.add_subplot(gs[0, 2])
bottom = np.zeros(len(TISSUES))
for g in KV4:
    vals = kv4_share.loc[g, TISSUES].values
    ax.bar(TISSUES, vals, bottom=bottom, color=VSPEC['Kv4.2'][0] if g == 'Kcnd2' else
           (VSPEC['Kv4.1'][0] if g == 'Kcnd1' else VSPEC['Kv4.3S'][0]), label=g.replace('Kcnd', 'Kv4.'), edgecolor='k', lw=0.4)
    for xi, (b, v) in enumerate(zip(bottom, vals)):
        if v > 0.06:
            ax.text(xi, b + v / 2, f'{v:.2f}', ha='center', va='center', fontsize=8)
    bottom += vals
ax.set_ylabel('share of total Kcnd CPM'); ax.set_ylim(0, 1)
ax.set_title('B  Kv4 paralog share (wild type)', loc='left'); ax.legend(fontsize=8, frameon=False, bbox_to_anchor=(1.0, 1.0), loc='upper left')
ax.tick_params(axis='x', rotation=30)

# C: KChIP membrane-association classes
ax = fig.add_subplot(gs[0, 3])
bottom = np.zeros(len(TISSUES))
for cls in list(KCHIP_CLASS) + ['unclassified']:
    vals = classes[classes.kchip_class == cls].set_index('tissue').loc[TISSUES, 'frac_wt'].values
    ax.bar(TISSUES, vals, bottom=bottom, color=CLASS_COLOR[cls], edgecolor='k', lw=0.4, label=cls)
    for xi, (b, v) in enumerate(zip(bottom, vals)):
        if v > 0.07:
            ax.text(xi, b + v / 2, f'{v:.2f}', ha='center', va='center', fontsize=8,
                    color='w' if cls == 'transmembrane' else 'k')
    bottom += vals
ax.set_ylabel('fraction of KChIP transcripts (CPM-weighted)'); ax.set_ylim(0, 1)
ax.set_title('C  KChIP N-terminal class (wild type)', loc='left'); ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc='upper left')
ax.tick_params(axis='x', rotation=30)

# D: per-gene variant composition, WT, tissue x gene
ax = fig.add_subplot(gs[1, :3])
xpos = 0; ticks = []; tlabels = []
legend_handles = {}
for g in AUX + ['Kcnd3']:
    for t in TISSUES:
        sub = fr[(fr.gene == g) & (fr.tissue == t)].sort_values('frac_wt', ascending=False)
        b = 0
        for _, r in sub.iterrows():
            if r.frac_wt < 0.005:
                continue
            face, hatch, hcol = VSPEC.get(r.vclass, ('#DDDDDD', 'solid', '#DDDDDD'))
            ax.bar(xpos, r.frac_wt, 0.8, bottom=b, color=face, hatch=HATCH[hatch], edgecolor=hcol if HATCH[hatch] else 'k', lw=0.4)
            legend_handles.setdefault(r.vclass, Patch(facecolor=face, hatch=HATCH[hatch], edgecolor=hcol if HATCH[hatch] else 'k'))
            b += r.frac_wt
        ticks.append(xpos); tlabels.append(t[:3] if t != 'cerebellum' else 'cb')
        if t == 'thalamus':
            ax.plot([xpos - 0.45, xpos + 0.45], [1.02, 1.02], color='#D6372E', lw=3, clip_on=False)
        xpos += 1
    ax.text(xpos - 2.5, 1.07, g, ha='center', fontsize=10, fontweight='bold')
    xpos += 0.8
ax.set_xticks(ticks); ax.set_xticklabels(tlabels, fontsize=7)
ax.set_ylim(0, 1); ax.set_ylabel('variant fraction (wild type)')
ax.set_title('D  N-terminal variant composition by tissue (red bar = thalamus)', loc='left')
ax.legend(legend_handles.values(), legend_handles.keys(), fontsize=6.5, frameon=False, ncol=8,
          bbox_to_anchor=(0.5, -0.12), loc='upper center')

# E: thalamus mutant vs WT — variant fraction differences
ax = fig.add_subplot(gs[1, 3])
th = fr[(fr.tissue == 'thalamus') & (fr.frac_wt > 0.02)].copy().sort_values(['gene', 'frac_wt'], ascending=[True, False])
ax.barh(np.arange(len(th)), th['diff'], color=['#D6372E' if p < 0.05 else '#AAAAAA' for p in th.p])
ax.set_yticks(np.arange(len(th))); ax.set_yticklabels(th.vclass, fontsize=7)
ax.axvline(0, color='k', lw=0.6); ax.invert_yaxis()
ax.set_xlabel('mutant − wild-type fraction'); ax.set_title('E  Thalamus: V404M vs WT\n(red = nominal p<0.05; none survive FDR)', loc='left', fontsize=10)
for i, (d, p) in enumerate(zip(th['diff'], th.p)):
    ax.text(d + (0.002 if d >= 0 else -0.002), i, f'p={p:.2f}', va='center', ha='left' if d >= 0 else 'right', fontsize=6)

fig.suptitle('V404M bulk RNA-seq: the Kv4 complex in thalamus vs the other tissues', fontsize=13, y=0.98)
fig.savefig(f'{OUT}/figures/thalamus_bulk_kv4_complex.png', dpi=160, bbox_inches='tight')
print(summary.round(3).T.to_string())
print()
print(classes.pivot(index='kchip_class', columns='tissue', values='frac_wt').round(3).to_string())
