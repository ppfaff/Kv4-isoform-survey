"""Hippocampal per-subclass variant composition — Kv4 / KChIP / DPP — in the accepted v17 scheme."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip'
DS = sys.argv[1] if len(sys.argv) > 1 else 'hip'
d, ORDER, REGION, SOURCE, CLS = load_dataset(DS)
OUT = f'{OUT}/{DS}'
import os; os.makedirs(OUT, exist_ok=True)
print('cells:', d.sample_name.nunique(), ' transcripts kept:', d.transcript.nunique())

MIN_CELLS = 15 if DS == 'hip' else 5   # subclasses with fewer cells are pooled out of the figure (kept in the CSV)
ncells = d.groupby('subclass')['sample_name'].nunique()
subs = [s for s in ORDER if ncells.get(s, 0) >= MIN_CELLS]
print('subclasses shown:', len(subs), ' dropped (<%d cells):' % MIN_CELLS, [s for s in ORDER if 0 < ncells.get(s, 0) < MIN_CELLS])

# pooled-read fractions per subclass x panel x variant (same definition as the cortex v17 figure)
agg = d.groupby(['subclass', 'cls', 'panel', 'vclass'])['reads'].sum().reset_index()
agg['panel_total'] = agg.groupby(['subclass', 'panel'])['reads'].transform('sum')
agg['frac'] = agg['reads'] / agg['panel_total'].replace(0, np.nan)
agg['n_cells'] = agg['subclass'].map(ncells)
agg.to_csv(f'{OUT}/{DS}_variant_fracs_by_type.csv', index=False)

# also the per-gene fractions (variant / gene) which is what the stats use
gagg = d.groupby(['subclass', 'cls', 'gene', 'vclass'])['reads'].sum().reset_index()
gagg['gene_total'] = gagg.groupby(['subclass', 'gene'])['reads'].transform('sum')
gagg['frac_of_gene'] = gagg['reads'] / gagg['gene_total'].replace(0, np.nan)
gagg['n_cells'] = gagg['subclass'].map(ncells)
gagg.to_csv(f'{OUT}/{DS}_variant_frac_of_gene_by_type.csv', index=False)

PANELS = ['Kv4', 'KChIP', 'DPP']
MIN_READS = 20
BAR_H = 0.75
GAP = 1.2
cls_list = [CLS(s) for s in subs]
ypos, y = [], 0.0
for i, s in enumerate(subs):
    if i > 0 and cls_list[i] != cls_list[i - 1]:
        y += GAP
    ypos.append(y); y += 1.0
ymax = max(ypos)
ypos = {s: ymax - yp for s, yp in zip(subs, ypos)}

fig, axes = plt.subplots(1, 4, figsize=(19, (0.32 if DS == 'hip' else 0.26) * (ymax + 2) + 3.6), sharey=True,
                          gridspec_kw={'width_ratios': [3, 3, 3, 1.1]})
bar_axes, ncell_ax = axes[:3], axes[3]
for ax, panel in zip(bar_axes, PANELS):
    sub = agg[agg['panel'] == panel]
    for s in subs:
        row = sub[sub['subclass'] == s]
        total = row['reads'].sum()
        yy = ypos[s]
        if total < MIN_READS or row.empty:
            ax.barh(yy, 1.0, color='#EEEEEE', height=BAR_H, zorder=1); continue
        left = 0.0
        for _, r in row[row['frac'] > 0].sort_values('frac', ascending=False).iterrows():
            base, kind, pcolor = VSPEC[r['vclass']]
            ax.barh(yy, r['frac'], left=left, color=base, edgecolor=pcolor, hatch=HATCH[kind],
                    height=BAR_H, linewidth=0.3, zorder=2)
            left += r['frac']
    ax.set_xlim(0, 1); ax.set_ylim(-0.75, ymax + 0.75)
    ax.set_title(panel, fontsize=13, fontweight='bold')
    ax.set_yticks([ypos[s] for s in subs]); ax.set_xticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
axes[0].set_yticklabels(subs, fontsize=9 if DS == 'hip' else 7)
for lbl, s in zip(axes[0].get_yticklabels(), subs):
    lbl.set_color(CLS_COLOR[CLS(s)])

nmax = int(ncells.reindex(subs).max())
for s in subs:
    n = int(ncells[s]); yy = ypos[s]
    ncell_ax.barh(yy, n, color='#3A6EA5', height=BAR_H, zorder=2)
    ncell_ax.text(n + nmax * 0.04, yy, f'{n:,}', va='center', ha='left', fontsize=8)
ncell_ax.set_xlim(0, nmax * 1.4); ncell_ax.set_ylim(-0.75, ymax + 0.75)
ncell_ax.set_title('n cells', fontsize=11, fontweight='bold'); ncell_ax.set_xticks([])
for sp in ncell_ax.spines.values(): sp.set_visible(False)

def legend_patches(vlist):
    pats = []
    for v in vlist:
        if v == BLANK:
            pats.append(mpatches.Patch(facecolor='none', edgecolor='none', label=' ')); continue
        base, kind, pcolor = VSPEC[v]
        pats.append(mpatches.Patch(facecolor=base, edgecolor=pcolor, hatch=HATCH[kind], label=v, linewidth=0.3))
    return pats
for ax, panel in zip(bar_axes, PANELS):
    groups = LEGEND_GROUPS[panel]
    xs = [0.5] if len(groups) == 1 else [0.28, 0.76]
    legs = [ax.legend(handles=legend_patches(v), loc='upper center', bbox_to_anchor=(x, -0.01), ncol=1,
                      fontsize=8, frameon=False, handlelength=3.0, handleheight=2.4) for x, v in zip(xs, groups)]
    for leg in legs[:-1]: ax.add_artist(leg)

fig.suptitle(f'{REGION}: per-type variant composition — Kv4 / KChIP / DPP\n'
             f'{SOURCE}, {d.sample_name.nunique():,} cells; pooled-read fractions, most-common segment first; '
             f'subclasses with ≥{MIN_CELLS} cells; grey bar = <{MIN_READS} reads on the family',
             fontsize=11)
fig.subplots_adjust(left=0.12, right=0.98, top=0.90, bottom=0.36, wspace=0.08)
fig.savefig(f'{OUT}/{DS}_variant_scheme_by_type.png', dpi=200)
print('done')
