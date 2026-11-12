"""Simplified variant-composition figure: hippocampus (Allen SSv4 HPF) vs cortex (Allen SSv4 VISp),
selected classes only, inhibitory and excitatory separated. Bars = mean per-neuron variant fraction
(each cell's variant reads / its family reads, cells with >=20 family reads), i.e. every neuron weighted equally."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip/compare'
CTX_SOURCE = sys.argv[1] if len(sys.argv) > 1 else 'visp'   # 'visp' (Allen SSv4) or 't16' (Tasic 2016)

hip = load_hip()
if CTX_SOURCE == 'visp':
    v = pd.read_pickle('/home/claude/kv4/visp/visp_per_cell.pkl')
    v['vclass'] = v['transcript'].map(final_map); v = v.dropna(subset=['vclass']).copy()
    v = v[~((v['biotype'] == 'retained_intron') & v['vclass'].str.contains('other|fragment', case=False))].copy()
    v['panel'] = v['gene'].map(PANEL_OF); v = v.rename(columns={'subclass_label': 'subclass'})
    ctx = v
    CTX = [('Vip', ['Vip']), ('Sst', ['Sst', 'Sst Chodl']), ('Pvalb', ['Pvalb']),
           ('L2/3 pyramidal', ['L2/3 IT CTX']), ('L5 IT pyramidal', ['L5 IT CTX']), ('L5 ET pyramidal', ['L5 PT CTX']),
           ('L6 CT pyramidal', ['L6 CT CTX']), ('L6 CC pyramidal', ['L6 IT CTX']), ('L6b pyramidal', ['L6b CTX'])]
    ctx_label = 'Cortex (V1) — Allen SSv4'
else:
    ctx = load_ctx()
    T = lambda p: [t for t in ctx['subclass'].unique() if t.startswith(p)]
    CTX = [('Vip', [t for t in T('Vip') if t != 'Vip Sncg']), ('Sst', T('Sst')), ('Pvalb', T('Pvalb')),
           ('L2/3 pyramidal', ['L2/3 Ptgs2', 'L2 Ngb']), ('L5 IT pyramidal', T('L5a')), ('L5 ET pyramidal', T('L5b')),
           ('L6 CT pyramidal', ['L6a Mgp', 'L6a Sla', 'L6a Car12']), ('L6 CC pyramidal', ['L6a Syt17']), ('L6b pyramidal', T('L6b'))]
    ctx_label = 'Cortex (V1) — Tasic 2016'
HIP = [('Vip', ['Vip']), ('Sst', ['Sst', 'Sst Chodl']), ('Pvalb', ['Pvalb']),
       ('DG', ['DG']), ('CA1', ['CA1-ProS']), ('CA3', ['CA3'])]
INHIB = {'Vip', 'Sst', 'Pvalb'}
PANEL_V = {
    'Kv4':  ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short', 'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich', 'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':  ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}
MIN_READS = 20

def mean_fracs(d, sel):
    d = d[d['subclass'].isin(sel)]
    n_cells = d['sample_name'].nunique()
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    out, n_used = {}, {}
    for panel, vs in PANEL_V.items():
        X = V.reindex(columns=vs).fillna(0.0); tot = X.sum(axis=1); ok = tot >= MIN_READS
        out[panel] = (X[ok].div(tot[ok], axis=0)).mean().to_dict() if ok.sum() else {v: 0 for v in vs}
        n_used[panel] = int(ok.sum())
    return out, n_cells, n_used

rows = [('Hippocampal formation — Allen SSv4', hip, HIP), (ctx_label, ctx, CTX)]
comp, table = [], []
for region, d, classes in rows:
    for label, sel in classes:
        fr, n, nu = mean_fracs(d, sel)
        comp.append((region, label, fr, n, nu))
        for panel in PANEL_V:
            for vc, f in fr[panel].items():
                table.append(dict(region=region, cell_class=label, n_cells=n, n_cells_used=nu[panel], panel=panel, variant=vc, mean_fraction=f))
pd.DataFrame(table).to_csv(f'{OUT}/simple_composition_{CTX_SOURCE}.csv', index=False)

# layout: hippocampus block on top, cortex below; inhibitory then excitatory within each; gaps between groups
ypos, y, labels, region_of = {}, 0.0, [], {}
for ri, (region, d, classes) in enumerate(rows):
    if ri > 0: y += 1.6
    prev_inh = None
    for label, sel in classes:
        inh = label in INHIB
        if prev_inh is not None and inh != prev_inh: y += 0.7
        ypos[(region, label)] = y; region_of[(region, label)] = region; y += 1.0; prev_inh = inh
ymax = y
fig, axes = plt.subplots(1, 3, figsize=(19, 0.5 * ymax + 1.2), sharey=True)
BH = 0.72
for ax, panel in zip(axes[:3], PANEL_V):
    for region, label, fr, n, nu in comp:
        yy = ymax - ypos[(region, label)]; left = 0.0
        for vc, f in sorted(fr[panel].items(), key=lambda kv: -kv[1]):
            if f <= 0.002: continue
            base, kind, pcolor = VSPEC[vc]
            ax.barh(yy, f, left=left, color=base, edgecolor=pcolor, hatch=HATCH[kind], height=BH, linewidth=0.3)
            left += f
    ax.set_xlim(0, 1); ax.set_ylim(0.2, ymax + 0.8); ax.set_xticks([])
    ax.set_title(panel, fontsize=15, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(False)
keys = [(r, l) for r, d, cl in rows for l, s in cl]
axes[0].set_yticks([ymax - ypos[k] for k in keys]); axes[0].set_yticklabels([k[1] for k in keys], fontsize=13, fontweight='bold')
for lbl, k in zip(axes[0].get_yticklabels(), keys):
    lbl.set_color(CLS_COLOR['GABA'] if k[1] in INHIB else CLS_COLOR['Glut'])
# region headers
for region, d, classes in rows:
    top = ymax - ypos[(region, classes[0][0])] + 0.75
    axes[0].text(-0.02, top, region, fontsize=12, fontweight='bold', ha='right', va='bottom', color='#333', transform=axes[0].get_yaxis_transform())
# legends
DISPLAY = {'Kcnip1-new_short': 'Kcnip1 short', 'Kcnip3-new-cys_rich': 'Kcnip3 Cys-rich'}
maxfrac = {}
for region, label, fr, n, nu in comp:
    for panel in PANEL_V:
        for vc, f in fr[panel].items(): maxfrac[vc] = max(maxfrac.get(vc, 0), f)
UNSEEN = {vc for vc, m in maxfrac.items() if m < 0.05}
def legend_patches(vlist):
    pats = []
    for vc in vlist:
        if vc == BLANK: pats.append(mpatches.Patch(facecolor='none', edgecolor='none', label=' ')); continue
        base, kind, pcolor = VSPEC[vc]
        pats.append(mpatches.Patch(facecolor=base, edgecolor=pcolor, hatch=HATCH[kind], label=DISPLAY.get(vc, vc), linewidth=0.3))
    return pats
for ax, panel in zip(axes[:3], PANEL_V):
    gl = LEGEND_GROUPS[panel]; xs = [0.5] if len(gl) == 1 else [0.28, 0.76]
    legs = []
    for x, vl in zip(xs, gl):
        vl = [v for v in vl if v in PANEL_V[panel] or v == BLANK]
        leg = ax.legend(handles=legend_patches(vl), loc='upper center', bbox_to_anchor=(x, -0.01),
                        ncol=1, fontsize=13, frameon=False, handlelength=2.6, handleheight=1.6, labelspacing=0.35)
        for txt, vc in zip(leg.get_texts(), vl):
            if vc in UNSEEN: txt.set_color('#D62728')
        legs.append(leg)
    for leg in legs[:-1]: ax.add_artist(leg)
fig.suptitle('Kv4 / KChIP / DPP variant composition — hippocampal formation vs cortex, selected classes\n'
             f'bar = mean per-neuron variant fraction (cells with ≥{MIN_READS} reads on the family); purple = inhibitory, blue = excitatory; '
             'both datasets Allen SMART-seq v4; red legend entries = variant below 5% in every class shown' if CTX_SOURCE == 'visp' else
             'Kv4 / KChIP / DPP variant composition — hippocampal formation (Allen SSv4) vs cortex (Tasic 2016), selected classes\n'
             f'bar = mean per-neuron variant fraction (cells with ≥{MIN_READS} reads on the family); purple = inhibitory, blue = excitatory; red legend entries = variant below 5% in every class shown; note the chemistry difference between datasets',
             fontsize=12)
fig.subplots_adjust(left=0.17, right=0.985, top=0.9, bottom=0.40, wspace=0.08)
fig.savefig(f'{OUT}/simple_composition_{CTX_SOURCE}.png', dpi=200)
print('written', CTX_SOURCE); print('below 5% everywhere:', sorted(UNSEEN))
for region, label, fr, n, nu in comp:
    print(f'{region[:12]:12s} {label:16s} n={n:5d}  used Kv4/KChIP/DPP = {nu["Kv4"]}/{nu["KChIP"]}/{nu["DPP"]}')
