"""Cortex (Tasic 2016, V1) vs hippocampal formation (Allen SSv4, Yao 2021): matched-subclass comparison of
Kv4 / KChIP / DPP variant composition, family expression (CPM) and depleted-cell fractions at matched depth."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip/compare'
import os; os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)

hip = load_hip(); ctx = load_ctx()
# ---- matched groups: name -> (cortex types, hip subclasses)
CTX_T = lambda prefix: [t for t in ctx['subclass'].unique() if t.startswith(prefix)]
GROUPS = {
    'Lamp5':            (['Ndnf Car4', 'Ndnf Cxcl14', 'Igtp'], ['Lamp5']),
    'Lamp5 Lhx6 / Smad3': (['Smad3'], None),   # hip side = Lamp5 Lhx6 clusters, handled below
    'Sncg':             (['Sncg', 'Vip Sncg'], ['Sncg']),
    'Vip':              ([t for t in CTX_T('Vip') if t != 'Vip Sncg'], ['Vip']),
    'Sst':              (CTX_T('Sst'), ['Sst', 'Sst Chodl']),
    'Pvalb':            (CTX_T('Pvalb'), ['Pvalb']),
    'L2/3 IT':          (['L2/3 Ptgs2', 'L2 Ngb'], ['L2/3 IT RHP', 'L2/3 IT PPP', 'L2/3 IT ENTl']),
    'L6 CT / L6a':      (CTX_T('L6a'), ['L6 CT CTX']),
    'L6b':              (CTX_T('L6b'), ['L6b CTX']),
    'all GABA':         (None, None),
    'all Glut (ctx) vs CA/DG/SUB (hip)': (None, None),
}
def group_cells(d, which, sel):
    if sel is None:
        if which == 'all GABA': return d[d.cls == 'GABA']
        if which.startswith('all Glut'): return d[d.cls == 'Glut']
    return d[d['subclass'].isin(sel)]

hip_lhx6 = hip[hip['cluster'].str.contains('Lamp5 Lhx6')]

def cells_of(name, region):
    ct, hs = GROUPS[name]
    if region == 'ctx':
        return group_cells(ctx, name, ct)
    if name == 'Lamp5 Lhx6 / Smad3':
        return hip_lhx6
    return group_cells(hip, name, hs)

# ---- per-cell tables: variant reads, family CPM
def percell(d):
    meta = d[['sample_name', 'subclass', 'cls', 'num_processed']].drop_duplicates().set_index('sample_name')
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0).reindex(meta.index).fillna(0.0)
    G = d[['sample_name', 'gene', 'gene_reads']].drop_duplicates().pivot(index='sample_name', columns='gene', values='gene_reads').reindex(meta.index).fillna(0.0)
    return meta, V, G
FAM = {'Kv4': ['Kcnd1', 'Kcnd2', 'Kcnd3'], 'KChIP': ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'], 'DPP': ['Dpp6', 'Dpp10']}
PANEL_V = {
    'Kv4':  ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short', 'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich', 'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':  ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}
HIP_DEPTH = hip[['sample_name', 'num_processed']].drop_duplicates()['num_processed'].median()
print('hip median depth (reads processed):', int(HIP_DEPTH))

rows, stat_rows = [], []
comp = {}   # (group, region) -> dict(panel -> {vclass: frac}), ncells, depleted
for name in GROUPS:
    for region in ['ctx', 'hip']:
        d = cells_of(name, region)
        if d.empty: continue
        meta, V, G = percell(d)
        n = len(meta)
        rec = dict(group=name, region=region, n=n, median_depth=meta['num_processed'].median())
        fr = {}
        for panel, vs in PANEL_V.items():
            tot = V.reindex(columns=vs).fillna(0).sum().sum()
            fr[panel] = {v: (V[v].sum() / tot if (v in V and tot > 0) else 0.0) for v in vs}
        for fam, genes in FAM.items():
            cpm = G[genes].sum(axis=1) / meta['num_processed'] * 1e6
            rec[f'{fam}_median_cpm'] = cpm.median(); rec[f'{fam}_pct_zero'] = 100 * (cpm == 0).mean()
            # gene shares within family
            gs = G[genes].sum(); rec.update({f'{g}_share': gs[g] / gs.sum() if gs.sum() else np.nan for g in genes})
        # depth-matched zero rates: thin cortex reads to the hip median depth
        if region == 'ctx':
            p = np.clip(HIP_DEPTH / meta['num_processed'], 0, 1)
            Gt = pd.DataFrame({g: rng.binomial(G[g].round().astype(int).values, p.values) for g in G.columns}, index=G.index)
            for fam, genes in FAM.items():
                rec[f'{fam}_pct_zero_thinned'] = 100 * (Gt[genes].sum(axis=1) == 0).mean()
        rows.append(rec); comp[(name, region)] = (fr, n)
        # per-cell variant fractions for stats (cells with >=20 reads on the gene)
        for panel, vs in PANEL_V.items():
            for v in vs:
                if v not in V: continue
                gene = [g for g, pv in PANEL_OF.items() if pv == panel]
                # fraction of the family's reads in this cell
                fam_reads = V.reindex(columns=vs).fillna(0).sum(axis=1)
                ok = fam_reads >= 20
                stat_rows.append(pd.DataFrame({'group': name, 'region': region, 'panel': panel, 'vclass': v,
                                               'frac': (V.loc[ok, v] / fam_reads[ok]).values}))
R = pd.DataFrame(rows); R.to_csv(f'{OUT}/region_summary.csv', index=False)
S = pd.concat(stat_rows)
# MWU per group x variant, cortex vs hip, on per-cell fractions
tests = []
for (g, panel, v), x in S.groupby(['group', 'panel', 'vclass']):
    a = x[x.region == 'ctx'].frac; b = x[x.region == 'hip'].frac
    if len(a) >= 5 and len(b) >= 5 and (a.mean() > 0.02 or b.mean() > 0.02):
        tests.append(dict(group=g, panel=panel, vclass=v, n_ctx=len(a), n_hip=len(b), mean_ctx=a.mean(), mean_hip=b.mean(),
                          diff_hip_minus_ctx=b.mean() - a.mean(), p=mannwhitneyu(a, b).pvalue))
T = pd.DataFrame(tests); T['p_BH'] = multipletests(T.p, method='fdr_bh')[1]
T = T.sort_values('p_BH'); T.to_csv(f'{OUT}/region_variant_tests.csv', index=False)

# ---- figure: paired stacked bars (ctx above hip) per matched group, 3 panels, + depleted/zero panel
groups = [g for g in GROUPS if (g, 'ctx') in comp and (g, 'hip') in comp]
fig, axes = plt.subplots(1, 4, figsize=(22, 1.05 * len(groups) + 3.5), gridspec_kw={'width_ratios': [3, 3, 3, 1.6]})
ypos = {}
y = 0.0
for g in groups[::-1]:
    ypos[g] = y; y += 1.0
BH = 0.36
for ax, panel in zip(axes[:3], PANEL_V):
    for g in groups:
        for k, region in enumerate(['ctx', 'hip']):
            fr, n = comp[(g, region)]
            yy = ypos[g] + (0.21 if region == 'ctx' else -0.21)
            left = 0.0
            for v, f in sorted(fr[panel].items(), key=lambda kv: -kv[1]):
                if f <= 0: continue
                base, kind, pcolor = VSPEC[v]
                ax.barh(yy, f, left=left, color=base, edgecolor=pcolor, hatch=HATCH[kind], height=BH, linewidth=0.3)
                left += f
            ax.text(1.01, yy, f'{"V1" if region == "ctx" else "HPF"} n={n:,}', va='center', fontsize=7, color='#555')
    ax.set_xlim(0, 1.25); ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.set_yticks([ypos[g] for g in groups]); ax.set_yticklabels(groups if ax is axes[0] else [], fontsize=9)
    ax.set_xticks([]); ax.set_title(panel, fontsize=13, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(False)
# zero / depleted panel
ax = axes[3]
Rz = R.set_index(['group', 'region'])
for g in groups:
    for k, (region, col, lab) in enumerate([('ctx', 'Kv4_pct_zero', 'V1 (full depth)'), ('ctx', 'Kv4_pct_zero_thinned', 'V1 thinned to HPF depth'), ('hip', 'Kv4_pct_zero', 'HPF')]):
        val = Rz.loc[(g, region), col]
        yy = ypos[g] + (0.26 - 0.26 * k)
        ax.barh(yy, val, height=0.24, color=['#9ecae1', '#3182bd', '#D62728'][k], label=lab if g == groups[0] else None)
        ax.text(val + 0.6, yy, f'{val:.0f}', va='center', fontsize=7)
ax.set_ylim(-0.6, len(groups) - 0.4); ax.set_yticks([]); ax.set_xlabel('% cells with zero Kv4 reads', fontsize=9)
ax.set_title('Kv4-null cells', fontsize=12, fontweight='bold'); ax.legend(fontsize=7.5, loc='lower right', frameon=False)
for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
# legends
def legend_patches(vlist):
    pats = []
    for v in vlist:
        if v == BLANK: pats.append(mpatches.Patch(facecolor='none', edgecolor='none', label=' ')); continue
        base, kind, pcolor = VSPEC[v]
        pats.append(mpatches.Patch(facecolor=base, edgecolor=pcolor, hatch=HATCH[kind], label=v, linewidth=0.3))
    return pats
for ax, panel in zip(axes[:3], PANEL_V):
    groups_l = LEGEND_GROUPS[panel]
    xs = [0.5] if len(groups_l) == 1 else [0.28, 0.76]
    legs = [ax.legend(handles=legend_patches([v for v in vl if v in PANEL_V[panel] or v == BLANK]), loc='upper center', bbox_to_anchor=(x, -0.01), ncol=1,
                      fontsize=7.5, frameon=False, handlelength=3.0, handleheight=2.2) for x, vl in zip(xs, groups_l)]
    for leg in legs[:-1]: ax.add_artist(leg)
fig.suptitle('Cortex (V1, Tasic 2016) vs hippocampal formation (Allen SSv4, Yao 2021): matched cell classes\n'
             'upper bar = V1, lower bar = HPF; pooled-read variant fractions, most-common segment first. '
             'Right: % Kv4-null cells, with V1 reads binomially thinned to the HPF median depth for a fair comparison', fontsize=11.5)
fig.subplots_adjust(left=0.11, right=0.985, top=0.9, bottom=0.42, wspace=0.12)
fig.savefig(f'{OUT}/region_comparison.png', dpi=190)

pd.set_option('display.width', 300)
cols = ['group', 'region', 'n', 'median_depth', 'Kv4_median_cpm', 'KChIP_median_cpm', 'DPP_median_cpm', 'Kv4_pct_zero', 'Kv4_pct_zero_thinned', 'KChIP_pct_zero', 'DPP_pct_zero']
print(R[cols].round(1).to_string(index=False))
print('\nsignificant variant differences (BH < 0.05), by |diff|:')
print(T[T.p_BH < 0.05].sort_values('diff_hip_minus_ctx', key=abs, ascending=False).round(3).head(40).to_string(index=False))
