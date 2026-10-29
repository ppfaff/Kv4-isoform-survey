"""Three-way comparison: Tasic 2016 VISp (100 bp SE) | Allen SSv4 VISp (50 bp PE) | Allen SSv4 hippocampal formation.
Same region, different chemistry  -> chemistry effect;   same chemistry, different region -> biology."""
import sys; sys.path.insert(0, '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip')
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import mannwhitneyu
from hip_common import *

OUT = '/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/hip/compare'
hip = load_hip(); ctx = load_ctx()
v = pd.read_pickle('/home/claude/kv4/visp/visp_per_cell.pkl')
v['vclass'] = v['transcript'].map(final_map); v = v.dropna(subset=['vclass']).copy()
v = v[~((v['biotype'] == 'retained_intron') & v['vclass'].str.contains('other|fragment', case=False))].copy()
v['panel'] = v['gene'].map(PANEL_OF); v['cls'] = v['subclass_label'].map(lambda s: 'GABA' if s in GABA + ['Meis2'] else ('NonN' if s in NONN + ['Endo', 'SMC-Peri', 'VLMC'] else 'Glut'))
v = v.rename(columns={'subclass_label': 'subclass', 'cluster_label': 'cluster'})

CTX_T = lambda prefix: [t for t in ctx['subclass'].unique() if t.startswith(prefix)]
GROUPS = {   # name: (Tasic2016 types, Allen VISp subclasses, Allen HPF subclasses)
    'Lamp5':   (['Ndnf Car4', 'Ndnf Cxcl14', 'Igtp', 'Smad3'], ['Lamp5'], ['Lamp5']),
    'Sncg':    (['Sncg', 'Vip Sncg'], ['Sncg'], ['Sncg']),
    'Vip':     ([t for t in CTX_T('Vip') if t != 'Vip Sncg'], ['Vip'], ['Vip']),
    'Sst':     (CTX_T('Sst'), ['Sst', 'Sst Chodl'], ['Sst', 'Sst Chodl']),
    'Pvalb':   (CTX_T('Pvalb'), ['Pvalb'], ['Pvalb']),
    'L2/3 IT': (['L2/3 Ptgs2', 'L2 Ngb'], ['L2/3 IT CTX'], ['L2/3 IT RHP', 'L2/3 IT PPP']),
    'L6 CT':   (['L6a Mgp', 'L6a Sla', 'L6a Car12'], ['L6 CT CTX'], ['L6 CT CTX']),
    'L6b':     (CTX_T('L6b'), ['L6b CTX'], ['L6b CTX']),
    # cortex-only classes (no hippocampal counterpart): chemistry comparison only
    'L4':      (CTX_T('L4'), ['L4/5 IT CTX'], []),
    'L5 IT':   (CTX_T('L5a'), ['L5 IT CTX'], []),
    'L5 PT':   (CTX_T('L5b'), ['L5 PT CTX'], []),
    'L5/6 NP': (['L5 Ucma'], ['L5/6 NP CTX'], []),
    'L6 IT':   (['L6a Syt17'], ['L6 IT CTX'], []),
}
GROUP_CLS = {g: ('GABA' if g in ['Lamp5', 'Sncg', 'Vip', 'Sst', 'Pvalb'] else 'Glut') for g in GROUPS}
DS = {'T16': ctx, 'A-VISp': v, 'A-HPF': hip}
DSCOL = {'T16': '#999999', 'A-VISp': '#2060B0', 'A-HPF': '#D62728'}
PANEL_V = {
    'Kv4':  ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short', 'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich', 'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':  ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}

def cell_fracs(d, sel):
    d = d[d['subclass'].isin(sel)]
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    out = {}
    for panel, vs in PANEL_V.items():
        X = V.reindex(columns=vs).fillna(0.0); tot = X.sum(axis=1); ok = tot >= 20
        out[panel] = (X[ok].div(tot[ok], axis=0), V.reindex(columns=vs).fillna(0).sum() / max(V.reindex(columns=vs).fillna(0).sum().sum(), 1))
    return out, d['sample_name'].nunique()

rows = []; pooled = {}
for g, sels in GROUPS.items():
    for ds, sel in zip(DS, sels):
        if not sel: continue
        fr, n = cell_fracs(DS[ds], sel)
        pooled[(g, ds)] = ({p: fr[p][1].to_dict() for p in PANEL_V}, n)
        for p in PANEL_V:
            F = fr[p][0]
            for vc in PANEL_V[p]:
                rows.append(dict(group=g, dataset=ds, panel=p, vclass=vc, n_cells=len(F), mean_frac=F[vc].mean() if len(F) else np.nan))
R = pd.DataFrame(rows)
W = R.pivot_table(index=['group', 'panel', 'vclass'], columns='dataset', values='mean_frac').reset_index()
W['chem_effect (A-VISp − T16)'] = W['A-VISp'] - W['T16']
W['region_effect (A-HPF − A-VISp)'] = W['A-HPF'] - W['A-VISp']
W.to_csv(f'{OUT}/threeway_variant_fracs.csv', index=False)

# significance of region effect within the same chemistry (per-cell MWU), BH over all tests
from statsmodels.stats.multitest import multipletests
tests = []
for g, sels in GROUPS.items():
    fa, _ = cell_fracs(v, sels[1]); ft, _ = cell_fracs(ctx, sels[0])
    fh = cell_fracs(hip, sels[2])[0] if sels[2] else None
    for p in PANEL_V:
        for vc in PANEL_V[p]:
            a, t = fa[p][0][vc], ft[p][0][vc]
            h = fh[p][0][vc] if fh is not None else pd.Series(dtype=float)
            if len(a) >= 5 and len(h) >= 5 and (a.mean() > 0.02 or h.mean() > 0.02):
                tests.append(dict(group=g, panel=p, vclass=vc, n_visp=len(a), n_hpf=len(h), mean_visp=a.mean(), mean_hpf=h.mean(),
                                  region_diff=h.mean() - a.mean(), p_region=mannwhitneyu(a, h).pvalue,
                                  mean_t16=t.mean() if len(t) else np.nan, chem_diff=(a.mean() - t.mean()) if len(t) >= 5 else np.nan,
                                  p_chem=mannwhitneyu(a, t).pvalue if len(t) >= 5 else np.nan))
T = pd.DataFrame(tests)
T['pBH_region'] = multipletests(T.p_region, method='fdr_bh')[1]
m = T.p_chem.notna(); T.loc[m, 'pBH_chem'] = multipletests(T.p_chem[m], method='fdr_bh')[1]
T.to_csv(f'{OUT}/threeway_tests.csv', index=False)

# ---- summary of consistency: variants whose chemistry effect has the same sign in >=6 of 8 groups
pd.set_option('display.width', 300)
print('=== variants with a consistent CHEMISTRY shift (same sign in >=6 of 8 matched groups; mean shift) ===')
cs = W.groupby('vclass')['chem_effect (A-VISp − T16)'].agg(['mean', lambda x: (np.sign(x) == np.sign(x.mean())).sum(), 'count'])
cs.columns = ['mean_shift', 'n_same_sign', 'n_groups']; cs = cs[(cs.n_same_sign >= 6) & (cs.mean_shift.abs() > 0.02)].sort_values('mean_shift')
print(cs.round(3).to_string())
print('\n=== variants with a consistent REGION shift (same chemistry, same sign in >=6 of 8) ===')
rs = W.groupby('vclass')['region_effect (A-HPF − A-VISp)'].agg(['mean', lambda x: (np.sign(x) == np.sign(x.mean())).sum(), 'count'])
rs.columns = ['mean_shift', 'n_same_sign', 'n_groups']; rs = rs[(rs.n_same_sign >= 6) & (rs.mean_shift.abs() > 0.02)].sort_values('mean_shift')
print(rs.round(3).to_string())
print('\n=== significant region differences at matched chemistry (BH<0.05, |diff|>0.1) ===')
print(T[(T.pBH_region < 0.05) & (T.region_diff.abs() > 0.1)].sort_values('region_diff', key=abs, ascending=False)[['group', 'vclass', 'n_visp', 'n_hpf', 'mean_t16', 'mean_visp', 'mean_hpf', 'region_diff', 'pBH_region']].round(3).to_string(index=False))

# ---- figure: three stacked bars per group (T16 / A-VISp / A-HPF)
groups = list(GROUPS)
fig, axes = plt.subplots(1, 3, figsize=(21, 1.15 * len(groups) + 4.0))
ypos = {g: i for i, g in enumerate(groups[::-1])}
BH = 0.27; OFF = {'T16': 0.3, 'A-VISp': 0.0, 'A-HPF': -0.3}
for ax, panel in zip(axes, PANEL_V):
    for g in groups:
        for ds in DS:
            if (g, ds) not in pooled: continue
            fr, n = pooled[(g, ds)]; yy = ypos[g] + OFF[ds]; left = 0.0
            for vc, f in sorted(fr[panel].items(), key=lambda kv: -kv[1]):
                if f <= 0: continue
                base, kind, pcolor = VSPEC[vc]
                ax.barh(yy, f, left=left, color=base, edgecolor=pcolor, hatch=HATCH[kind], height=BH, linewidth=0.3)
                left += f
            ax.text(1.01, yy, f'{ds} n={n}', va='center', fontsize=7, color=DSCOL[ds], fontweight='bold')
    ax.axhline(ypos['L6b'] - 0.5, color='#BBBBBB', linewidth=0.8)
    ax.set_xlim(0, 1.22); ax.set_ylim(-0.6, len(groups) - 0.4)
    ax.set_yticks([ypos[g] for g in groups]); ax.set_yticklabels(groups if ax is axes[0] else [], fontsize=14, fontweight='bold')
    if ax is axes[0]:
        for lbl, g in zip(ax.get_yticklabels(), groups): lbl.set_color(CLS_COLOR[GROUP_CLS[g]])
    ax.set_xticks([]); ax.set_title(panel, fontsize=13, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(False)
def legend_patches(vlist):
    pats = []
    for vc in vlist:
        if vc == BLANK: pats.append(mpatches.Patch(facecolor='none', edgecolor='none', label=' ')); continue
        base, kind, pcolor = VSPEC[vc]
        pats.append(mpatches.Patch(facecolor=base, edgecolor=pcolor, hatch=HATCH[kind], label=vc, linewidth=0.3))
    return pats
for ax, panel in zip(axes, PANEL_V):
    gl = LEGEND_GROUPS[panel]; xs = [0.5] if len(gl) == 1 else [0.28, 0.76]
    legs = [ax.legend(handles=legend_patches([x for x in vl if x in PANEL_V[panel] or x == BLANK]), loc='upper center', bbox_to_anchor=(x, -0.01),
                      ncol=1, fontsize=7.5, frameon=False, handlelength=3.0, handleheight=2.2) for x, vl in zip(xs, gl)]
    for leg in legs[:-1]: ax.add_artist(leg)
fig.suptitle('Variant composition in matched cell classes — three datasets per class\n'
             'grey label = Tasic 2016 VISp (100 bp single-end) · blue = Allen SSv4 VISp (50 bp paired-end) · red = Allen SSv4 hippocampal formation\n'
             'T16 vs A-VISp isolates the chemistry/read-length effect; A-VISp vs A-HPF is the regional comparison at matched chemistry. '
             'Purple = inhibitory, blue = excitatory classes; L4–L6 IT rows have no hippocampal counterpart', fontsize=11.5)
fig.subplots_adjust(left=0.08, right=0.985, top=0.93, bottom=0.2, wspace=0.1)
fig.savefig(f'{OUT}/threeway_variant_composition.png', dpi=190)
print('\nfigure written')
