#!/usr/bin/env python3
"""Activity state of Kv4-null vs Kv4-low vs Kv4-high DG granule cells.

Question: is the Kv4-null state (no mature Kcnd2 mRNA, burst state absent) an activity state? Kv4.2 is
activity-downregulated, so if it is, the null cells should be the HIGH-activity granule cells (high Arc/Fos/Egr1/Npas4,
ΔFosB), not the quiet ones.

Inputs:
  tables/allen_activity_markers.csv   per-cell counts for IEGs etc. from the Allen matrix (extract_markers.py, GENES=...)
  tables/dg_cells_markers.csv         the survey's 2,458 DG cells with Kv4 group
  tables/dg_fosb_junctions.tsv        (optional) Fosb-locus junctions from star_dg/*/SJ.out.tab for the 746 aligned cells
  tables/dg_readlevel_runlist.tsv     group labels for those 746
Outputs: tables/dg_activity_by_group.csv, tables/dg_activity_per_cell.csv, figures/dg_activity.png
Usage:   python3 scripts/dg_activity.py [TABLES_DIR] [OUT_DIR]
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)
GRPS = ['Kv4-null', 'Kv4-low', 'Kv4-high']; COL = {'Kv4-null': '#D6372E', 'Kv4-low': '#E8822D', 'Kv4-high': '#2E6FBA'}
IEG = ['Arc', 'Fos', 'Fosb', 'Egr1', 'Egr2', 'Egr3', 'Egr4', 'Npas4', 'Junb', 'Nr4a1', 'Nr4a2', 'Ier2', 'Dusp1', 'Bdnf', 'Homer1', 'Fosl2']
OTHER = ['Calb1', 'Penk', 'Kcnd2', 'Kcnip4', 'Dpp6']

A = pd.read_csv(f'{TABLES}/allen_activity_markers.csv')
D = pd.read_csv(f'{TABLES}/dg_cells_markers.csv')[['sample_name', 'grp', 'cluster_label', 'num_processed']]
M = D.merge(A, on='sample_name', how='inner')
tot = M['total_exon_counts'] if 'total_exon_counts' in M else M['num_processed']
genes = [g for g in IEG + OTHER if g in M.columns]
for g in genes:
    M[g + '_cpm'] = M[g] / tot * 1e6
M['ieg_index'] = np.log2(M[[g + '_cpm' for g in IEG if g in M.columns]] + 1).mean(axis=1)
print(f'{len(M)} DG cells with activity markers; groups: {dict(M.grp.value_counts())}')

rows = []
for g in genes + ['ieg_index']:
    col = g if g == 'ieg_index' else g + '_cpm'
    med = M.groupby('grp')[col].median().reindex(GRPS)
    det = M.groupby('grp')[col].apply(lambda x: (x > 0).mean() * 100).reindex(GRPS)
    p_nh = stats.mannwhitneyu(M[M.grp == 'Kv4-null'][col], M[M.grp == 'Kv4-high'][col]).pvalue
    p_nl = stats.mannwhitneyu(M[M.grp == 'Kv4-null'][col], M[M.grp == 'Kv4-low'][col]).pvalue
    rows.append((g, *med.values, *det.values, med['Kv4-null'] / med['Kv4-high'] if med['Kv4-high'] > 0 else np.nan, p_nh, p_nl))
R = pd.DataFrame(rows, columns=['gene', 'med_null', 'med_low', 'med_high', 'pct_det_null', 'pct_det_low', 'pct_det_high', 'null_over_high', 'p_null_vs_high', 'p_null_vs_low'])
R.to_csv(f'{OUT}/tables/dg_activity_by_group.csv', index=False)
M.to_csv(f'{OUT}/tables/dg_activity_per_cell.csv', index=False)
pd.set_option('display.width', 220)
print(R.round(3).to_string(index=False))

# ΔFosB from junctions (optional)
F = None
jp = f'{TABLES}/dg_fosb_junctions.tsv'
if os.path.exists(jp) and os.path.getsize(jp):
    J = pd.read_csv(jp, sep='\t', header=None, names=['sample_name', 'chr', 'start', 'end', 'strand', 'motif', 'annot', 'uniq', 'multi', 'overhang'])
    Rl = pd.read_csv(f'{TABLES}/dg_readlevel_runlist.tsv', sep='\t')[['sample_name', 'grp']]
    J = J.merge(Rl, on='sample_name')
    tot_j = J.groupby(['chr', 'start', 'end', 'annot']).uniq.sum().reset_index().sort_values('uniq', ascending=False)
    print('\nFosb-locus junctions (all 746 cells pooled):'); print(tot_j.head(12).to_string(index=False))
    # the ΔFosB junction is the short (~140-nt) internal junction inside exon 4; pick the most-used junction with span 100-200 nt
    tot_j['span'] = tot_j.end - tot_j.start + 1
    cand = tot_j[(tot_j.span.between(100, 200))]
    if len(cand):
        ds, de = cand.iloc[0].start, cand.iloc[0].end
        print(f'ΔFosB candidate junction: {cand.iloc[0].chr}:{ds}-{de} (span {de - ds + 1} nt, {int(cand.iloc[0].uniq)} reads pooled)')
        J['is_delta'] = (J.start == ds) & (J.end == de)
        F = J.groupby(['sample_name', 'grp']).apply(lambda d: pd.Series({'delta': d[d.is_delta].uniq.sum(), 'other': d[~d.is_delta].uniq.sum()})).reset_index()
        F['frac_delta'] = F.delta / (F.delta + F.other)
        print('\nΔFosB junction reads by group (cells aligned by run_dg_star.sh):')
        print(F.groupby('grp').agg(n=('delta', 'size'), delta_sum=('delta', 'sum'), other_sum=('other', 'sum'),
                                   pct_cells_delta=('delta', lambda x: (x > 0).mean() * 100)).reindex(GRPS).to_string())
        F.to_csv(f'{OUT}/tables/dg_fosb_delta_per_cell.csv', index=False)

# figure
show = [g for g in ['Arc', 'Fos', 'Fosb', 'Egr1', 'Npas4', 'Junb', 'Nr4a1', 'Bdnf', 'Calb1', 'Penk'] if g in genes]
fig, axes = plt.subplots(1, 2, figsize=(14, 4.2), gridspec_kw={'width_ratios': [3, 1]})
ax = axes[0]; x = np.arange(len(show)); w = 0.27
for i, grp in enumerate(GRPS):
    v = [np.log2(M[M.grp == grp][g + '_cpm'] + 1).mean() for g in show]
    e = [np.log2(M[M.grp == grp][g + '_cpm'] + 1).sem() for g in show]
    ax.bar(x + (i - 1) * w, v, w, yerr=e, color=COL[grp], label=f'{grp} (n={int((M.grp == grp).sum())})', capsize=1.5)
ax.set_xticks(x); ax.set_xticklabels(show, rotation=45, fontsize=9); ax.set_ylabel('mean log2(CPM+1)'); ax.legend(fontsize=8, frameon=False)
ax.set_title('Activity-regulated genes by Kv4 group (Allen exon counts, DG granule cells)')
ax = axes[1]
for i, grp in enumerate(GRPS):
    v = M[M.grp == grp].ieg_index
    ax.scatter(np.full(len(v), i) + np.random.uniform(-0.25, 0.25, len(v)), v, s=4, alpha=0.4, color=COL[grp])
    ax.plot([i - 0.3, i + 0.3], [v.median()] * 2, color='k', lw=1.5)
ax.set_xticks(range(3)); ax.set_xticklabels(GRPS, fontsize=8); ax.set_ylabel('IEG index (mean log2 CPM over %d genes)' % len([g for g in IEG if g in M.columns])); ax.set_title('Per-cell activity index')
fig.tight_layout(); fig.savefig(f'{OUT}/figures/dg_activity.png', dpi=200)
print('\ndone')
