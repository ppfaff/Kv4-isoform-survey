#!/usr/bin/env python3
"""Per-nucleus Kv4-complex N-terminal variant composition from the ThalamoSeq pools (tier 2, relay side).

Input:  tables/thalamoseq_isoform_per_sample.csv  (pipeline/collapse_thalamoseq.py; one row per transcript per pool)
        tables/thalamoseq_run_qc.csv
Outputs:
  tables/thalamoseq_variant_fracs_by_nucleus.csv   gene, vclass, nucleus, loop_group, n_pools, frac_mean, frac_sem, reads
  tables/thalamoseq_gene_level_by_nucleus.csv      gene TPM per nucleus (mean/sem), Kv4 paralog share, KChIP & DPP composition
  tables/thalamoseq_fo_vs_ho.csv                   first-order vs higher-order contrast per variant (Welch t on pool fractions)
  figures/thalamoseq_variants_by_nucleus.png       stacked variant composition per nucleus, Kv4/KChIP/DPP genes
  figures/thalamoseq_gene_level.png                Kcnd2/Kcnd3 share and partner composition per nucleus vs the 10x map

Pools are the replicate unit (3-9 per nucleus). Variant fractions are reads-weighted within a pool, then averaged
over pools; a pool contributes to a gene only if it has >= MIN_READS reads on that gene.
Usage: python3 thalamus/thalamoseq_variants.py [TABLES_DIR] [OUT_DIR] [--min-reads 200]
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import VSPEC, HATCH  # noqa: E402
from th_common import loop_group_of_nucleus  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('tables', nargs='?', default='tables'); ap.add_argument('out', nargs='?', default='.')
ap.add_argument('--min-reads', type=float, default=200)
args = ap.parse_args()
os.makedirs(f'{args.out}/tables', exist_ok=True); os.makedirs(f'{args.out}/figures', exist_ok=True)

A = pd.read_csv(f'{args.tables}/thalamoseq_isoform_per_sample.csv')
A['vclass'] = A.vclass.fillna('other/fragment')
# nucleus naming: VB = VPM/VPL; map to loop groups through th_common (first token)
NUC_ALIAS = {'VB': 'VPM'}
A['loop_group'] = A.nucleus.map(lambda n: loop_group_of_nucleus(NUC_ALIAS.get(n, n)))
A.loc[A.nucleus.isin(['SPA']), 'loop_group'] = 'ZI / LGv / IGL (non-relay)'
GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']

# ---------------- variant fractions per pool -> per nucleus
P = A.groupby(['sample_id', 'nucleus', 'loop_group', 'gene', 'vclass'], as_index=False).agg(reads=('reads', 'sum'))
P['gene_reads'] = P.groupby(['sample_id', 'gene']).reads.transform('sum')
P = P[P.gene_reads >= args.min_reads].copy()
P['frac'] = P.reads / P.gene_reads
# complete the grid so absent variants count as 0 in the pool mean
grid = P[['sample_id', 'nucleus', 'loop_group', 'gene']].drop_duplicates().merge(P[['gene', 'vclass']].drop_duplicates(), on='gene')
P = grid.merge(P[['sample_id', 'gene', 'vclass', 'frac', 'reads']], on=['sample_id', 'gene', 'vclass'], how='left').fillna({'frac': 0.0, 'reads': 0.0})
V = P.groupby(['gene', 'vclass', 'nucleus', 'loop_group'], as_index=False).agg(
    n_pools=('sample_id', 'nunique'), frac_mean=('frac', 'mean'), frac_sem=('frac', lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan),
    reads=('reads', 'sum'))
V.to_csv(f'{args.out}/tables/thalamoseq_variant_fracs_by_nucleus.csv', index=False)

# ---------------- gene level per nucleus
G = A.groupby(['sample_id', 'nucleus', 'loop_group', 'gene'], as_index=False).agg(tpm=('TPM', 'sum'), reads=('reads', 'sum'))
GN = G.groupby(['nucleus', 'loop_group', 'gene'], as_index=False).agg(n_pools=('sample_id', 'nunique'), tpm_mean=('tpm', 'mean'),
                                                                     tpm_sem=('tpm', lambda x: x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan))
W = GN.pivot(index='nucleus', columns='gene', values='tpm_mean').reindex(columns=GENES).fillna(0)
kv = W[['Kcnd1', 'Kcnd2', 'Kcnd3']].sum(axis=1); kc = W[['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']].sum(axis=1); dp = W[['Dpp6', 'Dpp10']].sum(axis=1)
S = W.copy()
for g in ['Kcnd1', 'Kcnd2', 'Kcnd3']:
    S[f'share_{g}'] = W[g] / kv
for g in ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']:
    S[f'KChIP_comp_{g}'] = W[g] / kc
for g in ['Dpp6', 'Dpp10']:
    S[f'DPP_comp_{g}'] = W[g] / dp
S['KChIP_per_Kv4'] = kc / kv; S['DPP_per_Kv4'] = dp / kv
S['loop_group'] = S.index.map(lambda n: loop_group_of_nucleus(NUC_ALIAS.get(n, n)))
S['n_pools'] = G.groupby('nucleus').sample_id.nunique()
S.to_csv(f'{args.out}/tables/thalamoseq_gene_level_by_nucleus.csv')

# ---------------- first-order vs higher-order contrast
FO = {'LGd', 'VB', 'VPMpc', 'VL', 'VA', 'VM', 'AV', 'AD', 'AM', 'LD'}      # primary sensory/motor/anterior relays
HO = {'PO', 'LP', 'MD', 'CM', 'CL', 'PCN', 'IMD', 'PVT', 'PF', 'RE', 'IAD', 'SPA'}
P['fo_ho'] = np.where(P.nucleus.isin(FO), 'FO', np.where(P.nucleus.isin(HO), 'HO', 'other'))
rows = []
for (g, v), sub in P[P.fo_ho != 'other'].groupby(['gene', 'vclass']):
    a, b = sub[sub.fo_ho == 'FO'].frac, sub[sub.fo_ho == 'HO'].frac
    if len(a) > 2 and len(b) > 2 and max(a.mean(), b.mean()) > 0.02:
        # Welch t, no scipy dependency
        se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
        t = (a.mean() - b.mean()) / se if se > 0 else np.nan
        rows.append(dict(gene=g, vclass=v, FO_mean=a.mean(), HO_mean=b.mean(), diff=a.mean() - b.mean(), welch_t=t, n_FO=len(a), n_HO=len(b)))
C = pd.DataFrame(rows).sort_values('welch_t', key=abs, ascending=False)
C.to_csv(f'{args.out}/tables/thalamoseq_fo_vs_ho.csv', index=False)

# ---------------- figure 1: variant composition per nucleus
nuc_order = sorted(S.index, key=lambda n: (['first-order sensory relay', 'motor relay', 'anterior / limbic', 'higher-order relay',
                                            'intralaminar / midline', 'ZI / LGv / IGL (non-relay)', 'other'].index(S.loc[n, 'loop_group'])
                                           if S.loc[n, 'loop_group'] in ['first-order sensory relay', 'motor relay', 'anterior / limbic', 'higher-order relay',
                                                                         'intralaminar / midline', 'ZI / LGv / IGL (non-relay)'] else 6, n))
show = ['Kcnd3', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Kcnip1', 'Dpp10']
fig, axes = plt.subplots(len(show), 1, figsize=(0.55 * len(nuc_order) + 3, 2.4 * len(show)), sharex=True)
handles = {}
for ax, g in zip(axes, show):
    sub = V[V.gene == g]
    for i, n in enumerate(nuc_order):
        s = sub[sub.nucleus == n].sort_values('frac_mean', ascending=False); b = 0
        for _, r in s.iterrows():
            if r.frac_mean < 0.01:
                continue
            face, hatch, hcol = VSPEC.get(r.vclass, ('#DDDDDD', 'solid', '#DDDDDD'))
            ax.bar(i, r.frac_mean, 0.8, bottom=b, color=face, hatch=HATCH[hatch], edgecolor=hcol if HATCH[hatch] else 'k', lw=0.4)
            handles.setdefault(r.vclass, Patch(facecolor=face, hatch=HATCH[hatch], edgecolor=hcol if HATCH[hatch] else 'k'))
            b += r.frac_mean
    ax.set_ylim(0, 1); ax.set_ylabel(g, fontsize=10)
    if g == 'Kcnd3':
        ax.set_title('ThalamoSeq pooled relay neurons: N-terminal variant composition by nucleus (mean of pools)', loc='left', fontsize=11)
axes[-1].set_xticks(range(len(nuc_order))); axes[-1].set_xticklabels([f'{n}\n{S.loc[n, "loop_group"].split()[0]}' for n in nuc_order], fontsize=8)
fig.legend(handles.values(), handles.keys(), fontsize=10, frameon=False, loc='lower center', ncol=8, bbox_to_anchor=(0.5, -0.02))
fig.savefig(f'{args.out}/figures/thalamoseq_variants_by_nucleus.png', dpi=160, bbox_inches='tight'); plt.close(fig)

# ---------------- figure 2: gene level vs 10x calibrated map
fig, axes = plt.subplots(1, 3, figsize=(16, 0.32 * len(nuc_order) + 2.5), sharey=True)
y = np.arange(len(nuc_order)); Sx = S.loc[nuc_order]
for ax, fam, cols, lab, title in [(axes[0], ['Kcnd1', 'Kcnd2', 'Kcnd3'], ['#4FA8D8', '#E8822D', '#B9A3D9'], 'share', 'Kv4 paralog share (TPM)'),
                                  (axes[1], ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'], ['#7B4FA3', '#4FA8D8', '#4CAF50', '#E8822D'], 'KChIP_comp', 'KChIP composition'),
                                  (axes[2], ['Dpp6', 'Dpp10'], ['#7B4FA3', '#E8822D'], 'DPP_comp', 'DPP composition')]:
    left = np.zeros(len(nuc_order))
    for g, c in zip(fam, cols):
        v = Sx[f'{lab}_{g}'].values
        ax.barh(y, v, left=left, color=c, label=g, edgecolor='k', lw=0.3); left += v
    ax.set_xlim(0, 1); ax.set_title(title, fontsize=10); ax.legend(fontsize=10, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=4)
axes[0].set_yticks(y); axes[0].set_yticklabels([f'{n}  (n={int(Sx.loc[n, "n_pools"])})' for n in nuc_order], fontsize=8); axes[0].invert_yaxis()
# overlay the 10x calibrated Kcnd2 share where the nucleus exists
try:
    T = pd.read_csv(f'{args.tables}/th10x_calibrated_nuclei.csv', index_col=0)
    alias10 = {'VB': 'VPM VPL LGd', 'VPMpc': 'VPMpc', 'PCN': 'PCN CL'}
    for i, n in enumerate(nuc_order):
        k = alias10.get(n, n)
        if k in T.index:
            axes[0].plot(T.loc[k, 'share_Kcnd2'] + T.loc[k, 'share_Kcnd1'], i, marker='|', color='k', ms=12, mew=2)
    axes[0].text(0.02, len(nuc_order) - 0.3, '| = 10x calibrated Kcnd1+Kcnd2 share', fontsize=7)
except Exception:
    pass
fig.suptitle('ThalamoSeq gene level (full-length pools) by nucleus', y=1.01)
fig.savefig(f'{args.out}/figures/thalamoseq_gene_level.png', dpi=160, bbox_inches='tight'); plt.close(fig)

pd.set_option('display.width', 220)
print(S[['loop_group', 'n_pools', 'share_Kcnd2', 'share_Kcnd3', 'KChIP_comp_Kcnip3', 'KChIP_comp_Kcnip4', 'DPP_comp_Dpp6']].round(2).sort_values('loop_group').to_string())
print('\nFO vs HO variant contrasts (|t| > 2):'); print(C[C.welch_t.abs() > 2].round(3).to_string(index=False))
