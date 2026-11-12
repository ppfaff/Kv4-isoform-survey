#!/usr/bin/env python3
"""Calibrate the ABC Atlas 10x thalamus summaries for gene-length (intronic) bias and make the main figure.

Anchor (2026-11-12): ThalamoSeq full-length relay-neuron TPM (tables/thalamoseq_gene_level_by_nucleus.csv) when present -
neurons only and the same relay populations - falling back to the bulk thalamus exon CPM otherwise.

The WMB 10x counts include intronic reads, so the Kv4-complex genes - all very long (Kcnd2 ~0.5 Mb, Kcnip4
~1.2 Mb, Dpp10 ~1.4 Mb, Dpp6 ~1 Mb; Kcnd1 ~10 kb) - are inflated relative to short genes and to each other.
Within-gene comparisons across cell types are unaffected; cross-gene quantities (Kv4 paralog share, KChIP and
DPP partner composition, subunit:channel ratios) are not. We calibrate each gene with one factor,
    f_gene = (bulk thalamus WT exon-level CPM) / (cell-weighted mean 10x CPM over thalamic neurons),
from the lab's V404M bulk RNA-seq (tables/v404m_kv4_genes.csv). This anchors the whole-thalamus cross-gene
ratios to mRNA and leaves the per-type contrasts exactly as the 10x data give them. Assumes the bulk signal
for these genes is neuron-dominated.

Usage: python3 thalamus/calibrate_th10x.py [TABLES_DIR] [OUT_DIR]
Writes tables/th10x_calibration_factors.csv, th10x_calibrated_loop.csv, th10x_calibrated_nuclei.csv,
       th10x_rt_clusters.csv, figures/th10x_kv4_calibrated.png
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)
KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']; KCHIP = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']; DPP = ['Dpp6', 'Dpp10']
GENES = KV4 + KCHIP + DPP

L = pd.read_csv(f'{TABLES}/th10x_loop_summary.csv', index_col=0)
N = pd.read_csv(f'{TABLES}/th10x_nucleus_summary.csv', index_col=0)
C = pd.read_csv(f'{TABLES}/th10x_cluster_means.csv')
# ---- anchor: ThalamoSeq full-length relay-neuron TPM (neurons only, 22 nuclei) when available, else bulk thalamus exon CPM.
# The 10x side of the anchor is the cell-weighted mean over the relay loop groups only (TRN/ZI/epithalamus excluded), so the
# factor compares like with like; TRN and the other non-relay groups are then scaled with the same per-gene factors.
RELAY = ['first-order sensory relay', 'higher-order relay', 'motor relay', 'intralaminar / midline', 'anterior / limbic']
ts_path = f'{TABLES}/thalamoseq_gene_level_by_nucleus.csv'
if os.path.exists(ts_path):
    T = pd.read_csv(ts_path, index_col=0)
    anchor = T[GENES].mean()                                  # mean TPM over nuclei (each nucleus = mean of its pools)
    anchor_name = 'thalamoseq_relay_tpm'
    L2 = L.loc[[g for g in RELAY if g in L.index]]
else:
    b = pd.read_csv(f'{TABLES}/v404m_kv4_genes.csv'); b = b[b.tissue == 'thalamus'].set_index('gene')
    anchor = 2.0 ** b.mean_logCPM_wt; anchor_name = 'bulk_exon_cpm_wt'
    L2 = L[L.index != 'other']
w = L2.n_cells / L2.n_cells.sum()
tenx = pd.Series({g: (L2[f'cpm_{g}'] * w).sum() for g in GENES})
cal = anchor[GENES] / tenx
pd.DataFrame({anchor_name: anchor[GENES], 'tenx_relay_mean_cpm': tenx, 'factor': cal}).to_csv(f'{OUT}/tables/th10x_calibration_factors.csv')
print('calibration anchor:', anchor_name)


def calibrate(df):
    out = pd.DataFrame({g: df[f'cpm_{g}'] * cal[g] for g in GENES}, index=df.index)
    kv = out[KV4].sum(axis=1); kc = out[KCHIP].sum(axis=1); dp = out[DPP].sum(axis=1)
    res = pd.DataFrame({'n_cells': df.n_cells.astype(int)})
    for g in GENES:
        res[f'{g}_cal_cpm'] = out[g]
    for g in KV4:
        res[f'share_{g}'] = out[g] / kv
    for g in KCHIP:
        res[f'KChIP_comp_{g}'] = out[g] / kc
    for g in DPP:
        res[f'DPP_comp_{g}'] = out[g] / dp
    res['KChIP_per_Kv4'] = kc / kv; res['DPP_per_Kv4'] = dp / kv
    # raw within-gene contrasts that need no calibration
    res['Kcnd3_raw_cpm'] = df['cpm_Kcnd3']; res['Kcnd2_raw_cpm'] = df['cpm_Kcnd2']
    res['det_Kcnd3'] = df['det_Kcnd3']
    for g in ['Camk2a', 'Slc32a1', 'Prox1']:
        if f'det_{g}' in df:
            res[f'det_{g}'] = df[f'det_{g}']
    return res


CL = calibrate(L); CL.to_csv(f'{OUT}/tables/th10x_calibrated_loop.csv')
CN = calibrate(N[N.n_cells >= 300]).sort_values('share_Kcnd2', ascending=False); CN.to_csv(f'{OUT}/tables/th10x_calibrated_nuclei.csv')
rt = C[C.subclass.str.contains('RT-ZI') & (C.n_cells >= 50)].copy()
rt.to_csv(f'{OUT}/tables/th10x_rt_clusters.csv', index=False)

# ---------------- figure
fig, axes = plt.subplots(1, 3, figsize=(17, 7.2), gridspec_kw={'width_ratios': [1.1, 1.2, 1.2]})
order = [i for i in CL.index if i != 'other']
y = np.arange(len(order))
cols = {'Kcnd1': '#4FA8D8', 'Kcnd2': '#E8822D', 'Kcnd3': '#B9A3D9'}
ax = axes[0]; left = np.zeros(len(order))
for g in KV4:
    v = CL.loc[order, f'share_{g}'].values
    ax.barh(y, v, left=left, color=cols[g], label=g.replace('Kcnd', 'Kv4.'), edgecolor='k', lw=0.3); left += v
for i, o in enumerate(order):
    ax.text(1.01, i, f'n={CL.loc[o, "n_cells"]:,}', va='center', fontsize=7)
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=9); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_xlabel('share of Kv4 mRNA (calibrated)'); ax.set_title('A  Kv4 paralog share by loop element', loc='left', fontsize=10)
ax.legend(fontsize=12, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3)
kc = {'Kcnip1': '#7B4FA3', 'Kcnip2': '#4FA8D8', 'Kcnip3': '#4CAF50', 'Kcnip4': '#E8822D'}
ax = axes[1]; left = np.zeros(len(order))
for g in KCHIP:
    v = CL.loc[order, f'KChIP_comp_{g}'].values
    ax.barh(y, v, left=left, color=kc[g], label=g.replace('Kcnip', 'KChIP'), edgecolor='k', lw=0.3); left += v
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_xlabel('share of KChIP mRNA (calibrated)'); ax.set_title('B  KChIP partner', loc='left', fontsize=10); ax.legend(fontsize=12, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=4)
ax = axes[2]; left = np.zeros(len(order))
for g, c in zip(DPP, ['#7B4FA3', '#E8822D']):
    v = CL.loc[order, f'DPP_comp_{g}'].values
    ax.barh(y, v, left=left, color=c, label=g.replace('Dpp', 'DPP'), edgecolor='k', lw=0.3); left += v
for i, o in enumerate(order):
    ax.text(1.01, i, f'KChIP/Kv4 {CL.loc[o, "KChIP_per_Kv4"]:.1f}, DPP/Kv4 {CL.loc[o, "DPP_per_Kv4"]:.1f}', va='center', fontsize=7)
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_xlabel('share of DPP mRNA (calibrated)'); ax.set_title('C  DPP partner and subunit:channel ratios', loc='left', fontsize=10); ax.legend(fontsize=12, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=2)
fig.suptitle('Kv4 complex across the thalamocortical loop — ABC Atlas 10x (183,553 thalamic neurons), cross-gene ratios calibrated to ThalamoSeq full-length relay-neuron TPM', y=1.0)
fig.savefig(f'{OUT}/figures/th10x_kv4_calibrated.png', dpi=160, bbox_inches='tight')
pd.set_option('display.width', 220)
print(CL[['n_cells', 'share_Kcnd2', 'share_Kcnd3', 'KChIP_comp_Kcnip1', 'KChIP_comp_Kcnip2', 'KChIP_comp_Kcnip3', 'KChIP_comp_Kcnip4', 'DPP_comp_Dpp6', 'KChIP_per_Kv4', 'DPP_per_Kv4']].round(2).to_string())
