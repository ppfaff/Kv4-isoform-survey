#!/usr/bin/env python3
"""TRN subnetworks of Hartley et al. (2024, Cell Rep 43:115037): Spp1+ (core; first-order-projecting) vs
Ecel1+ (shell; higher-order-projecting) reticular neurons, classified per cell in the ABC Atlas 10x data.

Spp1 mRNA is sparse in 10x data (detected in <1% of cells in most RT clusters), so the split is done per cell,
not per cluster: a cell is 'Spp1+' if Spp1 > 0, 'Ecel1+' if Ecel1 > 0 and Spp1 == 0, 'double' if both,
'neither' otherwise. Also classes by Pvalb/Sst (high = above the RT median) for comparison.

Runs where th10x_cells.parquet lives (Lambda). Usage: python3 trn_subnetworks.py TABLES_DIR OUT_DIR
Writes OUT/tables/th10x_trn_subnetworks.csv (per class: n, Kv4-complex means, detection, Pvalb/Sst) and
       OUT/tables/th10x_trn_cells.csv.gz (per RT cell: class labels + panel genes; ~9k rows) for plotting here.
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from th_common import KV4, KCHIP, DPP, read_cells  # noqa: E402

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True)
df = read_cells(f'{TABLES}/th10x_cells')
rt = df[df['subclass'].str.contains('RT-ZI', na=False)].copy()
rt = rt[rt['Gad2'] > 0]                                   # GABAergic RT-ZI cells only
spp = rt['Spp1'] > 0; ece = rt['Ecel1'] > 0
rt['subnetwork'] = np.select([spp & ~ece, ece & ~spp, spp & ece], ['Spp1+ (core/FO)', 'Ecel1+ (shell/HO)', 'Spp1+Ecel1+'], 'neither')
rt['pv_sst'] = np.select([(rt.Pvalb > rt.Pvalb.median()) & (rt.Sst <= rt.Sst.median()),
                          (rt.Sst > rt.Sst.median()) & (rt.Pvalb <= rt.Pvalb.median()),
                          (rt.Pvalb > rt.Pvalb.median()) & (rt.Sst > rt.Sst.median())],
                         ['Pvalb-high', 'Sst-high', 'Pvalb+Sst-high'], 'both-low')
genes = [g for g in KV4 + KCHIP + DPP + ['Pvalb', 'Sst', 'Spp1', 'Ecel1', 'Gnb3', 'Camk2a', 'Slc32a1', 'Cacna1g', 'Cacna1h', 'Cacna1i', 'Hcn2', 'Kcnc1', 'Kcnc2', 'Kcnn2'] if g in rt.columns]
keep = ['cluster', 'supertype', 'feature_matrix_label', 'subnetwork', 'pv_sst'] + genes
rt[keep].to_csv(f'{OUT}/tables/th10x_trn_cells.csv.gz')

rows = []
for key in ['subnetwork', 'pv_sst']:
    for lab, sub in rt.groupby(key):
        L = 2.0 ** sub[genes] - 1.0
        r = {'grouping': key, 'class': lab, 'n_cells': len(sub)}
        for g in genes:
            r[f'cpm_{g}'] = L[g].mean(); r[f'det_{g}'] = (sub[g] > 0).mean()
        kv = sum(r[f'cpm_{g}'] for g in KV4)
        r['share_Kcnd2_raw'] = r['cpm_Kcnd2'] / kv; r['Kcnd3_over_Kcnd2_raw'] = r['cpm_Kcnd3'] / r['cpm_Kcnd2']
        r['Kcnip1_over_Kcnip2_raw'] = r['cpm_Kcnip1'] / max(r['cpm_Kcnip2'], 1e-9)
        r['Dpp10_over_Dpp6_raw'] = r['cpm_Dpp10'] / r['cpm_Dpp6']
        rows.append(r)
S = pd.DataFrame(rows); S.to_csv(f'{OUT}/tables/th10x_trn_subnetworks.csv', index=False)
pd.set_option('display.width', 250)
print(f'{len(rt):,} GABAergic RT-ZI cells')
print(S[['grouping', 'class', 'n_cells', 'cpm_Kcnd2', 'cpm_Kcnd3', 'Kcnd3_over_Kcnd2_raw', 'cpm_Kcnip1', 'cpm_Kcnip2', 'cpm_Kcnip4',
         'cpm_Dpp6', 'cpm_Dpp10', 'cpm_Pvalb', 'cpm_Sst', 'cpm_Cacna1g', 'cpm_Cacna1h', 'cpm_Cacna1i']].round(1).to_string(index=False))
print('\nsubnetwork x Pvalb/Sst cross-tab:'); print(pd.crosstab(rt.subnetwork, rt.pv_sst))
