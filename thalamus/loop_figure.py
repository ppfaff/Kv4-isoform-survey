#!/usr/bin/env python3
"""The Kv4 complex around the thalamocortical loop: cortex -> TRN -> relay, one row per cell population.

Rows and sources (all gene-level TPM from full-length data except TRN/ZI, which are 10x calibrated to ThalamoSeq):
  cortex   Tasic 2016 V1 per-type TPM (tables/ctx_loop_types_gene_tpm.csv): L2/3 IT, L4, L5 IT, L5 ET, L6 CT, L6b,
           Pvalb, Sst, Vip, Lamp5/Ndnf;  Allen SSv4 adds CA1 and DG for reference
  TRN      ABC Atlas 10x RT-ZI cells, Spp1+ core vs rest (tables/th10x_trn_core_vs_shell.csv), calibrated per gene with
           the ThalamoSeq-anchored factors (tables/th10x_calibration_factors.csv); ZI/LGv likewise from th10x_calibrated_loop.csv
  relay    ThalamoSeq pools grouped by loop element (tables/thalamoseq_gene_level_by_nucleus.csv)

Panels: A Kv4 paralog share; B KChIP partner; C DPP partner; D auxiliary:Kv4 ratios; E V404M exposure index =
Kv4.2 share of Kv4 mRNA (a heterozygous dominant-negative subunit poisons tetramers in proportion to Kv4.2's share),
with Kv4 abundance (total Kv4 TPM, log scale) as the second axis. Rows ordered around the loop.
Usage: python3 thalamus/loop_figure.py [TABLES_DIR] [OUT_DIR]
Writes figures/kv4_thalamocortical_loop.png and tables/kv4_thalamocortical_loop.csv
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/figures', exist_ok=True); os.makedirs(f'{OUT}/tables', exist_ok=True)
KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']; KCHIP = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']; DPP = ['Dpp6', 'Dpp10']; GENES = KV4 + KCHIP + DPP

rows = []
# ---- cortex (Tasic) + hippocampus reference (Allen)
C = pd.read_csv(f'{TABLES}/ctx_loop_types_gene_tpm.csv')
CTX_ORDER = [('tasic2016_V1', 'L2/3 IT', 'L2/3 IT'), ('tasic2016_V1', 'L4', 'L4'), ('tasic2016_V1', 'L5 IT (L5a)', 'L5 IT'),
             ('tasic2016_V1', 'L5 ET (L5b, pyramidal tract)', 'L5 ET (→ HO thalamus)'), ('tasic2016_V1', 'L6 CT (corticothalamic)', 'L6 CT (→ relay + TRN)'),
             ('tasic2016_V1', 'L6b', 'L6b'), ('tasic2016_V1', 'Pvalb interneurons', 'Pvalb IN'), ('tasic2016_V1', 'Sst interneurons', 'Sst IN'),
             ('tasic2016_V1', 'Vip interneurons', 'Vip IN'), ('tasic2016_V1', 'Lamp5/Ndnf interneurons', 'Lamp5/Ndnf IN'),
             ('allen_ssv4', 'CA1-ProS', 'CA1 (ref.)'), ('allen_ssv4', 'DG', 'DG granule (ref.)')]
for ds, grp, lab in CTX_ORDER:
    r = C[(C.dataset == ds) & (C.group == grp)]
    if len(r) == 0:
        continue
    r = r.iloc[0]
    rows.append({'section': 'cortex' if ds == 'tasic2016_V1' else 'hippocampus (ref.)', 'population': lab, 'source': 'Tasic 2016 SMART-seq' if ds == 'tasic2016_V1' else 'Allen SSv4',
                 'n': int(r.n_cells), **{g: r[f'tpm_{g}'] for g in GENES}, 'det_Kcnd2': r['det_Kcnd2']})
# ---- TRN core / shell and ZI (10x, calibrated)
cal = pd.read_csv(f'{TABLES}/th10x_calibration_factors.csv', index_col=0)['factor']
T = pd.read_csv(f'{TABLES}/th10x_trn_core_vs_shell.csv')
for _, r in T.iterrows():
    lab = 'TRN core (Spp1+, → FO relay)' if r['class'].startswith('Spp1') else 'TRN shell / rest (Ecel1+, → HO relay)'
    rows.append({'section': 'TRN', 'population': lab, 'source': 'ABC Atlas 10x, calibrated', 'n': int(r.n),
                 **{g: r[f'cpm_{g}'] * cal[g] for g in GENES}, 'det_Kcnd2': r['det_Kcnd2']})
L10 = pd.read_csv(f'{TABLES}/th10x_calibrated_loop.csv', index_col=0)
zi = L10.loc['ZI / LGv / IGL (non-relay)']
rows.append({'section': 'TRN', 'population': 'ZI / LGv / IGL (ref.)', 'source': 'ABC Atlas 10x, calibrated', 'n': int(zi.n_cells),
             **{g: zi[f'{g}_cal_cpm'] for g in GENES}, 'det_Kcnd2': np.nan})
# ---- relay (ThalamoSeq), grouped by loop element
S = pd.read_csv(f'{TABLES}/thalamoseq_gene_level_by_nucleus.csv', index_col=0)
RELAY = [('first-order sensory relay', 'FO sensory relay (LGd, VB, VPMpc)'), ('motor relay', 'motor relay (VA, VL, VM)'),
         ('anterior / limbic', 'anterior / limbic (AD, AV, AM, LD, PVT)'), ('higher-order relay', 'HO relay (PO, LP, MD)'),
         ('intralaminar / midline', 'intralaminar / midline (CM, CL, PCN, PF, IMD, IAD, RE)')]
for key, lab in RELAY:
    sub = S[S.loop_group == key]
    rows.append({'section': 'relay', 'population': lab, 'source': f'ThalamoSeq, {len(sub)} nuclei', 'n': int(sub.n_pools.sum()),
                 **{g: sub[g].mean() for g in GENES}, 'det_Kcnd2': np.nan})
D = pd.DataFrame(rows)
kv = D[KV4].sum(axis=1); kc = D[KCHIP].sum(axis=1); dp = D[DPP].sum(axis=1)
for g in KV4:
    D[f'share_{g}'] = D[g] / kv
for g in KCHIP:
    D[f'KChIP_{g}'] = D[g] / kc
for g in DPP:
    D[f'DPP_{g}'] = D[g] / dp
D['KChIP_per_Kv4'] = kc / kv; D['DPP_per_Kv4'] = dp / kv; D['Kv4_total'] = kv
D['exposure'] = D['share_Kcnd2']            # fraction of Kv4 subunits that are Kv4.2 = fraction of tetramer slots a V404M allele can occupy
D.to_csv(f'{OUT}/tables/kv4_thalamocortical_loop.csv', index=False)

# ---------------- figure
order = list(D.index)
y = np.arange(len(order)); n = len(order)
fig, axes = plt.subplots(1, 5, figsize=(24, 0.5 * n + 3), gridspec_kw={'width_ratios': [1.1, 1.1, 0.9, 0.8, 1.0], 'wspace': 0.08})
sec_col = {'cortex': '#FFF4E5', 'hippocampus (ref.)': '#F2F2F2', 'TRN': '#FDE9E7', 'relay': '#E8F1FB'}
for ax in axes:
    for i, s in enumerate(D.section):
        ax.axhspan(i - 0.5, i + 0.5, color=sec_col[s], zorder=0)
ax = axes[0]; left = np.zeros(n)
for g, c, lab in [('Kcnd1', '#4FA8D8', 'Kv4.1'), ('Kcnd2', '#E8822D', 'Kv4.2'), ('Kcnd3', '#B9A3D9', 'Kv4.3')]:
    v = D[f'share_{g}'].values; ax.barh(y, v, left=left, color=c, edgecolor='k', lw=0.3, label=lab); left += v
for i in range(n):
    ax.text(D.share_Kcnd2.iloc[i] / 2 + D.share_Kcnd1.iloc[i], i, f'{D.share_Kcnd2.iloc[i]:.2f}', ha='center', va='center', fontsize=9, color='w', fontweight='bold')
ax.set_yticks(y); ax.set_yticklabels([f'{p}   (n={m:,})' for p, m in zip(D.population, D.n)], fontsize=11); ax.invert_yaxis(); ax.set_xlim(0, 1)
ax.set_title('A  Kv4 subunit share (number = Kv4.2 share)', loc='left', fontsize=11); ax.legend(fontsize=11, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=3)
ax = axes[1]; left = np.zeros(n)
for g, c, lab in [('Kcnip1', '#7B4FA3', 'KChIP1'), ('Kcnip2', '#4FA8D8', 'KChIP2'), ('Kcnip3', '#4CAF50', 'KChIP3'), ('Kcnip4', '#E8822D', 'KChIP4')]:
    v = D[f'KChIP_{g}'].values; ax.barh(y, v, left=left, color=c, edgecolor='k', lw=0.3, label=lab); left += v
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_title('B  KChIP partner', loc='left', fontsize=11)
ax.legend(fontsize=11, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=4)
ax = axes[2]; left = np.zeros(n)
for g, c, lab in [('Dpp6', '#7B4FA3', 'DPP6'), ('Dpp10', '#E8822D', 'DPP10')]:
    v = D[f'DPP_{g}'].values; ax.barh(y, v, left=left, color=c, edgecolor='k', lw=0.3, label=lab); left += v
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_title('C  DPP partner', loc='left', fontsize=11)
ax.legend(fontsize=11, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=2)
ax = axes[3]
ax.barh(y - 0.2, D.KChIP_per_Kv4.clip(upper=14.5), 0.38, color='#1B9E77', label='KChIP : Kv4'); ax.barh(y + 0.2, D.DPP_per_Kv4, 0.38, color='#D62728', label='DPP : Kv4')
for i in range(n):
    ax.text(min(D.KChIP_per_Kv4.iloc[i], 14.2) + 0.1, i - 0.2, f'{D.KChIP_per_Kv4.iloc[i]:.1f}' + (' ▶' if D.KChIP_per_Kv4.iloc[i] > 15 else ''), va='center', fontsize=8)
    ax.text(D.DPP_per_Kv4.iloc[i] + 0.1, i + 0.2, f'{D.DPP_per_Kv4.iloc[i]:.1f}', va='center', fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 15)
ax.set_title('D  auxiliary mRNA per Kv4 mRNA', loc='left', fontsize=11); ax.legend(fontsize=11, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.06), ncol=2)
ax = axes[4]
cols = ['#D6372E' if e >= 0.65 else '#E8822D' if e >= 0.5 else '#999999' for e in D.exposure]
ax.barh(y, D.exposure, color=cols, edgecolor='k', lw=0.3)
ax2 = ax.twiny()
ax2.plot(D.Kv4_total, y, 'k|', ms=12, mew=2); ax2.set_xscale('log'); ax2.set_xlabel('total Kv4 TPM (black ticks, log)', fontsize=9)
ax.set_yticks(y); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_xlim(0, 1); ax.axvline(0.5, color='k', lw=0.6, ls=':')
ax.set_xlabel('Kv4.2 share of Kv4 subunits'); ax.set_title('E  V404M exposure (red ≥0.65, orange ≥0.5)', loc='left', fontsize=11)
for i, s in enumerate(D.section):
    if i == 0 or D.section.iloc[i - 1] != s:
        axes[0].text(-1.02, i - 0.5, s.upper(), transform=axes[0].get_yaxis_transform(), ha='left', va='bottom', fontsize=9, color='#555', fontweight='bold')
fig.suptitle('The Kv4 complex around the thalamocortical loop: cortex (Tasic 2016 full-length) → TRN (ABC Atlas 10x, calibrated) → relay (ThalamoSeq full-length)\n'
             'A heterozygous dominant-negative Kv4.2 subunit affects a population in proportion to Kv4.2\'s share of its Kv4 pool; TRN and L6 CT have the least Kv4.3 backup', y=1.02, fontsize=12)
fig.savefig(f'{OUT}/figures/kv4_thalamocortical_loop.png', dpi=150, bbox_inches='tight')
pd.set_option('display.width', 220)
print(D[['section', 'population', 'n', 'share_Kcnd2', 'share_Kcnd3', 'KChIP_per_Kv4', 'DPP_per_Kv4', 'Kv4_total']].round(2).to_string(index=False))
