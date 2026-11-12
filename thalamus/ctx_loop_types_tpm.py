#!/usr/bin/env python3
"""Per-cell-type Kv4-complex gene TPM and coding-variant fractions for the cortical side of the thalamocortical loop,
from the full-length per-cell tables (Tasic 2016 V1 and Allen SSv4 cortex/HPF). Runs where those tables are
(CTX_DATA / HIP_DATA env vars or tables/*_isoform_per_cell_annotated.csv, as in scripts/hip_common.py).

Loop-relevant groups (pooled over the listed types; TPM averaged over cells, variant fractions reads-weighted):
  Tasic 2016 V1:   L6 CT (L6a Car12/Mgp/Sla/Syt17), L6b, L5 ET (L5b Cdh13/Chrna6/Tph2), L5 IT (L5a *), L4, L2/3,
                   Pvalb, Sst, Vip, Lamp5/Ndnf
  Allen SSv4:      L6 CT CTX, L6b CTX, L6 IT CTX, L4/5 IT CTX, L4 RSP-ACA, Pvalb, Sst, Vip, Lamp5, CA1-ProS, DG
Writes tables/ctx_loop_types_gene_tpm.csv  (dataset, group, n_cells, tpm_<gene> ... ) and
       tables/ctx_loop_types_variant_fracs.csv (dataset, group, gene, vclass, frac, reads)
Usage: python3 thalamus/ctx_loop_types_tpm.py [OUT_DIR]
"""
import os, sys
import pandas as pd, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
from hip_common import load_ctx, load_hip  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True)
GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']

CTX_GROUPS = {
    'L6 CT (corticothalamic)': lambda s: s.startswith('L6a'),
    'L6b': lambda s: s.startswith('L6b'),
    'L5 ET (L5b, pyramidal tract)': lambda s: s.startswith('L5b'),
    'L5 IT (L5a)': lambda s: s.startswith('L5a') or s == 'L5 Ucma',
    'L4': lambda s: s.startswith('L4'),
    'L2/3 IT': lambda s: s.startswith('L2'),
    'Pvalb interneurons': lambda s: s.startswith('Pvalb'),
    'Sst interneurons': lambda s: s.startswith('Sst'),
    'Vip interneurons': lambda s: s.startswith('Vip'),
    'Lamp5/Ndnf interneurons': lambda s: s.startswith('Ndnf') or s == 'Smad3',
}
HIP_GROUPS = {
    'L6 CT CTX': lambda s: s == 'L6 CT CTX', 'L6b CTX': lambda s: s == 'L6b CTX', 'L6 IT CTX': lambda s: s == 'L6 IT CTX',
    'L4/5 IT CTX': lambda s: s == 'L4/5 IT CTX', 'L4 RSP-ACA': lambda s: s == 'L4 RSP-ACA',
    'Pvalb': lambda s: s == 'Pvalb', 'Sst': lambda s: s in ('Sst', 'Sst Chodl'), 'Vip': lambda s: s == 'Vip', 'Lamp5': lambda s: s == 'Lamp5',
    'CA1-ProS': lambda s: s == 'CA1-ProS', 'DG': lambda s: s == 'DG',
}


def summarize(d, groups, dataset):
    rows, vrows = [], []
    for grp, sel in groups.items():
        sub = d[d.subclass.map(sel)]
        if sub.sample_name.nunique() < 8:
            continue
        cells = sub.sample_name.unique()
        # per-cell gene TPM (sum over transcripts), averaged over cells; cells lacking a gene count as 0
        gt = sub.groupby(['sample_name', 'gene']).TPM.sum().unstack('gene').reindex(index=cells, columns=GENES).fillna(0.0)
        r = {'dataset': dataset, 'group': grp, 'n_cells': len(cells)}
        for g in GENES:
            r[f'tpm_{g}'] = gt[g].mean(); r[f'det_{g}'] = (gt[g] > 0).mean()
        rows.append(r)
        v = sub.groupby(['gene', 'vclass']).reads.sum().reset_index()
        v['frac'] = v.reads / v.groupby('gene').reads.transform('sum')
        v['dataset'] = dataset; v['group'] = grp
        vrows.append(v)
    return pd.DataFrame(rows), pd.concat(vrows)


out, vout = [], []
for loader, groups, name in [(load_ctx, CTX_GROUPS, 'tasic2016_V1'), (load_hip, HIP_GROUPS, 'allen_ssv4')]:
    try:
        d = loader()
    except FileNotFoundError as e:
        print(f'skipping {name}: {e}', file=sys.stderr); continue
    if 'TPM' not in d.columns:
        sys.exit(f'{name}: per-cell table has no TPM column; re-run the collapse with TPM retained')
    a, b = summarize(d, groups, name); out.append(a); vout.append(b)
G = pd.concat(out); V = pd.concat(vout)
G.to_csv(f'{OUT}/tables/ctx_loop_types_gene_tpm.csv', index=False)
V.to_csv(f'{OUT}/tables/ctx_loop_types_variant_fracs.csv', index=False)
pd.set_option('display.width', 220)
print(G[['dataset', 'group', 'n_cells'] + [f'tpm_{g}' for g in GENES]].round(1).to_string(index=False))
