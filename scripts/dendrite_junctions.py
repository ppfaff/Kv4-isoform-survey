#!/usr/bin/env python3
"""Junction-level test of the two variant shifts that survived the cell-type-mixing check in the Farris laminar data.

Salmon's isoform estimates for Kcnd3 are confounded by the 4.8-kb 3'UTR that only the Kv4.3S model (Kcnd3-202)
carries, so the Kv4.3L vs S call is re-done here from STAR splice-junction counts (SJ.out.tab), which depend only on
reads crossing the L-specific cassette exon. The Kcnip3 cys-rich vs KChIP3a difference is a 4-nt alternative acceptor
on exon 3 and is read the same way.

Inputs (made on Lambda from star/*/SJ.out.tab; see README):
  tables/sj_kcnd3_kcnip3.tsv   run, chr, start, end, strand, motif, annotated, unique, multi, overhang (STAR columns)
  tables/dendrite_samples.tsv
Outputs: tables/dendrite_junction_psi.csv, figures/dendrite_junction_psi.png
Usage:   python3 scripts/dendrite_junctions.py [TABLES_DIR] [OUT_DIR]
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

# GRCm39 / GENCODE vM36 coordinates (STAR SJ.out.tab gives intron start/end, 1-based)
EVENTS = {
    'Kv4.3L exon inclusion': {           # Kcnd3-201/-203 (L) include chr3:105574278-105574334; -202 (S) skips it
        'inc': [(105572863, 105574277), (105574335, 105575395)], 'skip': [(105572863, 105575395)], 'inc_avg': True},
    'Kcnip3 cys-rich acceptor': {        # Kcnip3-201 (cys-rich) exon 3 ends ...894; -203 (KChIP3a) ends ...898 (minus strand)
        'inc': [(127307895, 127352748)], 'skip': [(127307899, 127352748)], 'inc_avg': False},
}

sj = pd.read_csv(f'{TABLES}/sj_kcnd3_kcnip3.tsv', sep='\t', header=None,
                 names=['run', 'chr', 'start', 'end', 'strand', 'motif', 'annot', 'uniq', 'multi', 'overhang'])
S = pd.read_csv(f'{TABLES}/dendrite_samples.tsv', sep='\t', dtype=str)
if 'extra' in S.columns:
    S = S[S.extra != '1']
S = S[S.dataset == 'farris']
meta = S.set_index('run_accession')[['region', 'compartment', 'unit']]


def junc(s, e):
    x = sj[(sj.start == s) & (sj.end == e)].set_index('run').uniq
    return x.reindex(meta.index).fillna(0)


rows, per_sample = [], []
for name, ev in EVENTS.items():
    inc = sum(junc(*j) for j in ev['inc']) / (len(ev['inc']) if ev['inc_avg'] else 1)
    skip = sum(junc(*j) for j in ev['skip'])
    t = meta.copy(); t['inc'] = inc.values; t['skip'] = skip.values; t['psi'] = t.inc / (t.inc + t.skip)
    t['event'] = name; per_sample.append(t.reset_index())
    p = t.pivot_table(index=['region', 'unit'], columns='compartment', values='psi').dropna()
    d = p.dendrite - p.soma
    try:
        pv = stats.wilcoxon(p.dendrite, p.soma).pvalue
    except ValueError:
        pv = np.nan
    rows.append((name, 'all', p.soma.mean(), p.dendrite.mean(), d.mean(), (d > 0).sum(), len(d), pv))
    for rg, g in t.groupby('region'):
        tot = g.groupby('compartment')[['inc', 'skip']].sum()
        rows.append((name, rg, tot.loc['soma', 'inc'] / tot.loc['soma'].sum(), tot.loc['dendrite', 'inc'] / tot.loc['dendrite'].sum(),
                     np.nan, np.nan, 3, np.nan))
R = pd.DataFrame(rows, columns=['event', 'region', 'psi_soma_layer', 'psi_dendrite_layer', 'mean_paired_diff', 'n_pairs_up', 'n_pairs', 'wilcoxon_p'])
R.to_csv(f'{OUT}/tables/dendrite_junction_psi.csv', index=False)
P = pd.concat(per_sample); P.to_csv(f'{OUT}/tables/dendrite_junction_psi_per_sample.csv', index=False)
print(R.round(3).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
for ax, (name, g) in zip(axes, P.groupby('event', sort=False)):
    for i, rg in enumerate(['CA1', 'CA2', 'CA3', 'DG']):
        for unit, u in g[g.region == rg].groupby('unit'):
            s = u.set_index('compartment').psi
            ax.plot([i - 0.15, i + 0.15], [s.get('soma', np.nan), s.get('dendrite', np.nan)], '-o', color='#555', ms=3, lw=0.8, alpha=0.8)
    ax.set_xticks(range(4)); ax.set_xticklabels(['CA1', 'CA2', 'CA3', 'DG']); ax.set_ylim(0, 1)
    ax.set_title(name, fontsize=10); ax.set_ylabel('fraction of junction reads')
    ax.text(0.02, 0.95, 'left = cell-body layer, right = dendritic layer', transform=ax.transAxes, fontsize=7, va='top', color='#555')
fig.suptitle('Farris 2019: splice-junction usage by hippocampal layer (STAR unique junction reads)', fontsize=10)
fig.tight_layout(); fig.savefig(f'{OUT}/figures/dendrite_junction_psi.png', dpi=200)
print('done')
