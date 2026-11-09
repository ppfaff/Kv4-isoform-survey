#!/usr/bin/env python3
"""Empirical calibration of the V404M bulk variant calls against full-length single-cell data.

Pools the Tasic 2016 (GSE71585, V1 neurons, full-length SMART-seq) per-type variant reads into a cortical
pseudo-bulk and sets it beside the wild-type cortex fractions from the 3'-biased polyA bulk libraries.
Variants that the bulk pipeline under-calls by more than a factor of two against the pseudo-bulk are flagged
as collapsed; those within 20% relative are flagged as reliable.

The disagreement is the 3' coverage collapse already shown by simulation for the DPP panel
(pipeline/simulate_variant_calibration.py): N-terminal variants distinguished only by 5' exons default to one
sink isoform per gene. This script extends the check to the KChIP panel without rerunning the simulator.

Usage: python3 scripts/bulk_vs_tasic_variants.py [TABLES_DIR] [OUT_DIR]
Writes OUT/tables/bulk_vs_tasic_variant_calibration.csv
"""
import os, sys
import numpy as np, pandas as pd

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True)

t = pd.read_csv(f'{TABLES}/ctx_variant_frac_of_gene_by_type.csv')
t = t[t.cls != 'NonN']
pool = t.groupby(['gene', 'vclass']).reads.sum(); pool = pool / pool.groupby(level=0).transform('sum')
glut = t[t.cls == 'Glut'].groupby(['gene', 'vclass']).reads.sum(); glut = glut / glut.groupby(level=0).transform('sum')
gaba = t[t.cls == 'GABA'].groupby(['gene', 'vclass']).reads.sum(); gaba = gaba / gaba.groupby(level=0).transform('sum')
b = pd.read_csv(f'{TABLES}/v404m_variant_fracs.csv')
bulk = b[b.tissue == 'cortex'].set_index(['gene', 'vclass'])

d = pd.DataFrame({'tasic_pseudobulk': pool, 'tasic_glut': glut, 'tasic_gaba': gaba,
                  'bulk_ctx_wt': bulk.frac_wt, 'bulk_ctx_mut': bulk.frac_mut}).fillna(0.0)
d = d[d[['tasic_pseudobulk', 'bulk_ctx_wt']].max(axis=1) > 0.01].copy()
d['ratio_bulk_over_tasic'] = d.bulk_ctx_wt / d.tasic_pseudobulk.replace(0, np.nan)


def verdict(r):
    if r.tasic_pseudobulk < 0.02 and r.bulk_ctx_wt < 0.02:
        return 'minor in both'
    if r.ratio_bulk_over_tasic < 0.5:
        return 'COLLAPSED in bulk (under-called >2x)'
    if r.ratio_bulk_over_tasic > 2.0:
        return 'SINK in bulk (over-called >2x)'
    if abs(r.ratio_bulk_over_tasic - 1) <= 0.2:
        return 'reliable'
    return 'shifted (within 2x)'


d['verdict'] = d.apply(verdict, axis=1)
d = d.reset_index().sort_values(['gene', 'tasic_pseudobulk'], ascending=[True, False])
d.to_csv(f'{OUT}/tables/bulk_vs_tasic_variant_calibration.csv', index=False)
pd.set_option('display.width', 220)
print(d.round(3).to_string(index=False))
print('\nBulk variant calls usable only for:', ', '.join(d[d.verdict == 'reliable'].vclass))
