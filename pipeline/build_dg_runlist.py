#!/usr/bin/env python3
"""Pick DG granule cells for the read-level Kcnd2 analysis and look up their NeMO FASTQ URLs.

Groups come from tables/dg_cells_markers.csv (survey): Kv4-null = zero salmon reads on Kcnd1/2/3; Kv4-low / Kv4-high =
below / above the median Kcnd2 among expressing cells. All null cells are taken, plus a seeded random sample of low and
high cells of matched size (default 150 each) as controls. Run on Lambda in ~/allen_meta (where hip_filtered_runlist.tsv
with the NeMO URLs lives) or pass that file's path.

Usage:  python3 build_dg_runlist.py [hip_filtered_runlist.tsv] [dg_cells_markers.csv] [N_CONTROL]
Writes: tables/dg_readlevel_runlist.tsv (same columns as the hip runlist, + grp) for run_dg_star.sh
"""
import sys, os
import pandas as pd

runlist = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser('~/allen_meta/hip_filtered_runlist.tsv')
markers = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser('~/Kv4-isoform-survey/tables/dg_cells_markers.csv')
n_ctrl = int(sys.argv[3]) if len(sys.argv) > 3 else 150

R = pd.read_csv(runlist, sep='\t', dtype=str)
D = pd.read_csv(markers)
sel = pd.concat([D[D.grp == 'Kv4-null'],
                 D[D.grp == 'Kv4-low'].sample(n=min(n_ctrl, (D.grp == 'Kv4-low').sum()), random_state=1),
                 D[D.grp == 'Kv4-high'].sample(n=min(n_ctrl, (D.grp == 'Kv4-high').sum()), random_state=1)])
M = R.merge(sel[['sample_name', 'grp', 'Kcnd2', 'Kcnd2_intronic']], on='sample_name', how='inner')
missing = set(sel.sample_name) - set(M.sample_name)
if missing:
    print(f'WARNING: {len(missing)} selected cells not in {runlist}', file=sys.stderr)
os.makedirs('tables', exist_ok=True)
M.to_csv('tables/dg_readlevel_runlist.tsv', sep='\t', index=False)
print(M.grp.value_counts().to_string(), file=sys.stderr)
print(f'{len(M)} cells -> tables/dg_readlevel_runlist.tsv', file=sys.stderr)
