#!/usr/bin/env python3
"""Does isoform-assignment uncertainty explain the apparent per-cell variant heterogeneity?

This is the test that settles what three earlier attempts could not.

  channel_noise_floor.py        models COUNTING noise given salmon's assignment, so it is blind to assignment
                                error by construction
  variant_distribution_shape.py found DPP6/DPP10/KChIP4 variants bimodal within single cell types
  variant_exclusivity.py        tried to separate real exclusive expression from winner-take-all EM collapse
                                and CANNOT: mutual exclusivity is the expected signature of both. A cell that
                                genuinely expresses only DPP6K is indistinguishable, by that test, from one
                                whose reads the EM dumped on DPP6K. It also flagged its own positive control.

Bootstrap replicates break the tie, because they re-run the assignment on resampled reads from the SAME cell:

    bootstrap spread HIGH   the EM is not determined by the data -- replicates of one cell land on different
                            variants, so a cell's apparent "choice" of variant is arbitrary and the
                            across-cell bimodality is an artefact
    bootstrap spread LOW    the assignment is well determined within each cell, so cells that differ really
                            do differ, and the across-cell spread is biological

The statistic is the variance ratio, directly analogous to a reliability coefficient:

    ICC-like = across-cell variance / (across-cell variance + mean within-cell bootstrap variance)

    -> near 1.0  assignment is reliable; essentially all the spread is between cells (biology)
    -> near 0.0  assignment is unreliable; the spread is assignment noise
    -> ~0.5      half and half; the heterogeneity claim would need to be stated as an upper bound

Also reports, per variant, how often a cell's modal variant changes between bootstrap replicates -- the most
direct readout of collapse, and the one that needs no variance algebra to interpret.

Usage:  python3 bootstrap_vs_biological.py <bootstrap tsv> [shape csv] [out_dir]
    bootstrap tsv: concatenated tables/bootstrap_uncertainty_<HOST>.tsv from the worker machines
    shape csv:     tables/variant_distribution_shape.csv, to put the two spreads side by side
"""
import sys, os
import numpy as np, pandas as pd

BOOT  = sys.argv[1] if len(sys.argv) > 1 else 'tables/bootstrap_uncertainty.tsv'
SHAPE = sys.argv[2] if len(sys.argv) > 2 else 'tables/variant_distribution_shape.csv'
OUT   = sys.argv[3] if len(sys.argv) > 3 else '.'

KEY = ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'KChIP4a', 'Kv4.2', 'Kv4.3S', 'Kv4.3L']


def main():
    B = pd.read_csv(BOOT, sep='\t')
    need = {'cell', 'panel', 'variant', 'replicate', 'frac'}
    if not need.issubset(B.columns):
        sys.exit(f'{BOOT} missing columns {need - set(B.columns)}')
    print(f'{B.cell.nunique()} cells, {B.replicate.nunique()} replicates, '
          f'{B.variant.nunique()} variants')

    # --- per cell x variant: the bootstrap (within-cell) spread, and the cell's point estimate
    g = B.groupby(['cell', 'panel', 'variant'])['frac']
    per_cell = g.agg(boot_mean='mean', boot_sd='std').reset_index()

    # --- modal-variant instability: does the cell's top variant change between replicates?
    inst = []
    for (cell, panel), sub in B.groupby(['cell', 'panel']):
        top = sub.loc[sub.groupby('replicate')['frac'].idxmax()]
        counts = top['variant'].value_counts(normalize=True)
        inst.append(dict(cell=cell, panel=panel, modal_variant=counts.index[0],
                         modal_stability=float(counts.iloc[0]), n_distinct_top=int(len(counts))))
    I = pd.DataFrame(inst)

    # --- variance decomposition per panel x variant, pooled over cells
    rows = []
    for (panel, variant), sub in per_cell.groupby(['panel', 'variant']):
        if len(sub) < 10:
            continue
        across = sub['boot_mean'].var(ddof=1)                 # between-cell variance of the point estimate
        within = np.nanmean(sub['boot_sd'] ** 2)              # mean within-cell assignment variance
        icc = across / (across + within) if (across + within) > 0 else np.nan
        rows.append(dict(panel=panel, variant=variant, n_cells=len(sub),
                         across_cell_sd=np.sqrt(across), mean_boot_sd=float(np.sqrt(within)),
                         icc=icc,
                         verdict=('biological' if icc > 0.8 else
                                  'ASSIGNMENT NOISE' if icc < 0.4 else 'mixed - state as upper bound')))
    V = pd.DataFrame(rows).sort_values('icc')

    print('\n' + '=' * 100)
    print('ASSIGNMENT RELIABILITY PER VARIANT')
    print('  across_cell_sd : SD of the per-cell point estimate between cells')
    print('  mean_boot_sd   : mean SD across bootstrap replicates WITHIN a cell (assignment uncertainty)')
    print('  icc            : across / (across + within). near 1 = reliable, near 0 = assignment noise')
    print('=' * 100)
    if not V.empty:
        print(f'{"panel":6} {"variant":10} {"n":>5} {"acrossSD":>9} {"bootSD":>8} {"icc":>6}  verdict')
        for _, r in V.iterrows():
            flag = ' <<<' if r.variant in KEY else ''
            print(f'{r.panel:6} {r.variant:10} {int(r.n_cells):5d} {r.across_cell_sd:9.4f} '
                  f'{r.mean_boot_sd:8.4f} {r.icc:6.3f}  {r.verdict}{flag}')

    print('\n' + '=' * 100)
    print('MODAL-VARIANT STABILITY  (fraction of replicates agreeing on the cell\'s top variant)')
    print('  1.00 everywhere = assignment is determined;  near 1/k = the EM is choosing at random')
    print('=' * 100)
    if not I.empty:
        s = I.groupby('panel')['modal_stability'].agg(['count', 'mean', 'median', 'min'])
        s['frac_cells_unstable'] = I.groupby('panel')['modal_stability'].apply(lambda x: (x < 0.9).mean())
        print(s.round(3).to_string())

    # --- side by side with the across-cell excess SD, if available
    if os.path.exists(SHAPE):
        S = pd.read_csv(SHAPE)
        if 'excess_sd' in S.columns:
            M = (S[S.variant.isin(KEY)]
                 .groupby(['panel', 'variant'])
                 .agg(shape_excess_sd=('excess_sd', 'median'),
                      shape_noise_sd=('median_noise_sd', 'median')).reset_index())
            J = V.merge(M, on=['panel', 'variant'], how='left')
            J['excess_over_boot'] = J['shape_excess_sd'] / J['mean_boot_sd']
            print('\n' + '=' * 100)
            print('THE COMPARISON THAT MATTERS')
            print('  shape_excess_sd : across-cell spread after removing COUNTING noise')
            print('  mean_boot_sd    : ASSIGNMENT uncertainty, which that analysis could not see')
            print('  ratio > ~3 -> the heterogeneity survives; ratio ~1 -> it was assignment noise')
            print('=' * 100)
            cols = ['panel', 'variant', 'n_cells', 'shape_excess_sd', 'shape_noise_sd',
                    'mean_boot_sd', 'excess_over_boot', 'icc', 'verdict']
            print(J[J.shape_excess_sd.notna()][cols].round(4).to_string(index=False))

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    V.to_csv(f'{OUT}/tables/bootstrap_reliability.csv', index=False)
    I.to_csv(f'{OUT}/tables/bootstrap_modal_stability.csv', index=False)
    per_cell.to_csv(f'{OUT}/tables/bootstrap_per_cell.csv', index=False)
    print(f'\nwrote {OUT}/tables/bootstrap_{{reliability,modal_stability,per_cell}}.csv')


if __name__ == '__main__':
    main()
