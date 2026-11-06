#!/usr/bin/env python3
"""Does the measured subunit excess quantitatively explain KChIP's short half-life?

THE THREE INPUTS, all independent of each other.

1. FUNCTIONAL (the lab's own rescue). DPP6 RNAi knocks down endogenous DPP6; co-expressing an
   RNAi-resistant rat DPP6 from a CMV promoter with a beta-globin intron restores the current - but
   only to NORMAL, not above it, despite a construct that should massively overexpress. So DPP6 is
   not the limiting component. Separately, overexpressing Kv4.2 does raise current amplitude. Taken
   together: Kv4.2 is limiting and the auxiliary subunits are in functional excess.

2. TRANSCRIPTIONAL (this project's single-cell data). Per-cell family totals give
   KChIP:Kv4 = 3.24 and DPP:Kv4 = 2.67 in cortical excitatory neurons (2.67 and 2.95 in GABAergic),
   and the per-cell `limiting` call already in ctx_acurrent_index_per_cell.csv names Kv4 the limiting
   family more often than either auxiliary family. Independent agreement with (1).

3. PROTEOMIC (Dorrbaum 2018). KChIP4 2.02 d, KChIP1 2.67 d, against Kv4.3 3.87 d in the same
   preparation.

THE TEST. If Kv4 is limiting, essentially every Kv4 subunit acquires a KChIP, so the KChIP flux that
assembles equals the Kv4 flux and the remainder is surplus. With a synthesis ratio r, a fraction 1/r
of KChIP synthesis assembles and (1 - 1/r) is surplus. Steady-state pools are flux x lifetime:

    f_free = [(1 - 1/r) * t_free] / [(1 - 1/r) * t_free + (1/r) * t_bound]

Under pure surplus quality control with NO exchange, bound KChIP lives as long as the complex, so
t_bound = t_complex. The question is then whether the resulting mixture, fitted as a single
exponential at Dorrbaum's harvest points, reproduces 2.02 d. If the predicted apparent half-life is
longer than observed, surplus alone is insufficient and something must additionally remove KChIP
from intact complexes.

This is a stronger test than the earlier free-parameter sweep, because r is now measured rather than
fitted, and it is measured two independent ways.

Usage:  python3 surplus_from_stoichiometry.py [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
LN2 = np.log(2)
T = np.array([0.0, 1.0, 3.0, 7.0])      # Dorrbaum's actual harvest timepoints
T_COMPLEX = 3.87                        # d, Kv4.3, same preparation
T_OBS = {'KChIP4': 2.02, 'KChIP1': 2.67}
R_MRNA = {'cortex Glut': 3.24, 'cortex GABA': 2.67}


def fit_half_life(y):
    k = -np.polyfit(T, np.log(y), 1)[0]
    pred = np.exp(-k * T)
    ss, st = np.sum((y - pred) ** 2), np.sum((y - y.mean()) ** 2)
    return LN2 / k, 1 - ss / st


def mixture(f_free, t_free, t_bound=T_COMPLEX):
    return f_free * np.exp(-LN2 / t_free * T) + (1 - f_free) * np.exp(-LN2 / t_bound * T)


def f_free_from_r(r, t_free, t_bound=T_COMPLEX):
    surplus, assembled = 1 - 1 / r, 1 / r
    return (surplus * t_free) / (surplus * t_free + assembled * t_bound)


def main():
    print('=' * 100)
    print('PREDICTED apparent KChIP half-life from the MEASURED synthesis excess')
    print(f'  bound KChIP assumed to live as long as the complex ({T_COMPLEX} d) - pure surplus, no exchange')
    print('=' * 100)
    rows = []
    for label, r in R_MRNA.items():
        print(f'\n{label}:  KChIP:Kv4 mRNA ratio r = {r:.2f}  '
              f'-> {(1-1/r)*100:.0f}% of KChIP synthesis is surplus')
        print(f'  {"free t1/2":>11} {"f_free":>8} {"predicted apparent t1/2":>25} {"R2":>7}')
        for h in (2, 6, 12, 24, 48):
            tf = h / 24
            f = f_free_from_r(r, tf)
            th, r2 = fit_half_life(mixture(f, tf))
            flag = ''
            if th <= T_OBS['KChIP4']:
                flag = '  <= observed KChIP4'
            elif th <= T_OBS['KChIP1']:
                flag = '  <= observed KChIP1'
            print(f'  {h:9d} h {f:8.3f} {th:22.2f} d {r2:7.3f}{flag}')
            rows.append(dict(cell_class=label, r_mrna=r, free_half_life_h=h, f_free=round(f, 3),
                             predicted_half_life_d=round(th, 2), r2=round(r2, 3)))

    print('\n' + '=' * 100)
    print('HOW BIG WOULD THE EXCESS HAVE TO BE?  Solving for r that reproduces the observation')
    print('=' * 100)
    print(f'  {"target":9} {"free t1/2":>11} {"r required":>12} {"vs measured r=3.24":>21}')
    inv = []
    for name, obs in T_OBS.items():
        for h in (6, 12, 24, 48):
            tf = h / 24
            lo, hi = 1.001, 400.0
            for _ in range(200):
                mid = (lo + hi) / 2
                th, _ = fit_half_life(mixture(f_free_from_r(mid, tf), tf))
                if th > obs:
                    lo = mid
                else:
                    hi = mid
            r_req = lo
            ok = 'WITHIN measured range' if r_req <= 3.5 else f'{r_req/3.24:.0f}x too high'
            print(f'  {name:9} {h:9d} h {r_req:12.1f} {ok:>21}')
            inv.append(dict(target=name, target_half_life_d=obs, free_half_life_h=h,
                            r_required=round(r_req, 2)))

    print('\n' + '=' * 100)
    print('VERDICT')
    print('=' * 100)
    best = min((row for row in rows if row['cell_class'] == 'cortex Glut'),
               key=lambda x: x['predicted_half_life_d'])
    print(f'  At the measured r = 3.24, the shortest apparent half-life surplus QC can produce over')
    print(f'  any plausible free-pool lifetime is {best["predicted_half_life_d"]:.2f} d '
          f'(at a {best["free_half_life_h"]} h free pool).')
    gap = best['predicted_half_life_d'] - T_OBS['KChIP4']
    if gap > 0.15:
        print(f'  Observed for KChIP4: {T_OBS["KChIP4"]:.2f} d. SURPLUS FALLS SHORT BY {gap:.2f} d.')
        print('  The measured excess is not large enough to account for the observed turnover, so')
        print('  something must additionally remove KChIP from intact complexes. The lab\'s rescue')
        print('  result, by pinning the excess at a modest value rather than leaving it free, is what')
        print('  closes the loophole that the four-timepoint curve shape could not.')
    else:
        print(f'  Observed for KChIP4: {T_OBS["KChIP4"]:.2f} d. Surplus QC is SUFFICIENT; no exchange')
        print('  need be invoked.')
    print('\n  CAVEATS THAT COULD REOPEN IT. r is an mRNA ratio standing in for a protein synthesis')
    print('  ratio, and differential translation efficiency would move it; the binding stoichiometry')
    print('  is taken as one KChIP per Kv4 subunit (4:4 per tetramer), and sub-saturated complexes')
    print('  would raise the effective surplus; and t_bound is Kv4.3 rather than Kv4.2, which was not')
    print('  quantified in the neuron-enriched preparation. A protein-level KChIP:Kv4 ratio would')
    print('  settle the first two, and that is the same co-IP that answers the question directly.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    pd.DataFrame(rows).to_csv(f'{OUT}/tables/surplus_predicted_half_life.csv', index=False)
    pd.DataFrame(inv).to_csv(f'{OUT}/tables/surplus_required_excess.csv', index=False)
    print(f'\nwrote {OUT}/tables/{{surplus_predicted_half_life,surplus_required_excess}}.csv')


if __name__ == '__main__':
    main()
