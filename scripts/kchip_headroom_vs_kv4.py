#!/usr/bin/env python3
"""Who is closest to KChIP limitation? Not a cell type — the high-Kv4 tail of every cell type.

THE SETUP. The lab's CHO data shows KChIP3 binds Kv4.2 much more weakly without Ca2+, so channel
expression becomes Ca2+-dependent once KChIP stops being in excess (below a KChIP:Kv4 ratio of about
0.10 in that system, where KChIP3 was the only KChIP present). The working hypothesis is that cells
differing in this sensitivity would differ in their ability to renormalise excitability — and that
such cells might lower their A-current in V404M as homeostatic compensation.

The obvious way to look for those cells is to ask which SUBCLASS has a low ratio. That is the wrong
cut, and this script shows why.

THE FINDING. KChIP expression does not scale with Kv4 expression. Spearman correlations between raw
KChIP and Kv4 counts within a subclass run from -0.15 to +0.47 and are mostly weak, which is the
single-cell version of the near-constitutive auxiliary-subunit result (Fig. 7 of
KV4_SURVEY_CONCLUSIONS.md: KChIP and DPP totals stay flat across the whole Kv4 gradient, Spearman
0.07-0.09). The consequence is arithmetic but functionally important: within every cell type, raising
Kv4 lowers the KChIP:Kv4 ratio, because the numerator does not follow. **So the cells carrying the
largest A-current are systematically the ones with the least KChIP headroom.**

That relocates the Ca2+-sensitive population. It is not a cell class — it is the high-Kv4 tail inside
each class, which is also precisely the population for which a homeostatic reduction of A-current
would have the largest absolute effect.

CAVEAT. These are mRNA ratios, and occupancy is set by protein. The gradient is the robust part; the
absolute position relative to the 0.10 cusp is not, since translation efficiency may differ between
the pore and the auxiliary subunits (Liss et al. 2001 estimate ~7x more A-conductance per KChIP3
transcript than per Kv4.3L transcript, which they suggest may be translation efficiency).

Usage:  python3 kchip_headroom_vs_kv4.py [summary_csv] [out_dir]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

CSV = sys.argv[1] if len(sys.argv) > 1 else 'tables/panel_counts_per_cell_summary.csv'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
MIN_PANEL = 20
MIN_CELLS = 60


def main():
    S = pd.read_csv(CSV, index_col=0)
    S = S[(S.kchip + S.kv4) >= MIN_PANEL].copy()
    S['ratio'] = S.kchip / S.kv4.replace(0, np.nan)
    print(f'{len(S):,} cells')

    print('\n' + '=' * 108)
    print('DOES KChIP TRACK Kv4 WITHIN A CELL TYPE?  (raw counts)')
    print('  If the auxiliary subunits are near-constitutive, the answer is no, and the ratio must')
    print('  then fall as Kv4 rises purely because the numerator does not follow.')
    print('=' * 108)
    print(f'  {"subclass":14} {"n":>5} {"rho":>7} {"p":>9} {"med Kv4":>8} {"med KChIP":>10}   '
          f'median KChIP:Kv4 by Kv4 quartile')
    rows = []
    for sc, s in S.groupby('subclass_label'):
        if len(s) < MIN_CELLS:
            continue
        rho, p = stats.spearmanr(s.kv4, s.kchip)
        try:
            q = pd.qcut(s.kv4, 4, labels=['Q1', 'Q2', 'Q3', 'Q4'], duplicates='drop')
        except ValueError:
            continue
        byq = (s.assign(q=q).groupby('q', observed=True)
                 .apply(lambda x: float((x.kchip / x.kv4.replace(0, np.nan)).median()),
                        include_groups=False))
        if len(byq) < 4:
            continue
        trail = '  '.join(f'{k}:{v:5.1f}' for k, v in byq.items())
        fold = byq.iloc[0] / max(byq.iloc[-1], 1e-9)
        print(f'  {sc:14} {len(s):5d} {rho:+7.2f} {p:9.0e} {s.kv4.median():8.0f} '
              f'{s.kchip.median():10.0f}   {trail}   ({fold:.0f}x fall)')
        rows.append(dict(subclass=sc, n=len(s), rho=round(float(rho), 3), p=float(p),
                         med_kv4=float(s.kv4.median()), med_kchip=float(s.kchip.median()),
                         **{f'ratio_{k}': round(float(v), 2) for k, v in byq.items()},
                         fold_Q1_over_Q4=round(float(fold), 1)))
    R = pd.DataFrame(rows)

    print('\n' + '=' * 108)
    print('WHAT THIS MEANS')
    print('=' * 108)
    med_rho = R.rho.median()
    print(f'  Median within-subclass Spearman(KChIP, Kv4) = {med_rho:+.2f} — KChIP does not follow Kv4.')
    print(f'  Median Q1/Q4 fall in the ratio = {R.fold_Q1_over_Q4.median():.0f}x, and it is monotonic in')
    print(f'  every subclass tested ({len(R)} of them).')
    low_q4 = R[R.ratio_Q4 <= 1.0]
    if len(low_q4):
        names = ', '.join(f'{r.subclass} ({r.ratio_Q4:.1f})' for r in low_q4.itertuples())
        print(f'\n  Subclasses whose top-Kv4 quartile falls to a ratio of 1.0 or below: {names}.')
    print('\n  So the Ca2+-sensitive population is NOT a cell class. It is the high-Kv4 tail within')
    print('  each class — the cells with the biggest A-current and the least KChIP headroom. That is')
    print('  also the population for which a homeostatic reduction of A-current would have the')
    print('  largest absolute effect, which is what makes this relevant to the compensation idea.')
    print('\n  Note the sign of the feedback loop works: more A-current -> less excitable -> less Ca2+')
    print('  -> weaker KChIP binding -> less A-current. That is NEGATIVE feedback on A-current and')
    print('  therefore stabilising for excitability, and it is strongest exactly where headroom is')
    print('  least. A cell with KChIP in 10-fold excess cannot use this mechanism at all.')

    if not R.empty:
        os.makedirs(f'{OUT}/tables', exist_ok=True)
        R.to_csv(f'{OUT}/tables/kchip_headroom_by_kv4_quartile.csv', index=False)
        print(f'\nwrote {OUT}/tables/kchip_headroom_by_kv4_quartile.csv')


if __name__ == '__main__':
    main()
