#!/usr/bin/env python3
"""Which neurons sit in the Ca2+-sensitive regime of the KChIP:Kv4 ratio?

THE LAB'S CHO RESULT, which supplies the threshold. Varying the KChIP3a:Kv4.2 cDNA ratio and
chelating Ca2+ with BAPTA-AM gives a current-vs-ratio curve whose Ca2+ sensitivity is itself
ratio-dependent: at LOW KChIP3 the induction ratio (normal / BAPTA) is ~6x, and at HIGH KChIP3 it
collapses to ~1.4x. Ca2+ therefore matters for channel expression only when KChIP3:Kv4.2 is below
roughly 10%. The interpretation is that KChIP3 binds Kv4.2 much more weakly without Ca2+: with KChIP
in excess, mass action still saturates the channel at low affinity; with KChIP limiting, weak binding
means poor occupancy and no functional channel.

Supporting this from the same system: a Zn2+-site mutant Kv4.2 gives no current, co-expressed KChIP3a
rescues it, and mutating KChIP3's EF-hands abolishes that rescue while leaving its ability to promote
WT Kv4.2 expression intact. So there are two separable activities - general chaperoning, which is
EF-hand-independent, and stabilisation of a compromised T1 interface, which needs Ca2+-loaded KChIP3.

THE QUESTION THIS SCRIPT ASKS. In neurons the KChIP:Kv4 ratio is ~3 (i.e. ~300%), far above the 10%
threshold, so most cells should be Ca2+-insensitive for channel expression. But that is the MEDIAN.
If a tail of cells sits near or below 10%, those cells' A-current would be activity-dependent in a way
the bulk is not - and they would be the cells whose A-current is most labile.

NOISE. A per-cell ratio of two small counts is badly behaved, and a cell with few reads of either
family can land anywhere. Two guards: a minimum-count filter, and a Dirichlet(counts + 0.5) posterior
so that a low ratio is only called when the counts support it, rather than taken at face value. The
unguarded fraction is reported alongside so the inflation is visible.

Usage:  CTX_DATA=... python3 kchip_ratio_ca_sensitivity.py [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
RNG = np.random.default_rng(0)
THRESH = [0.10, 0.20, 0.50]          # KChIP:Kv4 ratios; 0.10 is the lab's CHO threshold
MIN_KV4 = float(os.environ.get('MIN_KV4', 10))   # CPM floor for a usable ratio
N_DRAW = 2000


def posterior_below(kchip, kv4, thr, draws=N_DRAW):
    """P(ratio < thr) under independent Gamma(count+0.5) posteriors on the two family totals.

    Using the CPM values as pseudo-counts is approximate, but the point is only to stop a cell with
    almost no reads from being scored as a confident low-ratio cell."""
    a = np.maximum(np.asarray(kchip, float), 0) + 0.5
    b = np.maximum(np.asarray(kv4, float), 0) + 0.5
    ga = RNG.gamma(np.broadcast_to(a, (draws,) + a.shape))
    gb = RNG.gamma(np.broadcast_to(b, (draws,) + b.shape))
    return ((ga / np.maximum(gb, 1e-9)) < thr).mean(axis=0)


def main():
    path = os.environ.get('CTX_DATA', 'tables/ctx_acurrent_index_per_cell.csv.gz')
    if not os.path.exists(path):
        sys.exit(f'missing {path}')
    d = pd.read_csv(path)
    d = d[d.cls != 'NonN'].copy()
    d = d[d.Kv4 >= MIN_KV4].copy()
    d['ratio'] = d.KChIP / d.Kv4.replace(0, np.nan)
    print(f'{len(d):,} neurons with Kv4 >= {MIN_KV4:g} CPM')

    print('\n' + '=' * 96)
    print('DISTRIBUTION of the per-cell KChIP:Kv4 ratio')
    print(f'  the CHO data puts the Ca2+-sensitive regime below ~0.10')
    print('=' * 96)
    for cls, s in d.groupby('cls'):
        q = s.ratio.quantile([.01, .05, .10, .25, .50, .75, .90]).round(2)
        print(f'  {cls:5} n={len(s):4d}   p1 {q[.01]:5.2f}  p5 {q[.05]:5.2f}  p10 {q[.10]:5.2f}  '
              f'Q1 {q[.25]:5.2f}  med {q[.50]:5.2f}  Q3 {q[.75]:5.2f}  p90 {q[.90]:5.2f}')

    print('\n' + '=' * 96)
    print('FRACTION OF CELLS IN THE Ca2+-SENSITIVE REGIME')
    print('  naive = ratio below threshold at face value;  guarded = posterior P(ratio < thr) > 0.5')
    print('=' * 96)
    rows = []
    print(f'  {"class":6} {"thr":>6} {"naive %":>9} {"guarded %":>11} {"n guarded":>10}')
    for cls, s in d.groupby('cls'):
        for thr in THRESH:
            naive = float((s.ratio < thr).mean() * 100)
            p = posterior_below(s.KChIP.values, s.Kv4.values, thr)
            guarded = float((p > 0.5).mean() * 100)
            print(f'  {cls:6} {thr:6.2f} {naive:8.1f}% {guarded:10.1f}% {int((p>0.5).sum()):10d}')
            rows.append(dict(cls=cls, threshold=thr, naive_pct=round(naive, 2),
                             guarded_pct=round(guarded, 2), n_guarded=int((p > 0.5).sum()),
                             n_cells=len(s)))

    # ---- which subclasses carry the low-ratio cells?
    print('\n' + '=' * 96)
    print('WHICH TYPES CARRY THE LOW-RATIO CELLS?  (threshold 0.10, guarded, types with n >= 15)')
    print('=' * 96)
    d['p_low'] = posterior_below(d.KChIP.values, d.Kv4.values, 0.10)
    d['low'] = d.p_low > 0.5
    g = (d.groupby(['cls', 'subclass'])
           .agg(n=('low', 'size'), n_low=('low', 'sum'), med_ratio=('ratio', 'median'))
           .reset_index())
    g = g[g.n >= 15].copy()
    g['pct_low'] = g.n_low / g.n * 100
    g = g.sort_values('pct_low', ascending=False)
    print(g.head(14).round(2).to_string(index=False))
    if g.pct_low.max() == 0:
        print('  No type has any cell confidently below 0.10.')

    print('\n' + '=' * 96)
    print('WHAT THIS MEANS FOR THE RECORDINGS')
    print('=' * 96)
    tot_low = float(d.low.mean() * 100)
    print(f'  Across all neurons, {tot_low:.2f}% sit confidently below the 0.10 ratio where the CHO')
    print('  data says Ca2+ chelation costs ~6-fold of the current.')
    if tot_low < 2:
        print('  That is a small minority, so the lab\'s own reading looks right: with KChIP in ~3-fold')
        print('  excess, Ca2+-dependent KChIP occupancy should NOT make most neurons\' A-current')
        print('  activity-labile. The Ca2+ sensitivity is real but mostly buffered by mass action.')
        print('  The cells that would be labile are a tail, not a population - which also means the')
        print('  effect would show up as a few outlying recordings rather than a shifted mean.')
    else:
        print('  That is a substantial fraction and would make Ca2+-dependent occupancy a live factor')
        print('  in ordinary recordings.')
    print('\n  IMPORTANT ASYMMETRY. This uses TOTAL KChIP against TOTAL Kv4. The CHO experiment was')
    print('  KChIP3a with Kv4.2 specifically. In cortical excitatory neurons KChIP3 is 76.5% of KChIP')
    print('  so the two nearly coincide; in GABAergic neurons KChIP3 is only 8.4%, so a GABAergic cell')
    print('  with an adequate TOTAL ratio could still be Ca2+-sensitive if KChIP3 is the subunit that')
    print('  carries the Ca2+ dependence and KChIP1/KChIP4 do not. Whether the EF-hand requirement is')
    print('  shared across KChIPs is not established here and would change this answer for inhibitory')
    print('  neurons specifically.')

    print('\n  AND A TESTABLE CONSEQUENCE FOR THE TURNOVER QUESTION. If Ca2+ strengthens KChIP binding,')
    print('  then lowering Ca2+ should increase the free KChIP pool, and - since free cytosolic KChIP')
    print('  is degraded - should SHORTEN KChIP protein half-life while leaving Kv4 half-life alone.')
    print('  BAPTA-AM or TTX plus a cycloheximide chase, blotting KChIP3 and Kv4.2, tests exchange')
    print('  directly, in the same preparation that produced the figure.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    pd.DataFrame(rows).to_csv(f'{OUT}/tables/kchip_ratio_ca_regime.csv', index=False)
    g.to_csv(f'{OUT}/tables/kchip_ratio_low_by_subclass.csv', index=False)
    d[['sample_name', 'cls', 'subclass', 'Kv4', 'KChIP', 'ratio', 'p_low']].to_csv(
        f'{OUT}/tables/kchip_ratio_per_cell.csv', index=False)
    print(f'\nwrote {OUT}/tables/kchip_ratio_{{ca_regime,low_by_subclass,per_cell}}.csv')


if __name__ == '__main__':
    main()
