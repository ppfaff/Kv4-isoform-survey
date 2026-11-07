#!/usr/bin/env python3
"""The KChIP:Kv4 ratio done properly - on raw read counts, with a posterior that actually constrains.

WHY THIS REPLACES THE EARLIER VERSION. `kchip_ratio_ca_sensitivity.py` computed the same ratio from
CPM values and put a Dirichlet guard on them. That guard was inert: CPM is already normalised, so
using it as a pseudo-count asserts far more evidence than a shallow library contains, and the guarded
and naive percentages came out identical. salmon's NumReads is the quantity a posterior can be built
on, and `pipeline/extract_panel_counts.py` now supplies it per cell.

THE QUESTION. The lab's CHO data (BAPTA-AM x KChIP3a:Kv4.2 ratio) puts the Ca2+-sensitive regime
below a ratio of about 0.10: there, chelating Ca2+ costs ~6-fold of the current, because KChIP binds
Kv4 much more weakly without Ca2+ and only an excess of KChIP rescues occupancy by mass action. So:
what fraction of real neurons sit below that threshold, and does the answer survive honest error bars?

TWO RATIOS ARE COMPUTED, and they are not interchangeable.
  total KChIP : total Kv4   - right if the Ca2+ requirement is shared across the EF-hand family.
  Kcnip3      : Kcnd2       - the exact analogue of the CHO experiment, right if it is KChIP3-specific.
Both are reported because the choice changes the answer for GABAergic neurons by a large factor, and
nothing here settles which is correct.

THE POSTERIOR. Gamma(count + 0.5) on each family total, ratio of draws, then P(ratio < threshold). A
cell is called low only when the counts support it. Cells are also depth-stratified, because a cell
with almost no panel reads is uninformative about its own ratio however the point estimate falls.

Usage:  python3 kchip_ratio_rawcounts.py <panel_counts_per_cell.csv> [runlist.tsv] [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else 'panel_counts_per_cell.csv'
# second argument may be the runlist or the out_dir; decide by what it actually is
_a2 = sys.argv[2] if len(sys.argv) > 2 else None
RUNLIST = _a2 if (_a2 and os.path.isfile(_a2)) else None
OUT = sys.argv[3] if len(sys.argv) > 3 else (_a2 if (_a2 and os.path.isdir(_a2)) else '.')
RNG = np.random.default_rng(0)
THRESH = [0.10, 0.20, 0.50]
N_DRAW = 4000
MIN_PANEL = int(os.environ.get('MIN_PANEL', 20))     # total panel reads for a cell to be usable

KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']
KCHIP = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']


def p_below(num, den, thr, draws=N_DRAW):
    """P(num/den < thr) under independent Gamma(count + 0.5) posteriors."""
    a = np.asarray(num, float) + 0.5
    b = np.asarray(den, float) + 0.5
    ga = RNG.gamma(np.broadcast_to(a, (draws,) + a.shape))
    gb = RNG.gamma(np.broadcast_to(b, (draws,) + b.shape))
    return ((ga / gb) < thr).mean(axis=0)


def main():
    d = pd.read_csv(CSV)
    g = d.pivot_table(index='run', columns='gene', values='NumReads', aggfunc='sum').fillna(0.0)
    for c in KV4 + KCHIP:
        if c not in g.columns:
            g[c] = 0.0
    meta = d.groupby('run').agg(depth=('num_processed', 'first'), mapped=('percent_mapped', 'first'))
    X = pd.DataFrame({'kv4': g[KV4].sum(axis=1), 'kchip': g[KCHIP].sum(axis=1),
                      'kcnip3': g['Kcnip3'], 'kcnd2': g['Kcnd2']}).join(meta)
    X['panel_reads'] = g.sum(axis=1)
    print(f'{len(X):,} cells; median sequencing depth {X.depth.median():,.0f}; '
          f'median panel reads {X.panel_reads.median():,.0f}')

    if RUNLIST and os.path.exists(RUNLIST):
        rl = pd.read_csv(RUNLIST, sep='\t')
        key = next((c for c in ('sample_name', 'run', 'run_accession') if c in rl.columns), None)
        keep = [c for c in ('subclass_label', 'cluster_label', 'region_label') if c in rl.columns]
        if key and keep:
            X = X.join(rl.set_index(key)[keep])
            print(f'annotated {X.subclass_label.notna().sum():,} of {len(X):,} cells from {RUNLIST}')
    else:
        print('[no runlist supplied - reporting pooled results only; '
              'pass hip_filtered_runlist.tsv for per-type breakdown]')

    U = X[(X.kv4 + X.kchip) >= MIN_PANEL].copy()
    print(f'{len(U):,} cells with >= {MIN_PANEL} Kv4+KChIP reads (the rest cannot constrain a ratio)')

    for name, num, den in (('total KChIP : total Kv4', 'kchip', 'kv4'),
                           ('Kcnip3 : Kcnd2  (the CHO analogue)', 'kcnip3', 'kcnd2')):
        sub = U[U[den] + U[num] >= MIN_PANEL]
        r = sub[num] / sub[den].replace(0, np.nan)
        print('\n' + '=' * 92)
        print(f'{name}   n = {len(sub):,}')
        print('=' * 92)
        q = r.quantile([.05, .10, .25, .50, .75, .90]).round(3)
        print(f'  p5 {q[.05]}  p10 {q[.10]}  Q1 {q[.25]}  median {q[.50]}  Q3 {q[.75]}  p90 {q[.90]}')
        print(f'\n  {"thr":>6} {"naive %":>9} {"guarded %":>11} {"n guarded":>10}   '
              f'(guarded = posterior P(ratio < thr) > 0.5)')
        for thr in THRESH:
            naive = float((r < thr).mean() * 100)
            p = p_below(sub[num].values, sub[den].values, thr)
            guarded = float((p > 0.5).mean() * 100)
            print(f'  {thr:6.2f} {naive:8.1f}% {guarded:10.1f}% {int((p>0.5).sum()):10d}')
        # depth stratification on the 0.10 threshold
        p10 = p_below(sub[num].values, sub[den].values, 0.10)
        sub = sub.assign(low=p10 > 0.5)
        qs = pd.qcut(sub.panel_reads, 4, labels=['Q1 lowest', 'Q2', 'Q3', 'Q4 highest'])
        print('\n  by panel-read quartile (does the low-ratio call survive depth?):')
        for lab, s in sub.groupby(qs, observed=True):
            print(f'    {str(lab):11} n={len(s):5d}  median panel reads {s.panel_reads.median():7.0f}  '
                  f'low {s.low.mean()*100:5.1f}%')

    # per-type, if annotated
    if 'subclass_label' in X.columns and X.subclass_label.notna().any():
        U2 = U[U.subclass_label.notna()].copy()
        U2['low'] = p_below(U2.kchip.values, U2.kv4.values, 0.10) > 0.5
        t = (U2.groupby('subclass_label')
               .agg(n=('low', 'size'), n_low=('low', 'sum'),
                    med_ratio=('kchip', 'size')).reset_index())
        t['med_ratio'] = [float((U2[U2.subclass_label == s].kchip /
                                 U2[U2.subclass_label == s].kv4.replace(0, np.nan)).median())
                          for s in t.subclass_label]
        t['pct_low'] = t.n_low / t.n * 100
        t = t[t.n >= 15].sort_values('pct_low', ascending=False)
        print('\n' + '=' * 92)
        print('BY SUBCLASS (total ratio, threshold 0.10, guarded, n >= 15)')
        print('=' * 92)
        print(t.round(2).to_string(index=False))
        t.to_csv(f'{OUT}/tables/kchip_ratio_raw_by_subclass.csv', index=False)

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    X.to_csv(f'{OUT}/tables/panel_counts_per_cell_summary.csv')
    print(f'\nwrote {OUT}/tables/panel_counts_per_cell_summary.csv')


if __name__ == '__main__':
    main()
