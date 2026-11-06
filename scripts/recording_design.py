#!/usr/bin/env python3
"""How many recordings, and what fraction of cells genuinely lack an A-current?

Two questions that decide whether the V404M bimodality prediction is testable at the rig.

PART 1 - TRUE ABSENCE. The mRNA says 6.5% of L6a Sla cells have no detectable Kv4. That is an UPPER BOUND on
the fraction with no A-current, for three reasons that all push the same way: a cell can express below the
detection threshold; transcription bursts, so a cell silent at capture may hold protein made earlier (mature
Kcnd2 mRNA is a long-lived pool, per the Oct 30 read-level work); and the threshold itself is arbitrary.
Estimated two ways:
  (a) threshold extrapolation - the any-Kv4 fraction swept from 1 read down, since a population of genuine
      low-expressers keeps shrinking as the threshold falls while true nulls do not;
  (b) depth-stratified zero counts - the fraction of cells with EXACTLY ZERO Kv4 reads among the deepest
      quartile of libraries. A cell sequenced to several million reads with not one Kv4 read is as close to
      a demonstrated null as this data can give, and (b) is the more defensible of the two.

PART 2 - POWER. In a V404M heterozygote, the Kv4.2-negative cells (~19% of those with any Kv4) should be
functionally indistinguishable from wild-type, so the mutant A-current distribution should be a MIXTURE
while the wild-type one is not. Blind bimodality tests are famously underpowered at unbalanced mixing
proportions, so the design simulated here uses the wild-type recordings as a fixed reference component and
asks whether the mutant distribution needs a wild-type-like component at all:

    H0: every mutant cell is shifted          H1: a fraction p of mutant cells sits on the wild-type curve

tested by likelihood ratio. That is far more powerful than asking a single sample whether it looks bimodal,
and it matches the experiment actually available - wild-type littermates are recorded anyway.

Effect size is left as a parameter (d = shift / within-type SD) because the V404M shift in V1/2 inactivation
is not established here. Read the table at whatever d the lab's own recordings support.

Usage:
    HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 recording_design.py [out_dir]
    TYPES='L6a Sla,L6a Mgp,L2/3 Ptgs2' python3 recording_design.py     # restrict part 1
    SKIP_DATA=1 python3 recording_design.py                            # power analysis only
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats, optimize

RNG       = np.random.default_rng(0)
OUT       = sys.argv[1] if len(sys.argv) > 1 else '.'
MIN_CELLS = int(os.environ.get('MIN_CELLS', 20))
N_SIM     = int(os.environ.get('N_SIM', 400))
P_TRUE    = float(os.environ.get('P_TRUE', 0.194))   # expected Kv4.2-negative fraction of ALL cells patched
KV4_GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3']


# ---------------------------------------------------------------- part 1
def true_absence(d, types=None):
    g = d.groupby(['sample_name', 'gene'])['reads'].sum().unstack(fill_value=0.0)
    kv4 = g.reindex(columns=KV4_GENES).fillna(0).sum(axis=1)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    depth = None
    for c in ('num_processed', 'total_reads', 'n_reads'):
        if c in d.columns:
            depth = d.groupby('sample_name')[c].first(); break
    if depth is None:
        depth = d.groupby('sample_name')['reads'].sum()
    X = pd.DataFrame({'kv4': kv4, 'depth': depth.reindex(kv4.index)}).join(meta)
    lab = 'cluster' if X['cluster'].nunique() > X['subclass'].nunique() else 'subclass'
    rows = []
    for ct, sub in X.groupby(lab):
        if len(sub) < MIN_CELLS or (types and ct not in types):
            continue
        r = dict(cell_type=ct, n_cells=len(sub))
        for thr in (10, 5, 3, 2, 1):
            r[f'pct_no_Kv4_at_{thr}'] = float((sub.kv4 < thr).mean() * 100)
        r['pct_exactly_zero'] = float((sub.kv4 == 0).mean() * 100)
        deep = sub[sub.depth >= sub.depth.quantile(0.75)]
        r['n_deep'] = len(deep)
        r['pct_zero_in_deepest_quartile'] = float((deep.kv4 == 0).mean() * 100) if len(deep) else np.nan
        r['median_depth_deep'] = float(deep.depth.median()) if len(deep) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- part 2
def fit_mixture_lrt(x, mu0, sd0, iters=200, tol=1e-8):
    """LRT for 'does this sample need a component sitting on the wild-type curve?'

    H0: a single Gaussian (every cell shifted).  H1: p * N(mu0, sd0) + (1-p) * N(mu1, sd1), with the
    wild-type component FIXED from the wild-type recordings. Returns (p_hat, LR statistic).

    Fitted by EM rather than a general optimiser: with one component fixed this is a three-parameter
    problem with closed-form M-steps, which is both faster and more reliable than multi-start Nelder-Mead
    (the power curve needs ~10^5 fits, which a general optimiser makes impractical).
    """
    x = np.asarray(x, float)
    n = len(x)
    ll0 = np.sum(stats.norm.logpdf(x, x.mean(), max(x.std(), 1e-9)))
    best = (0.0, ll0)
    f0 = stats.norm.pdf(x, mu0, sd0)
    for p_init in (0.1, 0.3):
        p = p_init
        mu1, sd1 = x.mean() + 0.5 * x.std(), max(x.std(), 1e-9)
        prev = -np.inf
        for _ in range(iters):
            f1 = stats.norm.pdf(x, mu1, sd1)
            mix = p * f0 + (1 - p) * f1
            ll = np.sum(np.log(np.maximum(mix, 1e-300)))
            r = p * f0 / np.maximum(mix, 1e-300)          # responsibility of the wild-type component
            p = float(np.clip(r.mean(), 1e-6, 1 - 1e-6))
            w = 1 - r
            sw = max(w.sum(), 1e-9)
            mu1 = float((w * x).sum() / sw)
            sd1 = float(np.sqrt(max((w * (x - mu1) ** 2).sum() / sw, 1e-12)))
            if ll - prev < tol:
                break
            prev = ll
        if ll > best[1]:
            best = (p, ll)
    return best[0], max(2 * (best[1] - ll0), 0.0)


def power_curve(ns, ds, p=P_TRUE, n_sim=N_SIM, alpha=0.05):
    """Power to reject 'no wild-type-like component' in the mutant sample.

    Wild-type cells ~ N(0,1). Mutant: a fraction p unshifted (no Kv4.2, so genotype-insensitive), the rest
    shifted by d. The wild-type component is treated as known, which is what recording littermates buys.
    Critical value is calibrated by simulation under H0 rather than assumed chi2, because the mixing
    proportion is on a boundary at H0 and the asymptotic null does not hold.
    """
    rows = []
    for d in ds:
        for n in ns:
            null = []
            for _ in range(n_sim):
                x = RNG.normal(d, 1.0, n)                       # H0: every cell shifted
                null.append(fit_mixture_lrt(x, 0.0, 1.0)[1])
            crit = np.quantile(null, 1 - alpha)
            hits, phats = 0, []
            for _ in range(n_sim):
                k = RNG.binomial(n, p)
                x = np.concatenate([RNG.normal(0.0, 1.0, k), RNG.normal(d, 1.0, n - k)])
                ph, st = fit_mixture_lrt(x, 0.0, 1.0)
                hits += st > crit
                phats.append(ph)
            rows.append(dict(d=d, n=n, power=hits / n_sim, crit=crit,
                             median_p_hat=float(np.median(phats))))
    return pd.DataFrame(rows)


def main():
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    types = [t.strip() for t in os.environ['TYPES'].split(',')] if os.environ.get('TYPES') else None

    if not os.environ.get('SKIP_DATA'):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from hip_common import load_ctx
        try:
            A = true_absence(load_ctx(), types)
        except Exception as e:
            print(f'[warn] could not load Tasic 2016: {e}')
            A = pd.DataFrame()
        if not A.empty:
            A.to_csv(f'{OUT}/tables/true_absence.csv', index=False)
            print('=' * 112)
            print('PART 1: WHAT FRACTION GENUINELY HAS NO Kv4 A-CURRENT?  (Tasic 2016)')
            print('  The >=5 read figure is an upper bound. "zero in deepest quartile" is the defensible')
            print('  estimate: a cell sequenced to millions of reads with not one Kv4 read.')
            print('=' * 112)
            cols = ['cell_type', 'n_cells', 'pct_no_Kv4_at_5', 'pct_no_Kv4_at_1', 'pct_exactly_zero',
                    'n_deep', 'pct_zero_in_deepest_quartile']
            print(A[cols].sort_values('pct_zero_in_deepest_quartile').to_string(index=False))

    print('\n' + '=' * 112)
    print(f'PART 2: RECORDINGS NEEDED TO DETECT THE MIXTURE  (expected unshifted fraction p = {P_TRUE:.3f})')
    print('  d = shift in V1/2 (or any parameter) divided by the within-type SD of that parameter.')
    print('  Power to reject "every mutant cell is shifted", using wild-type littermates as the reference')
    print('  component. n is the number of MUTANT cells; record a comparable number of wild-type.')
    print('=' * 112)
    ns = [20, 30, 40, 60, 80, 120, 200]
    ds = [1.0, 1.5, 2.0, 3.0]
    P = power_curve(ns, ds)
    P.to_csv(f'{OUT}/tables/recording_power.csv', index=False)
    piv = P.pivot_table(index='n', columns='d', values='power')
    print('\npower by n (rows) and effect size d (columns):')
    print(piv.round(2).to_string())
    print('\nn for 80% power:')
    for d in ds:
        s = P[(P.d == d) & (P.power >= 0.8)]
        print(f'  d = {d:.1f}  ->  ' + (f'n = {int(s.n.min())}' if len(s) else 'not reached by n = 200'))
    print('\nThe estimate of p itself is the more useful output: median p_hat by n and d')
    print(P.pivot_table(index='n', columns='d', values='median_p_hat').round(3).to_string())
    print(f'\nwrote {OUT}/tables/{{true_absence,recording_power}}.csv')


if __name__ == '__main__':
    main()
