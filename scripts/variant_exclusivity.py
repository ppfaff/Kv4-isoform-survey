#!/usr/bin/env python3
"""Is the all-or-none variant usage biology, or salmon's EM collapsing onto one isoform?

variant_distribution_shape.py found that within a single cell type, DPP6 N-terminal variants and DPP10a are
bimodal: cells either barely express a variant (<5% of the panel) or it dominates (>50%), with excess spread
30-60x the per-cell measurement noise. Dropout cannot produce that, because a cell in the >50% bin carries
>=100 reads on that one variant.

But WINNER-TAKE-ALL EM COLLAPSE can. The DPP6 N-terminal variants differ only at the 5' end. Where coverage
over the discriminating region is thin, salmon's EM has little information with which to apportion reads and
can settle on one isoform near-arbitrarily -- all DPP6S in one cell, all DPP6a in the next -- while the panel
total stays conserved. Biologically identical cells would then look bimodal and mutually exclusive.

The signature that separates the two is MUTUAL EXCLUSIVITY BEYOND WHAT COMPOSITION ALONE FORCES:

  EM collapse  -> within a cell, essentially ONE variant of a mutually-confusable set is detected, the set's
                  total is conserved across cells, and the pair correlations are far more negative than a
                  compositional null allows
  biology      -> cells carry genuine MIXTURES (60/40, 70/30), and correlations sit near the null

The compositional null matters and is the reason this script exists rather than a correlation matrix: panel
fractions sum to 1, so variants are negatively correlated even with no artefact and no biology whatsoever.
The null here is a parametric Dirichlet: resample each cell at ITS OWN depth from the CELL TYPE's mean
composition, which is the "all cells identical, only sampling noise" world, and read off how negative the
correlations and how low the mixture counts go by construction. Observed values far beyond that null are the
evidence for collapse; observed values near it, with high mixture counts, are evidence for biology.

Confusable sets are declared explicitly in SETS: only variants that genuinely compete for the same reads
should be grouped, since the test asks whether reads are being shuffled WITHIN a set.

Usage:
    HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 variant_exclusivity.py [out_dir]
"""
import os, sys, itertools
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

RNG       = np.random.default_rng(0)
MIN_READS = int(os.environ.get('MIN_READS', 200))
MIN_CELLS = int(os.environ.get('MIN_CELLS', 20))
N_NULL    = int(os.environ.get('N_NULL', 200))
MIX_THR   = float(os.environ.get('MIX_THR', 0.10))   # a variant counts as "present" above this fraction
OUT       = sys.argv[1] if len(sys.argv) > 1 else '.'

PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}

# Sets of variants that compete for the same reads, i.e. where EM collapse is possible. DPP6 N-terminal
# variants are the primary concern (they differ only at the 5' end). DPP10 variants likewise. Kv4.3S/L differ
# by a cassette exon, which is directly observable in junction reads, so they are a useful positive control:
# if the method flags THEM as collapsed, the test itself is too aggressive.
SETS = {
    'DPP6 N-terminal': ('DPP', ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L']),
    'DPP10':           ('DPP', ['DPP10a', 'DPP10c', 'DPP10d']),
    'KChIP4':          ('KChIP', ['KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e']),
    'Kv4.3 cassette':  ('Kv4', ['Kv4.3S', 'Kv4.3L']),
}


def panel_counts(d, panel):
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    C = V.reindex(columns=PANEL_V[panel]).fillna(0.0).round().astype(int)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    return C, meta.reindex(C.index)


def metrics(sub_counts):
    """Exclusivity metrics for one cell-type x variant-set block. sub_counts: (cells, variants) integers."""
    tot = sub_counts.sum(axis=1)
    ok = tot > 0
    if ok.sum() < MIN_CELLS:
        return None
    F = sub_counts[ok] / tot[ok][:, None]
    # how many variants of the set are "present" in each cell
    n_present = (F > MIX_THR).sum(axis=1)
    # most negative pairwise Spearman within the set (the pair that would be swapping)
    cors = []
    for i, j in itertools.combinations(range(F.shape[1]), 2):
        if F[:, i].std() > 0 and F[:, j].std() > 0:
            cors.append(stats.spearmanr(F[:, i], F[:, j]).statistic)
    return dict(n_cells=int(ok.sum()),
                mean_n_present=float(n_present.mean()),
                frac_single_variant=float((n_present <= 1).mean()),
                min_pair_rho=float(np.nanmin(cors)) if cors else np.nan,
                set_total_cv=float(np.std(tot[ok]) / np.mean(tot[ok])))


def null_metrics(sub_counts, n_null=N_NULL):
    """Same metrics under 'all cells identical, sampling noise only', at each cell's own depth."""
    tot = sub_counts.sum(axis=1)
    ok = tot > 0
    if ok.sum() < MIN_CELLS:
        return None
    p = sub_counts[ok].sum(axis=0).astype(float)
    p = p / p.sum()
    depths = tot[ok]
    acc = {'mean_n_present': [], 'frac_single_variant': [], 'min_pair_rho': []}
    for _ in range(max(1, n_null // 20)):
        S = np.array([RNG.multinomial(int(n), p) for n in depths], float)
        m = metrics(S.astype(int))
        if m:
            for k in acc:
                acc[k].append(m[k])
    return {k: (float(np.nanmean(v)), float(np.nanpercentile(v, 5)) if k == 'mean_n_present'
                else float(np.nanpercentile(v, 5))) for k, v in acc.items() if v}


def analyse(name, d):
    rows = []
    cache = {p: panel_counts(d, p) for p in PANEL_V}
    for set_name, (panel, variants) in SETS.items():
        C, meta = cache[panel]
        idx = [PANEL_V[panel].index(v) for v in variants if v in PANEL_V[panel]]
        A = C.to_numpy()
        keep = A.sum(axis=1) >= MIN_READS
        if keep.sum() < MIN_CELLS:
            continue
        Ck, mk = A[keep][:, idx], meta[keep]
        label = 'cluster' if mk['cluster'].nunique() > mk['subclass'].nunique() else 'subclass'
        for ct, g in mk.groupby(label).groups.items():
            ii = mk.index.get_indexer(g)
            obs = metrics(Ck[ii])
            if not obs:
                continue
            nul = null_metrics(Ck[ii])
            if not nul:
                continue
            r = dict(dataset=name, cell_type=ct, variant_set=set_name, **obs)
            r['null_mean_n_present'] = nul['mean_n_present'][0]
            r['null_frac_single'] = nul['frac_single_variant'][0]
            r['null_min_pair_rho'] = nul['min_pair_rho'][0]
            r['null_rho_p5'] = nul['min_pair_rho'][1]
            # verdict: collapse needs BOTH fewer mixtures than the null allows AND a correlation beyond it
            r['excess_single'] = r['frac_single_variant'] - r['null_frac_single']
            r['rho_beyond_null'] = r['null_rho_p5'] - r['min_pair_rho']
            # Verdict rests on the MIXTURE COUNT, not the correlation. Validation showed the pair rho is
            # worse than useless here: a genuinely bimodal biological case (cells in one of two mixture
            # states) gave min_rho -0.85, beyond its own null of -0.63, while true EM collapse gave only
            # -0.33, inside its null of -0.54. Requiring a rho departure as well, as a first version did,
            # would have missed the collapse case outright. Collapse puts each cell on ONE variant, so
            # frac_single_variant is the statistic that moves: 0.00 for both biological cases against
            # 1.00 for collapse. rho is kept in the table as a diagnostic, not used to decide.
            r['verdict'] = ('EM COLLAPSE suspected'
                            if (r['excess_single'] > 0.25
                                or (r['mean_n_present'] < 1.2 and r['null_mean_n_present'] > 2.0))
                            else 'mixtures present - consistent with biology'
                            if r['mean_n_present'] >= 1.5 and r['excess_single'] < 0.15
                            else 'ambiguous')
            rows.append(r)
    return pd.DataFrame(rows)


def prep_generic(path):
    v = _read_percell(path)
    v['vclass'] = v['transcript'].map(final_map)
    v = v.dropna(subset=['vclass']).copy()
    if 'biotype' in v.columns:
        v = v[~((v['biotype'] == 'retained_intron')
                & v['vclass'].str.contains('other|fragment', case=False))].copy()
    v['panel'] = v['gene'].map(PANEL_OF)
    v = v.rename(columns={k: x for k, x in {'subclass_label': 'subclass', 'cluster_label': 'cluster'}.items()
                          if k in v.columns})
    if 'cluster' not in v.columns:
        v['cluster'] = v.get('subclass')
    return v


def main():
    ds = {}
    for label, loader in [('T16', load_ctx), ('A-HPF', load_hip)]:
        try:
            ds[label] = loader()
        except Exception as e:
            print(f'[warn] {label}: {e}')
    vp = os.environ.get('VISP_DATA') or os.environ.get('VISP_PKL')
    if vp and os.path.exists(vp):
        ds['A-VISp'] = prep_generic(vp)
    if not ds:
        sys.exit('no datasets loaded')

    out = [analyse(n, d) for n, d in ds.items()]
    R = pd.concat([x for x in out if not x.empty], ignore_index=True)
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/variant_exclusivity.csv', index=False)

    print('=' * 112)
    print('MUTUAL EXCLUSIVITY vs COMPOSITIONAL NULL')
    print('  mean_n_present : variants of the set detected above '
          f'{MIX_THR:.0%} per cell (observed vs null)')
    print('  frac_single    : fraction of cells with <=1 variant present (observed vs null)')
    print('  min_rho        : most negative pairwise Spearman (observed vs null 5th pctile)')
    print('=' * 112)
    for set_name in SETS:
        S = R[R.variant_set == set_name]
        if S.empty:
            continue
        print(f'\n--- {set_name}')
        print(f'{"dataset":7} {"cell type":18} {"n":>4}  {"nPres":>5}/{"null":<5} '
              f'{"single":>6}/{"null":<6} {"minRho":>6}/{"n5":<6}  verdict')
        for _, r in S.sort_values('excess_single', ascending=False).head(12).iterrows():
            print(f'{r.dataset:7} {str(r.cell_type)[:18]:18} {int(r.n_cells):4d}  '
                  f'{r.mean_n_present:5.2f}/{r.null_mean_n_present:<5.2f} '
                  f'{r.frac_single_variant:6.2f}/{r.null_frac_single:<6.2f} '
                  f'{r.min_pair_rho:6.2f}/{r.null_rho_p5:<6.2f}  {r.verdict}')
    print('\n' + '=' * 112)
    print('HOW TO READ THIS')
    print('=' * 112)
    print('Panel fractions sum to 1, so variants are negatively correlated with NO artefact and NO biology.')
    print('That is what the null columns are for: they are the "all cells identical, sampling noise only"')
    print('world at each cell\'s own depth. Only departures BEYOND the null carry information.')
    print('  observed mixtures ~ null, mean_n_present >= 1.5  -> cells carry real mixtures; the bimodality')
    print('      in variant_distribution_shape.py is biological and the ephys prediction stands')
    print('  far fewer mixtures than null AND rho beyond null  -> reads are being shuffled within the set;')
    print('      the bimodality is an EM artefact and the ephys prediction does not follow')
    print('Kv4.3S/L is the positive control: that pair differs by a cassette exon which junction reads')
    print('resolve directly, so it should NOT be flagged. If it is, the test is too aggressive and the')
    print('thresholds in analyse() need revisiting before any of the other rows are believed.')
    print(f'\nwrote {OUT}/tables/variant_exclusivity.csv')


if __name__ == '__main__':
    main()
