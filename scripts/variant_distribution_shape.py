#!/usr/bin/env python3
"""Within one cell type, is a functionally decisive variant graded across cells or all-or-none?

This is the question that decides what to expect from electrophysiology, and it is NOT answered by the
within-type variance ratio (channel_noise_floor.py). A population mean of 20% DPP6K is consistent with

  (a) every cell at ~20%                 -> narrow V1/2-inactivation distribution, no ephys surprise
  (b) 20% of cells near 100%, rest at 0  -> BIMODAL V1/2, spanning most of the DPP6 shift (~20-25 mV)

Both give the same mean and both can give the same variance ratio, but they predict completely different
recordings: (b) says a single pyramidal population contains discrete biophysical subpopulations, findable by
patch-clamp alone without any sequencing; (a) says composition is a population-level average property and
per-cell measurement would add little.

For each (cell type x panel x variant) with enough cells this reports:
  - the per-cell fraction distribution (mean, SD, IQR, and the fraction of cells below 5% / above 50%)
  - 1- vs 2-component Gaussian mixture by BIC, as a bimodality test, plus the bimodality coefficient
  - EXCESS SD: the observed SD corrected for per-cell measurement uncertainty, via the same
    Dirichlet(counts + alpha) posterior used by channel_noise_floor.py. A variant whose apparent spread is
    all measurement error has excess_sd ~ 0 and predicts nothing, however bimodal the raw histogram looks --
    sparse count data produces spurious zero-vs-nonzero bimodality very easily, which is the main trap here.

Then, for the subunits whose biophysical effect is documented, it converts the excess spread into a predicted
spread in the recorded parameter under an explicitly stated and deliberately crude mapping (see EFFECTS). The
mapping assumes the parameter shifts linearly with the subunit's share of its panel, from 0 at 0% to the full
literature shift at 100%. That is certainly wrong in detail -- Kv4 channels are tetramers with up to four
auxiliary subunits, so occupancy is quantised and the dose-response probably saturates -- so the number is an
order-of-magnitude expectation, not a prediction. It is reported because "expect a few mV" and "expect 15 mV"
lead to different experiments.

Usage:
    HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 variant_distribution_shape.py [out_dir]
    TYPES='L6a Sla,L6a Mgp,L6a Car12' python3 variant_distribution_shape.py     # restrict to some types
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats
from sklearn.mixture import GaussianMixture

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

RNG       = np.random.default_rng(0)
ALPHA     = float(os.environ.get('ALPHA', 0.5))
N_DRAWS   = int(os.environ.get('N_DRAWS', 48))
MIN_READS = int(os.environ.get('MIN_READS', 200))   # per panel, per cell
MIN_CELLS = int(os.environ.get('MIN_CELLS', 15))
OUT       = sys.argv[1] if len(sys.argv) > 1 else '.'

PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}

# Documented effect of going from none to full occupancy of this subunit, and which recorded parameter it
# moves. Signs and magnitudes are approximate and from the literature on heterologous co-expression, which
# is where the crudeness of this mapping mostly comes from. Edit freely -- these are the assumptions, and
# they are meant to be argued with rather than trusted.
EFFECTS = {
    'DPP6S':   ('V1/2 inactivation (mV)', -25.0),
    'DPP6K':   ('V1/2 inactivation (mV)', -25.0),
    'DPP6L':   ('V1/2 inactivation (mV)', -25.0),
    'DPP6a':   ('V1/2 inactivation (mV)', -25.0),
    'DPP10a':  ('V1/2 inactivation (mV)', -20.0),
    'KChIP4a': ('inactivation completeness (frac)', -0.6),   # KIS domain suppresses inactivation
}


def panel_counts(d, panel):
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    C = V.reindex(columns=PANEL_V[panel]).fillna(0.0).round().astype(int)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    return C, meta.reindex(C.index)


def posterior_sd(counts_row, idx, alpha=ALPHA, draws=N_DRAWS):
    """SD of one variant's fraction under Dirichlet(counts+alpha) -- the per-cell measurement uncertainty."""
    a = np.asarray(counts_row, float) + alpha
    G = RNG.gamma(shape=np.broadcast_to(a, (draws, len(a))), scale=1.0)
    P = G / G.sum(axis=1, keepdims=True)
    return P[:, idx].std()


def shape_stats(frac, noise_sd):
    """Distribution shape of one variant's per-cell fractions, with measurement error removed."""
    n = len(frac)
    obs_sd = frac.std(ddof=1)
    # variance decomposition: observed = biological + measurement (independent)
    excess_var = max(obs_sd ** 2 - np.mean(noise_sd ** 2), 0.0)
    excess_sd = np.sqrt(excess_var)
    out = dict(n_cells=n, mean=frac.mean(), median=np.median(frac), obs_sd=obs_sd,
               median_noise_sd=float(np.median(noise_sd)), excess_sd=excess_sd,
               iqr=float(np.subtract(*np.percentile(frac, [75, 25]))),
               frac_below_5pct=float((frac < 0.05).mean()), frac_above_50pct=float((frac > 0.50).mean()))
    # bimodality coefficient: >5/9 is the conventional flag
    if obs_sd > 0:
        g = stats.skew(frac); k = stats.kurtosis(frac, fisher=True)
        out['bimodality_coef'] = (g ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))) if n > 3 else np.nan
    else:
        out['bimodality_coef'] = np.nan
    # 1 vs 2 component Gaussian mixture by BIC
    X = frac.reshape(-1, 1)
    try:
        b1 = GaussianMixture(1, random_state=0).fit(X).bic(X)
        b2 = GaussianMixture(2, random_state=0, n_init=3).fit(X).bic(X)
        out['bic_1comp'], out['bic_2comp'], out['delta_bic'] = b1, b2, b1 - b2
        # NB: named shape_call, not shape -- 'shape' collides with the pandas Series attribute, so
        # row.shape returns the shape tuple and never the column.
        out['shape_call'] = 'bimodal' if (b1 - b2) > 10 else ('graded' if out['bimodality_coef'] and
                                                              out['bimodality_coef'] < 0.555 else 'unclear')
    except Exception:
        out['bic_1comp'] = out['bic_2comp'] = out['delta_bic'] = np.nan
        out['shape_call'] = 'unclear'
    return out


def analyse(name, d, types=None):
    rows = []
    for panel, variants in PANEL_V.items():
        C, meta = panel_counts(d, panel)
        tot = C.sum(axis=1).to_numpy()
        keep = tot >= MIN_READS
        if keep.sum() < MIN_CELLS:
            continue
        Ck, mk = C[keep], meta[keep]
        Pk = Ck.to_numpy() / Ck.to_numpy().sum(axis=1, keepdims=True)
        label = 'cluster' if mk['cluster'].nunique() > mk['subclass'].nunique() else 'subclass'
        for ct, idx in mk.groupby(label).groups.items():
            if types and ct not in types:
                continue
            ii = mk.index.get_indexer(idx)
            if len(ii) < MIN_CELLS:
                continue
            for vi, v in enumerate(variants):
                if v not in EFFECTS and panel != 'DPP':
                    continue                      # keep the table focused on the decisive subunits
                frac = Pk[ii, vi]
                if frac.max() < 0.02:
                    continue                      # variant essentially absent from this type
                nsd = np.array([posterior_sd(Ck.to_numpy()[i], vi) for i in ii])
                s = shape_stats(frac, nsd)
                s.update(dataset=name, cell_type=ct, panel=panel, variant=v)
                if v in EFFECTS:
                    param, full = EFFECTS[v]
                    s['predicted_param'] = param
                    s['predicted_spread'] = abs(full) * s['excess_sd']        # 1 SD
                    s['predicted_range_p5_p95'] = abs(full) * float(
                        np.subtract(*np.percentile(frac, [95, 5])))
                rows.append(s)
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
    types = [t.strip() for t in os.environ['TYPES'].split(',')] if os.environ.get('TYPES') else None
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

    out = []
    for name, d in ds.items():
        r = analyse(name, d, types)
        print(f'{name}: {len(r)} (type x variant) combinations with >= {MIN_CELLS} cells at '
              f'>= {MIN_READS} reads/panel')
        out.append(r)
    R = pd.concat(out, ignore_index=True)
    if R.empty:
        sys.exit('nothing passed the thresholds -- lower MIN_READS or MIN_CELLS')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/variant_distribution_shape.csv', index=False)

    print('\n' + '=' * 108)
    print('ALL-OR-NONE vs GRADED, for the subunits with documented biophysical effects')
    print('=' * 108)
    E = R[R['variant'].isin(EFFECTS)].copy()
    if not E.empty:
        E = E.sort_values('predicted_spread', ascending=False)
        print(f'{"dataset":7} {"cell type":16} {"variant":8} {"n":>4} {"mean":>6} {"excSD":>6} '
              f'{"noise":>6} {"<5%":>5} {">50%":>5} {"dBIC":>7} {"shape":9} predicted 1SD spread')
        for _, r in E.head(40).iterrows():
            print(f'{r.dataset:7} {r.cell_type[:16]:16} {r.variant:8} {int(r.n_cells):4d} '
                  f'{r["mean"]:6.3f} {r.excess_sd:6.3f} {r.median_noise_sd:6.3f} '
                  f'{r.frac_below_5pct:5.2f} {r.frac_above_50pct:5.2f} {r.delta_bic:7.1f} '
                  f'{r["shape_call"]:9} {r.predicted_spread:5.1f} {r.predicted_param}')
    print('\n' + '=' * 108)
    print('HOW TO READ THIS')
    print('=' * 108)
    print('excess_sd is the per-cell spread AFTER removing measurement uncertainty. It, not the raw SD,')
    print('is what predicts anything. A variant with excess_sd ~ 0 predicts no ephys variance however')
    print('bimodal its histogram looks -- sparse counts manufacture zero-vs-nonzero bimodality easily.')
    print('frac_below_5pct and frac_above_50pct together describe the all-or-none pattern directly: a')
    print('variant with both high is present in some cells and absent in others, which is the case that')
    print('predicts discrete biophysical subpopulations findable by patch-clamp WITHOUT sequencing.')
    print('delta_bic > 10 favours two components over one. Predicted spreads assume a linear subunit-share')
    print('to parameter mapping (see EFFECTS in the source) and are order-of-magnitude only: Kv4 channels')
    print('are tetramers, so occupancy is quantised and the real dose-response likely saturates.')
    print(f'\nwrote {OUT}/tables/variant_distribution_shape.csv')


if __name__ == '__main__':
    main()
