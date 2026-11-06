#!/usr/bin/env python3
"""Is per-cell variance in Kv4-complex composition biological, or just sampling noise?

Diagnostic to run BEFORE any grouping analysis, because it decides whether the question is answerable in a
given dataset at all. The first run of channel_phenotype_groups.py returned "split" for 24 of 28 groupings
including Allen subclasses against their own clusters, on R^2 values of 0.01-0.22 -- i.e. significance driven
by sample size. The unresolved question underneath was whether the large within-group scatter is real
heterogeneity or measurement error, and nothing in that script could tell the two apart.

METHOD. Split each cell's reads for a panel into two independent halves (binomial split per variant), take
the CLR of each half, and measure the distance between them. That is a technical-replicate distance for that
exact cell at that exact depth, assuming nothing about cell types. Scaling: if each half has error e_h about
the cell's true composition and the two errors are independent, E|a-b|^2 = 2*e_h^2, so |a-b| = sqrt(2)*e_h;
and error falls as 1/sqrt(n), so the full-depth error is e_h/sqrt(2) = |a-b|/2. The reported noise floor is
therefore half the split distance.

Compare that against the within-group spread (distance from each cell to its group's centroid):

    spread / noise  ~ 1   the scatter IS measurement error; cells of this type are homogeneous and no
                          grouping verdict can be drawn at this depth
    spread / noise >> 1   genuine biological heterogeneity within the type

The ratio is dimensionless, so unlike a thinning control it is comparable across datasets of different depth
-- which matters here because Tasic 2016 is ~12x deeper than the Allen SSv4 sets.

It also sweeps the per-panel read threshold, reporting the ratio and the surviving cell count at each, which
answers directly how deep a library has to be for this question and how many cells each dataset can afford.
Panels are gated independently: requiring all three at once dropped whole neuronal subclasses to 0%.

Usage:
    HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 channel_noise_floor.py [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

RNG    = np.random.default_rng(0)
PSEUDO = float(os.environ.get('PSEUDO', 0.5))
OUT    = sys.argv[1] if len(sys.argv) > 1 else '.'
SWEEP  = [20, 50, 100, 200, 300, 500, 1000, 2000]
N_SPLIT = int(os.environ.get('N_SPLIT', 3))      # independent split replicates per cell, averaged

PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}


def clr(A, pseudo=PSEUDO):
    A = np.asarray(A, float) + pseudo
    A = A / A.sum(axis=1, keepdims=True)
    L = np.log(A)
    return L - L.mean(axis=1, keepdims=True)


def panel_counts(d, panel):
    """cells x variants integer count matrix for one panel, plus per-cell group labels."""
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    C = V.reindex(columns=PANEL_V[panel]).fillna(0.0)
    C = C.round().astype(int)                       # salmon reads are fractional; splitting needs integers
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    return C, meta.reindex(C.index)


def noise_floor(C):
    """Per-cell technical-replicate distance, halved to full-depth error (see module docstring)."""
    A = C.to_numpy()
    out = np.zeros(len(A))
    for _ in range(N_SPLIT):
        h1 = RNG.binomial(A, 0.5)
        h2 = A - h1
        out += np.linalg.norm(clr(h1) - clr(h2), axis=1)
    return (out / N_SPLIT) / 2.0


def within_spread(X, labels):
    """Distance from each cell to its own group centroid."""
    out = np.zeros(len(X))
    labels = np.asarray(labels)
    for g in pd.unique(labels):
        i = np.where(labels == g)[0]
        if len(i) < 2:
            out[i] = np.nan; continue
        out[i] = np.linalg.norm(X[i] - X[i].mean(axis=0), axis=1)
    return out


def analyse(name, d):
    rows, sweep_rows = [], []
    for panel in PANEL_V:
        C, meta = panel_counts(d, panel)
        tot = C.sum(axis=1).to_numpy()
        for thr in SWEEP:
            keep = tot >= thr
            if keep.sum() < 20:
                sweep_rows.append(dict(dataset=name, panel=panel, threshold=thr,
                                       n_cells=int(keep.sum()), noise=np.nan, spread=np.nan, ratio=np.nan))
                continue
            Ck = C[keep]
            nf = np.median(noise_floor(Ck))
            X = clr(Ck.to_numpy())
            sp = np.nanmedian(within_spread(X, meta.loc[Ck.index, 'subclass']))
            sweep_rows.append(dict(dataset=name, panel=panel, threshold=thr, n_cells=int(keep.sum()),
                                   noise=nf, spread=sp, ratio=sp / nf if nf > 0 else np.nan))
        # per-subclass detail at a usable threshold
        thr = int(os.environ.get('DETAIL_THR', 200))
        keep = tot >= thr
        if keep.sum() >= 20:
            Ck = C[keep]; mk = meta.loc[Ck.index]
            nf = noise_floor(Ck); X = clr(Ck.to_numpy())
            sp = within_spread(X, mk['subclass'].to_numpy())
            for sc, idx in pd.Series(range(len(mk)), index=mk['subclass'].to_numpy()).groupby(level=0):
                i = idx.to_numpy()
                if len(i) < 10:
                    continue
                rows.append(dict(dataset=name, panel=panel, subclass=sc, n_cells=len(i),
                                 median_total_reads=float(np.median(Ck.sum(axis=1).to_numpy()[i])),
                                 noise=float(np.median(nf[i])), spread=float(np.nanmedian(sp[i])),
                                 ratio=float(np.nanmedian(sp[i]) / np.median(nf[i])) if np.median(nf[i]) > 0 else np.nan))
    return pd.DataFrame(rows), pd.DataFrame(sweep_rows)


def prep_generic(path, label):
    v = _read_percell(path)
    v['vclass'] = v['transcript'].map(final_map)
    v = v.dropna(subset=['vclass']).copy()
    if 'biotype' in v.columns:
        v = v[~((v['biotype'] == 'retained_intron')
                & v['vclass'].str.contains('other|fragment', case=False))].copy()
    v['panel'] = v['gene'].map(PANEL_OF)
    ren = {'subclass_label': 'subclass', 'cluster_label': 'cluster'}
    v = v.rename(columns={k: x for k, x in ren.items() if k in v.columns})
    if 'subclass' not in v.columns:
        sys.exit(f'{label}: no subclass column; found {list(v.columns)}')
    if 'cluster' not in v.columns:
        v['cluster'] = v['subclass']
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
        try:
            ds['A-VISp'] = prep_generic(vp, 'A-VISp')
        except SystemExit as e:
            print(f'[warn] {e}')
    if not ds:
        sys.exit('no datasets loaded')

    detail, sweep = [], []
    for name, d in ds.items():
        print(f'\n=== {name}: {d["sample_name"].nunique():,} cells')
        a, b = analyse(name, d)
        detail.append(a); sweep.append(b)
    D = pd.concat(detail, ignore_index=True); S = pd.concat(sweep, ignore_index=True)

    print('\n' + '=' * 96)
    print('READ-DEPTH SWEEP:  spread/noise ~1 means the scatter is sampling noise, not biology')
    print('=' * 96)
    for name in ds:
        print(f'\n{name}')
        for panel in PANEL_V:
            s = S[(S.dataset == name) & (S.panel == panel)]
            if s.empty:
                continue
            cells = '  '.join(f'{int(r.threshold)}:{r.n_cells}' for _, r in s.iterrows())
            rat = '  '.join(f'{int(r.threshold)}:{r.ratio:.2f}' if np.isfinite(r.ratio) else f'{int(r.threshold)}:--'
                            for _, r in s.iterrows())
            print(f'  {panel:6} cells  {cells}')
            print(f'  {"":6} ratio  {rat}')

    thr = int(os.environ.get('DETAIL_THR', 200))
    print('\n' + '=' * 96)
    print(f'PER-SUBCLASS at >= {thr} reads/panel  (ratio = within-subclass spread / noise floor)')
    print('=' * 96)
    if not D.empty:
        piv = D.pivot_table(index=['dataset', 'subclass'], columns='panel', values='ratio')
        n = D.pivot_table(index=['dataset', 'subclass'], columns='panel', values='n_cells')
        print(piv.round(2).join(n.round(0), rsuffix='_n').to_string())
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    D.to_csv(f'{OUT}/tables/channel_noise_by_subclass.csv', index=False)
    S.to_csv(f'{OUT}/tables/channel_noise_sweep.csv', index=False)
    print(f'\nwrote {OUT}/tables/channel_noise_{{by_subclass,sweep}}.csv')
    print('\nCalibration (synthetic, 17-variant panel): pure sampling noise gives ratio 0.99 at depth 200 and')
    print('0.98 at depth 2000 -- the metric is depth-independent under the null. But real MILD heterogeneity')
    print('reads only 1.29 at depth 200 against 2.63 at depth 2000, so sensitivity does depend on depth:')
    print('  ratio > ~1.5   real biological heterogeneity, trustworthy at any depth')
    print('  ratio ~ 1.0-1.3 at low depth   INCONCLUSIVE, not evidence of homogeneity -- mild real variance')
    print('                 is masked here, so no grouping verdict should be drawn')
    print('  ratio ~ 1.0 at high depth      genuinely homogeneous within that group')
    print('Read the cell counts beside the ratios: a threshold that reaches a usable ratio with 30 cells left')
    print('has not made the question answerable, it has just made it unanswerable a different way.')


if __name__ == '__main__':
    main()
