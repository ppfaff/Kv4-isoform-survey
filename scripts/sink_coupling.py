#!/usr/bin/env python3
"""Is the per-cell variant bimodality biology, or cell-to-cell variation in RNA integrity?

simulate_variant_calibration.py established the mechanism. With full-length coverage salmon assigns the DPP
variants perfectly (100% on the diagonal for all seven). Under 3' bias it collapses into two DEFAULT BINS,
one per gene: DPP6K and DPP6L flow into DPP6S (which absorbs 129-166% of other variants' reads), and DPP10a
flows into DPP10d (77%). Applying that confusion matrix to the 2012 qRT-PCR composition of cerebellar granule
cells predicts DPP6S at 50-57% against the 53.7% actually observed in the V404M bulk, so the artefact is
real and quantitatively sufficient.

That leaves one way for the single-cell bimodality to be artefactual too. RNA integrity varies from cell to
cell. A cell with intact RNA gives full-length coverage and reports DPP6K correctly; a degraded cell is
effectively 3'-biased and reports the same molecule as DPP6S. Cells would then split into "expresses DPP6K"
and "doesn't" purely by RNA quality -- the exact all-or-none pattern variant_distribution_shape.py found,
with the sink variant as the dominant one.

THE TEST. Dpp6 and Dpp10 have INDEPENDENT sinks. A shared technical cause -- that cell's coverage profile --
pushes both at once, so sink usage should be CORRELATED ACROSS CELLS: cells high in DPP6S/Dpp6 should also
be high in DPP10d/Dpp10. Biology offers no reason for a cell's choice between Dpp6 N-terminal variants to be
coupled to its choice between Dpp10 variants; they are different genes on different chromosomes with
different promoters.

    strong positive correlation -> shared technical driver; the bimodality is an RNA-quality artefact
    no correlation             -> the two genes vary independently; a technical sink cannot explain it

Controls, because a raw correlation is not enough:
  - compositional null. Sink shares are bounded fractions, so some correlation arises by construction;
    the null resamples each cell at its own depth from the cell type's mean composition.
  - a NON-sink pair (DPP6a vs DPP10c) is reported alongside, but it is NOT a negative control and the
    verdict does not use it. Validation showed why: compositional closure means that when a sink's share
    rises, every other variant of that gene falls, so DPP6a and DPP10c move together almost as strongly as
    the sinks do under a shared driver (+0.727 against +0.749). Closure acts WITHIN a gene, so neither pair
    correlates across genes without a shared driver; the control therefore adds no discrimination. It is
    kept only because the sinks coupling while the control does not would be strange and worth seeing.

WHAT COUPLING DOES AND DOES NOT PROVE. Coupling establishes a shared CELL-LEVEL driver acting on both genes.
That could be RNA integrity, which is the technical hypothesis the simulation makes concrete -- or a
splicing factor coordinating both genes, which would be biology. To separate those, correlate the sink share
against a DIRECT per-cell 3'-bias metric (gene-body coverage from an alignment); a technical driver predicts
that correlation, a splicing factor does not. Absence of coupling is the cleaner result: it rules the
technical mechanism out, since one coverage profile per cell cannot move two genes independently.
  - the test runs WITHIN cell type, so it cannot be driven by type composition differences.

Usage:  HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 sink_coupling.py [out_dir]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

RNG       = np.random.default_rng(0)
MIN_READS = int(os.environ.get('MIN_READS', 50))     # per GENE, per cell
MIN_CELLS = int(os.environ.get('MIN_CELLS', 25))
N_NULL    = int(os.environ.get('N_NULL', 200))
OUT       = sys.argv[1] if len(sys.argv) > 1 else '.'

DPP6  = ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L']
DPP10 = ['DPP10a', 'DPP10c', 'DPP10d']
PAIRS = [('sink', 'DPP6S', 'DPP10d'),        # the two default bins
         ('control', 'DPP6a', 'DPP10c')]     # robust at every tau in the simulation


def shares(d):
    """Per cell: each variant's share of its OWN gene, so Dpp6 and Dpp10 are measured independently."""
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    out = {}
    for gene, vs in [('Dpp6', DPP6), ('Dpp10', DPP10)]:
        C = V.reindex(columns=vs).fillna(0.0)
        tot = C.sum(axis=1)
        out[gene] = (C, tot)
    return out, meta.reindex(V.index)


def null_rho(c1, t1, c2, t2, i1, i2, n_null=N_NULL):
    """Correlation expected from the compositional constraint alone, at each cell's own depth."""
    p1 = c1.sum(axis=0) / max(c1.sum(), 1e-9)
    p2 = c2.sum(axis=0) / max(c2.sum(), 1e-9)
    out = []
    for _ in range(max(1, n_null // 10)):
        s1 = np.array([RNG.multinomial(int(n), p1) if n > 0 else np.zeros(len(p1)) for n in t1], float)
        s2 = np.array([RNG.multinomial(int(n), p2) if n > 0 else np.zeros(len(p2)) for n in t2], float)
        f1 = s1[:, i1] / np.maximum(s1.sum(axis=1), 1e-9)
        f2 = s2[:, i2] / np.maximum(s2.sum(axis=1), 1e-9)
        if f1.std() > 0 and f2.std() > 0:
            out.append(stats.spearmanr(f1, f2).statistic)
    return (float(np.nanmean(out)), float(np.nanpercentile(out, 95))) if out else (np.nan, np.nan)


def analyse(name, d):
    g, meta = shares(d)
    (C6, T6), (C10, T10) = g['Dpp6'], g['Dpp10']
    keep = (T6 >= MIN_READS) & (T10 >= MIN_READS)
    if keep.sum() < MIN_CELLS:
        return pd.DataFrame()
    C6, T6, C10, T10, mk = C6[keep], T6[keep], C10[keep], T10[keep], meta[keep]
    label = 'cluster' if mk['cluster'].nunique() > mk['subclass'].nunique() else 'subclass'
    rows = []
    for ct, grp in mk.groupby(label).groups.items():
        ii = mk.index.get_indexer(grp)
        if len(ii) < MIN_CELLS:
            continue
        a6, b6 = C6.to_numpy()[ii], T6.to_numpy()[ii]
        a10, b10 = C10.to_numpy()[ii], T10.to_numpy()[ii]
        for kind, v6, v10 in PAIRS:
            i1, i2 = DPP6.index(v6), DPP10.index(v10)
            f1 = a6[:, i1] / np.maximum(b6, 1e-9)
            f2 = a10[:, i2] / np.maximum(b10, 1e-9)
            if f1.std() == 0 or f2.std() == 0:
                continue
            rho = stats.spearmanr(f1, f2)
            nm, n95 = null_rho(a6, b6, a10, b10, i1, i2)
            rows.append(dict(dataset=name, cell_type=ct, kind=kind, pair=f'{v6} vs {v10}',
                             n_cells=len(ii), rho=rho.statistic, p=rho.pvalue,
                             null_rho=nm, null_rho_p95=n95,
                             rho_beyond_null=rho.statistic - n95))
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

    R = pd.concat([analyse(n, d) for n, d in ds.items()], ignore_index=True)
    if R.empty:
        sys.exit(f'no cell type had {MIN_CELLS}+ cells with {MIN_READS}+ reads on BOTH genes')
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/sink_coupling.csv', index=False)

    print('=' * 104)
    print('DO THE TWO SINKS MOVE TOGETHER WITHIN CELLS?')
    print('  sink pair    = DPP6S vs DPP10d, the two default bins under 3\' bias')
    print('  control pair = DPP6a vs DPP10c, robust at every tau in the simulation')
    print('=' * 104)
    print(f'{"dataset":7} {"cell type":18} {"kind":8} {"n":>4} {"rho":>7} {"null95":>7} '
          f'{"beyond":>7}  p')
    for _, r in R.sort_values(['dataset', 'cell_type', 'kind']).iterrows():
        print(f'{r.dataset:7} {str(r.cell_type)[:18]:18} {r.kind:8} {int(r.n_cells):4d} '
              f'{r.rho:7.3f} {r.null_rho_p95:7.3f} {r.rho_beyond_null:7.3f}  {r.p:.2g}')

    print('\n' + '=' * 104)
    print('SUMMARY')
    print('=' * 104)
    for kind in ('sink', 'control'):
        S = R[R.kind == kind]
        if S.empty:
            continue
        print(f'{kind:8} median rho {S.rho.median():+.3f}   median beyond null {S.rho_beyond_null.median():+.3f}'
              f'   {(S.rho_beyond_null > 0.1).mean():.0%} of types exceed their null')
    s, c = R[R.kind == 'sink'], R[R.kind == 'control']
    if not s.empty and not c.empty:
        print()
        if s.rho_beyond_null.median() > 0.15:
            print('VERDICT: the sinks move together across cells, beyond the compositional null.')
            print('  A shared CELL-LEVEL driver is shifting both genes at once. RNA integrity is the')
            print('  candidate the simulation makes concrete, and on that reading the per-cell variant')
            print('  bimodality is an artefact and the ephys predictions built on it do not follow.')
            print('  Before concluding that, rule out the alternative shared driver -- a splicing factor')
            print('  acting on both genes, which is biology: correlate the sink share against a direct')
            print('  per-cell 3\'-bias metric from an alignment. Technical predicts it; splicing does not.')
        elif s.rho_beyond_null.median() < 0.1:
            print('VERDICT: no coupling beyond the compositional null.')
            print('  Dpp6 and Dpp10 variant usage varies independently across cells, which a shared coverage')
            print('  artefact cannot produce. The bimodality survives this test.')
        else:
            print('VERDICT: ambiguous -- partial coupling, or the control pair behaves similarly.')
            print('  Compare the sink and control rows directly before concluding anything.')
    print(f'\nwrote {OUT}/tables/sink_coupling.csv')


if __name__ == '__main__':
    main()
