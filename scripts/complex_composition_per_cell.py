#!/usr/bin/env python3
"""What would a patch pipette find? Per-cell Kv4-complex composition, as recordable categories.

Translates the per-cell isoform data into the categories an electrophysiologist records in:

  1  no Kv4 channel at all            (Kcnd1 + Kcnd2 + Kcnd3 undetected)
  2  Kv4 present but no Kv4.2         (Kcnd2 undetected, Kcnd1 or Kcnd3 present)
  3  Kv4 alone, no KChIP, no DPP      -> fast-inactivating, ER-retained, little surface current
  4  Kv4 + DPP, no KChIP              -> hyperpolarised V1/2, fast kinetics, poor surface expression
  5  Kv4 + KChIP, no DPP              -> surface-expressed, slowed inactivation, depolarised V1/2
  6  Kv4 + KChIP + DPP                -> the full ternary complex

plus the N-type question (DPP6a and DPP10a carry the conserved MNQTA N-terminus that confers fast
voltage-independent N-type inactivation; the other variants do not), and the per-cell distribution of DPP6K
as a share of all DPP.

THE DOMINANT CAVEAT IS DETECTION, NOT BIOLOGY. Every "absent" call is really "undetected", and detection
depends on library depth. If the cells called KChIP-negative are simply the shallow ones, the category
fractions are a sequencing artefact. So for every zero-dependent call this reports the median depth of the
negative cells against the positive cells: a ratio near 1.0 means the call is not depth-driven, and a ratio
well below 1.0 means it is. It also sweeps the read threshold, since a conclusion that changes between
">=1 read" and ">=5 reads" is not a conclusion.

TWO CAVEATS IT CANNOT ADDRESS, both of which make these fractions UNDERESTIMATE the protein-level presence:
  - transcriptional bursting. A gene off at this instant may have been on an hour ago, and the companion
    Oct 30 analysis showed mature Kcnd2 mRNA is a long-lived accumulated pool fed by bursts, so a cell with
    no detectable mRNA can still hold channel protein.
  - mRNA is not protein. Translation, trafficking, assembly stoichiometry and protein half-life all
    intervene, and the complex assembles somatically, which integrates composition over time and should
    SMOOTH what a recording sees relative to what the transcript says.
Read every fraction below as "what the mRNA snapshot says", not "what the pipette will find".

Usage:  HIP_DATA=... CTX_DATA=... VISP_DATA=... python3 complex_composition_per_cell.py [out_dir]
        TYPES='L6a Sla,L6a Mgp,L2/3 Ptgs2' python3 complex_composition_per_cell.py
"""
import os, sys, itertools
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

THRESHOLDS = [int(x) for x in os.environ.get('THRESHOLDS', '1,3,5,10').split(',')]
MAIN_THR   = int(os.environ.get('MAIN_THR', 5))
MIN_CELLS  = int(os.environ.get('MIN_CELLS', 20))
OUT        = sys.argv[1] if len(sys.argv) > 1 else '.'

KV4_GENES   = ['Kcnd1', 'Kcnd2', 'Kcnd3']
KCHIP_GENES = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']
DPP_GENES   = ['Dpp6', 'Dpp10']
NTYPE_V     = ['DPP6a', 'DPP10a']      # conserved MNQTA N-terminus -> N-type inactivation
DPP_V       = ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d']


def per_cell(d):
    """Gene-level reads per cell, plus the DPP variant breakdown and a depth measure."""
    g = d.groupby(['sample_name', 'gene'])['reads'].sum().unstack(fill_value=0.0)
    v = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'), cluster=('cluster', 'first'))
    depth = None
    for c in ('num_processed', 'total_reads', 'n_reads'):
        if c in d.columns:
            depth = d.groupby('sample_name')[c].first()
            break
    if depth is None:                      # fall back to total Kv4-complex reads as a relative depth proxy
        depth = d.groupby('sample_name')['reads'].sum()
    X = pd.DataFrame(index=g.index)
    X['kv4'] = g.reindex(columns=KV4_GENES).fillna(0).sum(axis=1)
    X['kcnd2'] = g.reindex(columns=['Kcnd2']).fillna(0).sum(axis=1)
    X['kchip'] = g.reindex(columns=KCHIP_GENES).fillna(0).sum(axis=1)
    X['dpp'] = g.reindex(columns=DPP_GENES).fillna(0).sum(axis=1)
    V = v.reindex(columns=DPP_V).fillna(0.0)
    X['dpp_variant_total'] = V.sum(axis=1)
    X['ntype'] = V.reindex(columns=NTYPE_V).fillna(0).sum(axis=1)
    X['dpp6k'] = V.reindex(columns=['DPP6K']).fillna(0).sum(axis=1)
    X['depth'] = depth.reindex(X.index)
    return X.join(meta)


def categorise(X, thr):
    kv4 = X.kv4 >= thr
    k2 = X.kcnd2 >= thr
    kc = X.kchip >= thr
    dp = X.dpp >= thr
    cat = pd.Series(index=X.index, dtype=object)
    cat[~kv4] = '1 no Kv4'
    cat[kv4 & ~k2] = '2 Kv4 but no Kv4.2'
    cat[kv4 & k2 & ~kc & ~dp] = '3 Kv4 only'
    cat[kv4 & k2 & ~kc & dp] = '4 Kv4+DPP, no KChIP'
    cat[kv4 & k2 & kc & ~dp] = '5 Kv4+KChIP, no DPP'
    cat[kv4 & k2 & kc & dp] = '6 Kv4+KChIP+DPP'
    return cat


CATS = ['1 no Kv4', '2 Kv4 but no Kv4.2', '3 Kv4 only', '4 Kv4+DPP, no KChIP',
        '5 Kv4+KChIP, no DPP', '6 Kv4+KChIP+DPP']


def depth_check(X, mask_neg, label):
    """Are the cells called negative simply the shallow ones?"""
    if mask_neg.sum() == 0 or (~mask_neg).sum() == 0:
        return dict(call=label, n_neg=int(mask_neg.sum()), depth_ratio=np.nan, p=np.nan)
    dn, dp_ = X.depth[mask_neg], X.depth[~mask_neg]
    p = stats.mannwhitneyu(dn, dp_).pvalue if len(dn) > 2 and len(dp_) > 2 else np.nan
    return dict(call=label, n_neg=int(mask_neg.sum()),
                median_depth_neg=float(np.median(dn)), median_depth_pos=float(np.median(dp_)),
                depth_ratio=float(np.median(dn) / max(np.median(dp_), 1e-9)), p=p)


def analyse(name, d, types=None):
    X = per_cell(d)
    label = 'cluster' if X['cluster'].nunique() > X['subclass'].nunique() else 'subclass'
    rows, dist, checks = [], [], []
    for ct, sub in X.groupby(label):
        if len(sub) < MIN_CELLS or (types and ct not in types):
            continue
        for thr in THRESHOLDS:
            c = categorise(sub, thr).value_counts(normalize=True)
            r = dict(dataset=name, cell_type=ct, n_cells=len(sub), threshold=thr)
            for k in CATS:
                r[k] = float(c.get(k, 0.0)) * 100
            # N-type: share of DPP reads carrying the MNQTA N-terminus, and how many cells have any
            ok = sub.dpp_variant_total >= thr
            nt = (sub.ntype[ok] / sub.dpp_variant_total[ok]) if ok.sum() else pd.Series(dtype=float)
            r['n_cells_dpp_measurable'] = int(ok.sum())
            r['pct_cells_any_Ntype'] = float((nt > 0.05).mean() * 100) if len(nt) else np.nan
            r['Ntype_frac_median'] = float(nt.median()) if len(nt) else np.nan
            r['Ntype_frac_iqr'] = float(nt.quantile(.75) - nt.quantile(.25)) if len(nt) else np.nan
            k6 = (sub.dpp6k[ok] / sub.dpp_variant_total[ok]) if ok.sum() else pd.Series(dtype=float)
            r['DPP6K_median'] = float(k6.median()) if len(k6) else np.nan
            r['DPP6K_mean'] = float(k6.mean()) if len(k6) else np.nan
            r['DPP6K_pct_cells_below5'] = float((k6 < .05).mean() * 100) if len(k6) else np.nan
            r['DPP6K_pct_cells_above50'] = float((k6 > .50).mean() * 100) if len(k6) else np.nan
            rows.append(r)
            if thr == MAIN_THR and ok.sum():
                for cell, val in nt.items():
                    dist.append(dict(dataset=name, cell_type=ct, cell=cell,
                                     ntype_frac=val, dpp6k_frac=k6.get(cell, np.nan)))
        # detection diagnostics at the main threshold
        for lab, m in [('no Kv4', sub.kv4 < MAIN_THR), ('no Kv4.2', sub.kcnd2 < MAIN_THR),
                       ('no KChIP', sub.kchip < MAIN_THR), ('no DPP', sub.dpp < MAIN_THR)]:
            c = depth_check(sub, m, lab); c.update(dataset=name, cell_type=ct, n_cells=len(sub))
            checks.append(c)
    return pd.DataFrame(rows), pd.DataFrame(dist), pd.DataFrame(checks)


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
    for lab, loader in [('T16', load_ctx), ('A-HPF', load_hip)]:
        try:
            ds[lab] = loader()
        except Exception as e:
            print(f'[warn] {lab}: {e}')
    vp = os.environ.get('VISP_DATA') or os.environ.get('VISP_PKL')
    if vp and os.path.exists(vp):
        ds['A-VISp'] = prep_generic(vp)
    if not ds:
        sys.exit('no datasets loaded')

    R, D, C = [], [], []
    for n, d in ds.items():
        a, b, c = analyse(n, d, types)
        R.append(a); D.append(b); C.append(c)
    R = pd.concat(R, ignore_index=True); D = pd.concat(D, ignore_index=True); C = pd.concat(C, ignore_index=True)
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/complex_categories.csv', index=False)
    D.to_csv(f'{OUT}/tables/complex_percell_fracs.csv', index=False)
    C.to_csv(f'{OUT}/tables/complex_detection_checks.csv', index=False)

    M = R[R.threshold == MAIN_THR].copy()
    print('=' * 118)
    print(f'WHAT THE mRNA SAYS A PIPETTE WOULD FIND  (>= {MAIN_THR} reads per gene group to call "present")')
    print('  % of cells in each category. Read the detection table below before trusting any "no X" column.')
    print('=' * 118)
    show = M.sort_values(['dataset', 'cell_type'])
    print(f'{"dataset":7} {"cell type":18} {"n":>4} ' + ' '.join(f'{c.split()[0]:>6}' for c in CATS)
          + '   legend: ' + '; '.join(CATS))
    for _, r in show.iterrows():
        print(f'{r.dataset:7} {str(r.cell_type)[:18]:18} {int(r.n_cells):4d} '
              + ' '.join(f'{r[c]:6.1f}' for c in CATS))

    print('\n' + '=' * 118)
    print('N-TYPE INACTIVATION (DPP6a + DPP10a share of DPP) and DPP6K')
    print('=' * 118)
    print(f'{"dataset":7} {"cell type":18} {"nDPP":>5} {"%cells N-type":>13} {"N-frac med":>11} '
          f'{"IQR":>6}  {"DPP6K med":>9} {"mean":>6} {"<5%":>6} {">50%":>6}')
    for _, r in show.iterrows():
        if not np.isfinite(r.get('Ntype_frac_median', np.nan)):
            continue
        print(f'{r.dataset:7} {str(r.cell_type)[:18]:18} {int(r.n_cells_dpp_measurable):5d} '
              f'{r.pct_cells_any_Ntype:13.1f} {r.Ntype_frac_median:11.3f} {r.Ntype_frac_iqr:6.3f}  '
              f'{r.DPP6K_median:9.3f} {r.DPP6K_mean:6.3f} {r.DPP6K_pct_cells_below5:6.1f} '
              f'{r.DPP6K_pct_cells_above50:6.1f}')

    print('\n' + '=' * 118)
    print('DETECTION CHECK: are the cells called negative just the shallow ones?')
    print('  depth_ratio = median library depth of negative cells / positive cells.')
    print('  ~1.0 means the call is not depth-driven.  <<1.0 means it is, and that category is an artefact.')
    print('=' * 118)
    agg = (C[C.n_neg > 0].groupby('call')
           .agg(types=('cell_type', 'nunique'), median_depth_ratio=('depth_ratio', 'median'),
                frac_types_depth_driven=('depth_ratio', lambda x: float((x < 0.7).mean()))))
    print(agg.round(3).to_string())

    # cross-type comparison: are the distributions type-specific or interchangeable?
    print('\n' + '=' * 118)
    print('ARE TYPES DISTINGUISHABLE?  pairwise tests on the per-cell N-type and DPP6K distributions')
    print('=' * 118)
    for ds_name, sub in D.groupby('dataset'):
        cts = [c for c, g in sub.groupby('cell_type') if len(g) >= MIN_CELLS]
        if len(cts) < 2:
            continue
        print(f'\n{ds_name}:')
        res = []
        for a, b in itertools.combinations(cts, 2):
            x = sub[sub.cell_type == a]; y = sub[sub.cell_type == b]
            for col in ('ntype_frac', 'dpp6k_frac'):
                xv, yv = x[col].dropna(), y[col].dropna()
                if len(xv) < 10 or len(yv) < 10:
                    continue
                ks = stats.ks_2samp(xv, yv)
                res.append(dict(a=a, b=b, measure=col, ks=ks.statistic, p=ks.pvalue,
                                med_a=xv.median(), med_b=yv.median()))
        if not res:
            continue
        P = pd.DataFrame(res)
        from statsmodels.stats.multitest import multipletests
        P['q'] = multipletests(P.p, method='fdr_bh')[1]
        sig = P[P.q < 0.05].sort_values('q')
        print(f'  {len(sig)} of {len(P)} pairwise comparisons differ at FDR 5%')
        for _, r in sig.head(12).iterrows():
            print(f'    {r.measure:11} {str(r.a)[:16]:16} vs {str(r.b)[:16]:16} '
                  f'D={r.ks:.2f} q={r.q:.1e}  medians {r.med_a:.3f} / {r.med_b:.3f}')
        P.to_csv(f'{OUT}/tables/complex_type_comparisons_{ds_name}.csv', index=False)

    print('\n' + '=' * 118)
    print('THRESHOLD SENSITIVITY: category 6 (full ternary complex), % of cells')
    print('=' * 118)
    piv = R.pivot_table(index=['dataset', 'cell_type'], columns='threshold', values='6 Kv4+KChIP+DPP')
    print(piv.round(1).to_string())
    print('\nA category that moves a lot across thresholds is a detection statement, not a biological one.')
    print(f'\nwrote {OUT}/tables/complex_{{categories,percell_fracs,detection_checks}}.csv')


if __name__ == '__main__':
    main()
