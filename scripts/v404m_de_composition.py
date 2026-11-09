#!/usr/bin/env python3
"""Per-tissue differential expression with tissue-composition covariates.

Motivation (2026-11-09): the V404M 'cortex' DE is dominated by olfactory-bulb tissue in some mutant blocks
(Omp, S100a5, Stoml3, Doc2g, Tbx21, Th, Slc6a3 all up) and the hippocampal DE's large-fold tail is choroid
plexus (Ttr, Kcnj13, Slc4a5, Cldn3, Otop1, Atp6v1b1), with a small striatal imbalance in cortex on top.
Here every library gets composition indices (log2 mean signature CPM relative to Snap25) and the DE is
rerun with and without them in the design, per tissue:

    without:  ~ sex + genotype
    with:     ~ sex + genotype + <composition indices for that tissue>

Engine: pydeseq2 if importable (matches the original analysis), otherwise an OLS model on log2(CPM+1)
with a per-gene t-test on the genotype coefficient and Benjamini-Hochberg FDR (limma-style, no shrinkage).
Both engines run the same two designs so the before/after comparison is like for like.

Runs where the counts matrix lives (tables/v404m_star_counts.csv.gz); kv4py env is enough for the OLS fallback.
Usage: python3 scripts/v404m_de_composition.py [TABLES_DIR] [OUT_DIR] [--engine ols|pydeseq2]
Writes OUT/tables/v404m_composition_indices.csv, v404m_de_adjusted_<tissue>.csv, v404m_de_adjusted_summary.csv,
       v404m_de_adjusted_keygenes.csv
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v404m_panels import load_counts, load_meta  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('tables', nargs='?', default='tables')
ap.add_argument('out', nargs='?', default='.')
ap.add_argument('--engine', default='auto', choices=['auto', 'ols', 'pydeseq2'])
ap.add_argument('--min-cpm', type=float, default=1.0)
args = ap.parse_args()
os.makedirs(f'{args.out}/tables', exist_ok=True)

SIG = {
    'striatum': ['Gpr88', 'Adora2a', 'Ppp1r1b', 'Penk', 'Drd1', 'Drd2', 'Rgs9', 'Pde10a'],
    'olfactory_bulb': ['Omp', 'S100a5', 'Stoml3', 'Doc2g', 'Cdhr1', 'Tbx21', 'Cnga2', 'Gng13', 'Th', 'Slc6a3'],
    'choroid_plexus': ['Ttr', 'Kcnj13', 'Folr1', 'Clic6', 'Slc4a5', 'Kl', 'Otx2', 'Prlr', 'Cldn3', 'Atp6v1b1'],
    'cortex': ['Slc17a7', 'Satb2', 'Cux2', 'Rorb', 'Fezf2', 'Foxp2', 'Syt6'],
    'hippocampus': ['Prox1', 'C1ql2', 'Dsp', 'Zbtb20', 'Amigo2'],
    'thalamus': ['Tcf7l2', 'Prkcd', 'Gbx2', 'Slc17a6', 'Ntng1'],
    'cerebellum': ['Pcp2', 'Gabra6', 'Cbln1', 'Zic1'],
}
# which contaminants can plausibly vary within each dissection
COVARS = {
    'cortex': ['striatum', 'olfactory_bulb'],
    'hippocampus': ['choroid_plexus', 'cortex'],
    'thalamus': ['striatum', 'cortex', 'hippocampus'],
    'cerebellum': ['choroid_plexus'],
}
KEY = {
    'IEG': ['Fos', 'Fosb', 'Egr1', 'Egr2', 'Egr4', 'Npas4', 'Arc', 'Fosl2', 'Junb', 'Nr4a1'],
    'DG activity program': ['Penk', 'Pdyn', 'Drd1', 'Npy', 'Ecel1', 'Igfbpl1', 'Rprm', 'Cck', 'Rasd1', 'Bdnf'],
    'reactive astrocyte': ['Gfap', 'Serpina3n', 'Timp1', 'Ptx3', 'Vim', 'Osmr', 'Cd44', 'C3', 'Aqp4', 'Kcnj10'],
    'OPC': ['Pdgfra', 'Cspg4', 'Sox10', 'Olig2', 'Gpr17'],
    'Kv4 complex': ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10'],
    'olfactory bulb': ['Omp', 'S100a5', 'Stoml3', 'Tbx21'],
    'choroid plexus': ['Ttr', 'Kcnj13', 'Cldn3', 'Otop1'],
}

meta = load_meta(args.tables)
C, src = load_counts(args.tables)
common = [s for s in C.columns if s in meta.index]
C, meta = C[common].round().astype(int), meta.loc[common].copy()
meta['genotype'] = np.where(meta.genotype.str.startswith('mut'), 'mutant', 'wildtype')
cpm = C / C.sum() * 1e6
snap = np.log2(cpm.loc['Snap25'] + 1)
idx = pd.DataFrame({k: np.log2(cpm.loc[[g for g in v if g in cpm.index]] + 1).mean() - snap for k, v in SIG.items()})
idx.join(meta).to_csv(f'{args.out}/tables/v404m_composition_indices.csv')
print(f'counts: {src}; {len(common)} libraries')

engine = args.engine
if engine == 'auto':
    try:
        import pydeseq2  # noqa: F401
        engine = 'pydeseq2'
    except ImportError:
        engine = 'ols'
print('engine:', engine)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); r = np.empty(n)
    r[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(r, 1.0)


def run_ols(sub_counts, design):
    """OLS on log2(CPM+1); returns DataFrame gene, baseMean, log2FoldChange, pvalue, padj for the genotype term."""
    Y = np.log2(sub_counts / sub_counts.sum() * 1e6 + 1).values  # genes x samples
    X = design.values.astype(float); n, k = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    B = Y @ X @ XtXi.T                       # genes x k
    resid = Y - B @ X.T
    s2 = (resid ** 2).sum(axis=1) / (n - k)
    j = list(design.columns).index('genotype_mutant')
    se = np.sqrt(s2 * XtXi[j, j])
    t = B[:, j] / se
    p = 2 * stats.t.sf(np.abs(t), n - k)
    out = pd.DataFrame({'gene': sub_counts.index, 'baseMean': sub_counts.mean(axis=1).values,
                        'log2FoldChange': B[:, j], 'pvalue': p})
    out['padj'] = bh(out.pvalue)
    return out


def run_pydeseq2(sub_counts, design_df, factors):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    dds = DeseqDataSet(counts=sub_counts.T, metadata=design_df, design_factors=factors,
                       continuous_factors=[f for f in factors if f not in ('sex', 'genotype')], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=['genotype', 'mutant', 'wildtype'], quiet=True); st.summary()
    r = st.results_df.reset_index().rename(columns={'index': 'gene'})
    return r[['gene', 'baseMean', 'log2FoldChange', 'pvalue', 'padj']]


summary, key_rows = [], []
for tissue in ['cortex', 'hippocampus', 'thalamus', 'cerebellum']:
    samp = meta.index[meta.tissue == tissue]
    sub = C[samp]
    keep = (cpm[samp].mean(axis=1) >= args.min_cpm)
    sub = sub[keep]
    md = meta.loc[samp, ['sex', 'genotype']].join(idx.loc[samp, COVARS[tissue]])
    for cv in COVARS[tissue]:                           # center the covariates
        md[cv] = md[cv] - md[cv].mean()
    res = {}
    for label, covs in [('without', []), ('with', COVARS[tissue])]:
        if engine == 'pydeseq2':
            res[label] = run_pydeseq2(sub, md, ['sex', 'genotype'] + covs)
        else:
            X = pd.DataFrame({'intercept': 1.0, 'sex_male': (md.sex == 'male').astype(float),
                              'genotype_mutant': (md.genotype == 'mutant').astype(float)}, index=samp)
            for cv in covs:
                X[cv] = md[cv]
            res[label] = run_ols(sub, X)
    a, b = res['without'].set_index('gene'), res['with'].set_index('gene')
    both = a.join(b, lsuffix='_without', rsuffix='_with')
    both['tissue'] = tissue
    both.to_csv(f'{args.out}/tables/v404m_de_adjusted_{tissue}.csv')
    for thr in (0.05, 0.1):
        sa, sb = set(a.index[a.padj < thr]), set(b.index[b.padj < thr])
        summary.append(dict(tissue=tissue, fdr=thr, n_without=len(sa), n_with=len(sb), n_both=len(sa & sb),
                            lost=len(sa - sb), gained=len(sb - sa)))
    # covariate genotype imbalance, for the record (same on both FDR rows)
    for cv in COVARS[tissue]:
        v = idx.loc[samp, cv]; m = meta.loc[samp, 'genotype'] == 'mutant'
        for row in summary[-2:]:
            row[f'{cv}_mut_minus_wt'] = round(v[m].mean() - v[~m].mean(), 3)
    for grp, genes in KEY.items():
        for g in genes:
            if g in both.index:
                r = both.loc[g]
                key_rows.append(dict(tissue=tissue, group=grp, gene=g,
                                     lfc_without=r.log2FoldChange_without, padj_without=r.padj_without,
                                     lfc_with=r.log2FoldChange_with, padj_with=r.padj_with))

S = pd.DataFrame(summary); S.to_csv(f'{args.out}/tables/v404m_de_adjusted_summary.csv', index=False)
K = pd.DataFrame(key_rows); K.to_csv(f'{args.out}/tables/v404m_de_adjusted_keygenes.csv', index=False)
pd.set_option('display.width', 220)
print('\nDE gene counts, genotype term, before vs after composition covariates:')
print(S.to_string(index=False))
print('\nkey genes (lfc / padj without -> with covariates):')
K['without'] = K.apply(lambda r: f'{r.lfc_without:+.2f} ({r.padj_without:.3f})', axis=1)
K['with'] = K.apply(lambda r: f'{r.lfc_with:+.2f} ({r.padj_with:.3f})', axis=1)
for grp in KEY:
    print(f'\n-- {grp}')
    print(K[K.group == grp].pivot(index='gene', columns='tissue', values=['without', 'with'])
          .reindex(columns=pd.MultiIndex.from_product([['without', 'with'], ['cortex', 'hippocampus', 'thalamus', 'cerebellum']]))
          .to_string())
