#!/usr/bin/env python3
"""Does Kv4.2-V404M expression activate an integrated stress response or cell death?

Differential expression of V404M vs wild-type littermates within each tissue (hippocampus, cortex, thalamus,
cerebellum), then two complementary panel tests per pathway:

  competitive  -- are the panel's Wald statistics shifted relative to all other expressed genes? (Mann-Whitney).
                  This is the conservative test: it asks whether the pathway stands out against the rest of the
                  transcriptome, so a global shift in the data cannot create a positive.
  self-contained -- is the per-sample mean z-score of the panel genes different between genotypes? This is more
                  sensitive but will also fire if the whole transcriptome drifts, so the two are read together.

Both are reported with Benjamini-Hochberg FDR across panels within a tissue. A pathway is only called activated
when the competitive test agrees with the self-contained one and the direction is consistent across tissues or
the effect is large in the tissue where the phenotype lives.

Usage:  python3 scripts/v404m_stress.py [TABLES_DIR] [OUT_DIR]
Outputs: tables/v404m_de_<tissue>.csv, tables/v404m_panel_tests.csv, tables/v404m_panel_scores.csv,
         figures/v404m_panels.png, figures/v404m_stress_genes.png, figures/v404m_qc.png
"""
import sys, os, warnings
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v404m_panels import PANELS, STRESS_PANELS, SEX_GENES, load_counts, load_meta, cpm_log  # noqa: E402

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)
MIN_COUNT, MIN_SAMPLES = 10, 3
TISSUE_ORDER = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']

C, source = load_counts(TABLES)
meta = load_meta(TABLES)
common = [s for s in C.columns if s in meta.index]
C, meta = C[common].round().astype(int), meta.loc[common]
print(f'counts from {source}: {C.shape[0]} genes x {C.shape[1]} samples')
print(meta.groupby(['tissue', 'genotype']).size().unstack(fill_value=0).to_string())

# ---------------------------------------------------------------- QC
L = cpm_log(C)
qc = pd.DataFrame({'lib_size_M': C.sum() / 1e6}).join(meta)
for g in SEX_GENES:
    if g in L.index:
        qc[g] = L.loc[g].round(2)
qc['sex_called'] = np.where(qc.get('Xist', pd.Series(0, index=qc.index)) > qc.get('Ddx3y', pd.Series(0, index=qc.index)), 'female', 'male')
qc['sex_mismatch'] = qc.sex_called.str.lower() != qc.sex.str.lower()
qc.to_csv(f'{OUT}/tables/v404m_qc.csv')
if qc.sex_mismatch.any():
    print('\nWARNING: sex mismatch between the sample sheet and Xist/Ddx3y expression:')
    print(qc[qc.sex_mismatch][['tissue', 'genotype', 'sex', 'sex_called'] + [g for g in SEX_GENES if g in qc]].to_string())
    # Xist and the Y genes are unambiguous, so the recorded sex is corrected for the model. Sex is only a
    # covariate here, but a wrong label adds noise to every gene; USE_SHEET_SEX=1 keeps the sheet's values.
    if os.environ.get('USE_SHEET_SEX') != '1':
        meta = meta.copy()
        meta.loc[qc.index[qc.sex_mismatch], 'sex'] = qc.loc[qc.sex_mismatch, 'sex_called']
        print(f'  -> using the RNA-derived sex for these {int(qc.sex_mismatch.sum())} samples '
              f'(set USE_SHEET_SEX=1 to keep the sheet values)')
        print('  group composition after correction:')
        print(meta.groupby(['tissue', 'genotype', 'sex']).size().unstack(fill_value=0).to_string())
else:
    print('\nsex check: all samples match the sample sheet')
print('library sizes: median %.1fM, range %.1f-%.1f' % (qc.lib_size_M.median(), qc.lib_size_M.min(), qc.lib_size_M.max()))


def de_tissue(tissue):
    """DESeq2 (pydeseq2) within one tissue; falls back to a moderated t-test on log-CPM if pydeseq2 is absent."""
    idx = meta.index[meta.tissue == tissue]
    sub = C[idx]
    keep = (sub >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
    sub = sub.loc[keep]
    md = meta.loc[idx].copy()
    both_sexes = md.groupby(['genotype', 'sex']).size().unstack(fill_value=0)
    use_sex = (both_sexes > 0).all().all() if both_sexes.shape[1] > 1 else False
    design = '~sex + genotype' if use_sex else '~genotype'
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        try:
            dds = DeseqDataSet(counts=sub.T, metadata=md, design=design, refit_cooks=True, quiet=True)
        except TypeError:
            dds = DeseqDataSet(counts=sub.T, metadata=md,
                               design_factors=['sex', 'genotype'] if use_sex else ['genotype'],
                               refit_cooks=True, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=['genotype', 'mutant', 'wildtype'], quiet=True)
        st.summary()
        res = st.results_df.copy()
        res['stat'] = res['stat'] if 'stat' in res else res['log2FoldChange'] / res['lfcSE']
        method = f'pydeseq2 {design}'
    except Exception as e:                                    # pragma: no cover
        print(f'  pydeseq2 unavailable or failed ({e}); using log-CPM t-test', file=sys.stderr)
        lg = cpm_log(sub)
        a = lg[md.index[md.genotype == 'mutant']]; b = lg[md.index[md.genotype == 'wildtype']]
        t, p = stats.ttest_ind(a, b, axis=1)
        res = pd.DataFrame({'log2FoldChange': a.mean(axis=1) - b.mean(axis=1), 'stat': t, 'pvalue': p}, index=lg.index)
        res['padj'] = stats.false_discovery_control(res.pvalue.fillna(1))
        method = 'log-CPM t-test ~genotype'
    res['tissue'] = tissue
    res.to_csv(f'{OUT}/tables/v404m_de_{tissue}.csv')
    n_sig = int((res.padj < 0.05).sum())
    print(f'\n{tissue}: {sub.shape[0]} genes tested ({method}), {n_sig} with FDR < 0.05')
    top = res.dropna(subset=['padj']).nsmallest(12, 'padj')[['log2FoldChange', 'pvalue', 'padj']]
    if len(top):
        print(top.round(4).to_string())
    return res, sub, md


def panel_tests(res, sub, md, tissue):
    """Competitive (Wald-stat shift) and self-contained (per-sample z-score) tests for every panel."""
    lg = cpm_log(sub)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0)
    mut = md.index[md.genotype == 'mutant']; wt = md.index[md.genotype == 'wildtype']
    rows, scores = [], []
    for name, genes in PANELS.items():
        present = [g for g in genes if g in res.index]
        if len(present) < 4:
            continue
        s_in = res.loc[present, 'stat'].dropna()
        s_out = res.loc[res.index.difference(present), 'stat'].dropna()
        if len(s_in) < 4 or len(s_out) < 50:
            continue
        u_p = stats.mannwhitneyu(s_in, s_out).pvalue
        zz = z.loc[[g for g in present if g in z.index]].dropna(how='all')
        per_sample = zz.mean(axis=0)
        sc_p = stats.mannwhitneyu(per_sample[mut], per_sample[wt]).pvalue if len(mut) >= 3 and len(wt) >= 3 else np.nan
        rows.append((tissue, name, len(s_in), res.loc[present, 'log2FoldChange'].median(),
                     s_in.median(), s_out.median(), u_p,
                     per_sample[mut].mean() - per_sample[wt].mean(), sc_p,
                     int((res.loc[present, 'padj'] < 0.05).sum())))
        for s, v in per_sample.items():
            scores.append((tissue, name, s, md.at[s, 'genotype'], md.at[s, 'sex'], v))
    P = pd.DataFrame(rows, columns=['tissue', 'panel', 'n_genes', 'median_log2FC', 'median_stat_panel',
                                    'median_stat_rest', 'p_competitive', 'z_diff_mut_minus_wt',
                                    'p_self_contained', 'n_genes_FDR05'])
    for col, out in [('p_competitive', 'fdr_competitive'), ('p_self_contained', 'fdr_self_contained')]:
        ok = P[col].notna()
        P.loc[ok, out] = stats.false_discovery_control(P.loc[ok, col])
    return P, pd.DataFrame(scores, columns=['tissue', 'panel', 'sample', 'genotype', 'sex', 'z'])


tissues = [t for t in TISSUE_ORDER if t in set(meta.tissue)] + [t for t in sorted(set(meta.tissue)) if t not in TISSUE_ORDER]
all_res, all_panels, all_scores = {}, [], []
for t in tissues:
    md_n = meta[meta.tissue == t].groupby('genotype').size()
    if md_n.min() < 2 or len(md_n) < 2:
        print(f'\n{t}: skipped (genotype counts {md_n.to_dict()})'); continue
    res, sub, md = de_tissue(t)
    all_res[t] = res
    P, S = panel_tests(res, sub, md, t)
    all_panels.append(P); all_scores.append(S)

PT = pd.concat(all_panels, ignore_index=True) if all_panels else pd.DataFrame()
SC = pd.concat(all_scores, ignore_index=True) if all_scores else pd.DataFrame()
PT.to_csv(f'{OUT}/tables/v404m_panel_tests.csv', index=False)
SC.to_csv(f'{OUT}/tables/v404m_panel_scores.csv', index=False)

# ---------------------------------------------------------------- headline
print('\n' + '=' * 110)
print('PANEL SUMMARY — V404M vs wild type (positive = higher in mutant)')
print('=' * 110)
show = PT[PT.panel.isin(STRESS_PANELS + ['Immediate-early / activity', 'Kv4 complex', 'Neuronal integrity'])]
print(show.pivot_table(index='panel', columns='tissue', values='median_log2FC')
      .reindex([p for p in PANELS if p in set(show.panel)])[[t for t in tissues if t in set(show.tissue)]]
      .round(3).to_string())
print('\ncompetitive FDR (panel vs rest of transcriptome):')
print(show.pivot_table(index='panel', columns='tissue', values='fdr_competitive')
      .reindex([p for p in PANELS if p in set(show.panel)])[[t for t in tissues if t in set(show.tissue)]]
      .round(4).to_string())
hits = PT[(PT.fdr_competitive < 0.05) & (PT.fdr_self_contained < 0.05)]
print('\npanels significant by BOTH tests:')
print(hits[['tissue', 'panel', 'n_genes', 'median_log2FC', 'fdr_competitive', 'z_diff_mut_minus_wt',
            'fdr_self_contained', 'n_genes_FDR05']].round(4).to_string(index=False) if len(hits) else '  (none)')

# ---------------------------------------------------------------- figures
if len(PT):
    panels_plot = [p for p in PANELS if p in set(PT.panel)]
    tiss = [t for t in tissues if t in set(PT.tissue)]
    fig, ax = plt.subplots(figsize=(1.6 * len(tiss) + 5, 0.42 * len(panels_plot) + 2))
    for i, p in enumerate(panels_plot):
        for j, t in enumerate(tiss):
            r = PT[(PT.panel == p) & (PT.tissue == t)]
            if not len(r):
                continue
            lfc = float(r.median_log2FC.iloc[0]); fdr = float(r.fdr_competitive.iloc[0])
            ax.scatter(j, len(panels_plot) - 1 - i, s=30 + 220 * min(-np.log10(max(fdr, 1e-10)) / 4, 1.6),
                       c=[lfc], cmap='RdBu_r', vmin=-0.4, vmax=0.4, edgecolors='k' if fdr < 0.05 else 'none', linewidths=1.0)
    ax.set_xticks(range(len(tiss))); ax.set_xticklabels(tiss, rotation=20)
    ax.set_yticks(range(len(panels_plot))); ax.set_yticklabels(panels_plot[::-1], fontsize=8)
    sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=plt.Normalize(-0.4, 0.4)); sm.set_array([])
    plt.colorbar(sm, ax=ax, label='median log2FC (mutant / WT)', shrink=0.7)
    ax.set_title('V404M vs wild type: pathway panels\ncircle size = competitive significance, black edge = FDR < 0.05', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/figures/v404m_panels.png', dpi=200); plt.close(fig)

    key = ['ISR (ATF4 output)', 'UPR chaperones (ATF6)', 'Apoptosis (intrinsic)', 'Astrogliosis',
           'Microglia / neuroinflammation', 'Immediate-early / activity', 'Kv4 complex']
    genes = [g for p in key for g in PANELS[p]]
    mat = pd.DataFrame({t: all_res[t].reindex(genes).log2FoldChange for t in tissues if t in all_res}).dropna(how='all')
    if len(mat):
        fig, ax = plt.subplots(figsize=(1.3 * mat.shape[1] + 3, 0.17 * len(mat) + 2))
        im = ax.imshow(mat.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)
        ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=20, fontsize=8)
        ax.set_yticks(range(len(mat))); ax.set_yticklabels(mat.index, fontsize=5)
        plt.colorbar(im, ax=ax, label='log2FC (mutant / WT)', shrink=0.6)
        ax.set_title('Stress, death, gliosis, activity and Kv4 genes', fontsize=10)
        fig.tight_layout(); fig.savefig(f'{OUT}/figures/v404m_stress_genes.png', dpi=200); plt.close(fig)

    kp = [p for p in ['ISR (ATF4 output)', 'UPR chaperones (ATF6)', 'Apoptosis (intrinsic)', 'Astrogliosis',
                      'Microglia / neuroinflammation', 'Immediate-early / activity'] if p in set(SC.panel)]
    fig, axes = plt.subplots(1, len(kp), figsize=(2.5 * len(kp) + 1, 3.6), sharey=False)
    axes = np.atleast_1d(axes)
    for ax, p in zip(axes, kp):
        d = SC[SC.panel == p]
        for j, t in enumerate(tiss):
            for gt, dx, col in [('wildtype', -0.17, '#4C72B0'), ('mutant', 0.17, '#C44E52')]:
                v = d[(d.tissue == t) & (d.genotype == gt)].z
                ax.scatter(np.full(len(v), j + dx) + np.random.uniform(-0.05, 0.05, len(v)), v, s=12, color=col,
                           label=gt if j == 0 else None)
                if len(v): ax.plot([j + dx - 0.1, j + dx + 0.1], [v.mean()] * 2, color='k', lw=1.2)
        ax.set_xticks(range(len(tiss))); ax.set_xticklabels([t[:4] for t in tiss], fontsize=8)
        ax.set_title(p, fontsize=8); ax.axhline(0, color='#999', lw=0.6)
    axes[0].set_ylabel('per-sample panel z-score'); axes[0].legend(fontsize=7, frameon=False)
    fig.suptitle('Per-sample pathway scores by genotype (each point one animal)', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{OUT}/figures/v404m_qc.png', dpi=200); plt.close(fig)
print('\ndone')
