#!/usr/bin/env python3
"""Dendritic vs somatic expression of the Kv4-complex subunits and their variants.

Inputs (from pipeline/collapse_dendrite.py, and optionally pipeline/run_dendrite_star.sh), read from TABLES:
  dendrite_isoform_per_sample.csv, dendrite_gene_tpm.csv.gz,
  star_exon_intron_counts.tsv + star_intron_depth.tsv + star_ReadsPerGene.tsv  (Farris only; skipped if absent)

Questions answered
  1. QC: do known dendritic / somatic markers separate the compartments as expected? (per dataset)
  2. Gene level: log2(dendrite/soma) TPM for the 9 genes — Farris per subfield (paired by animal),
     Middleton per neuron (paired) — with a paired test and a global-normalisation alternative.
  3. Variant level: fraction of each gene's reads per variant, by compartment; paired dendrite-soma
     differences (Kv4.2 : Kv4.3S : Kv4.3L pool composition, DPP6a/S/L/K, KChIP N-terminal forms).
  4. Exonic vs intronic reads (Farris, STAR): exon:intron ratio per gene per compartment, and per-intron
     Kcnd2 depth relative to exon depth in cell-body vs dendritic layers.

Outputs: tables/dendrite_*.csv and figures/dendrite_*.png next to the repo (or OUT=dir).
Usage:  python3 scripts/dendrite_compartment.py [TABLES_DIR] [OUT_DIR]
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import VSPEC, HATCH, LEGEND_GROUPS, BLANK, PANEL_OF  # noqa: E402

try:
    from scipy import stats
except ImportError:  # pragma: no cover
    stats = None

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)

GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']
PANELS = ['Kv4', 'KChIP', 'DPP']
# literature dendritic / somatic-nuclear marker sets used only as a QC on compartment separation
DEND_MARKERS = ['Camk2a', 'Map2', 'Dlg4', 'Shank1', 'Ddn', 'Arc', 'Bdnf', 'Map1b', 'Calm1', 'Actb']
SOMA_MARKERS = ['Malat1', 'Meg3', 'Snhg11', 'Neat1', 'Xist', 'Kcnq1ot1', 'Rpl13a', 'Rps27a', 'Hist1h1c', 'Gapdh']
MIN_GENE_READS = 20     # a sample's variant fractions for a gene are used only above this
PSEUDO = 0.1            # TPM pseudocount for log ratios

D = pd.read_csv(f'{TABLES}/dendrite_isoform_per_sample.csv')
G = pd.read_csv(f'{TABLES}/dendrite_gene_tpm.csv.gz', index_col=0)
meta = D[['dataset', 'sample_id', 'region', 'compartment', 'unit']].drop_duplicates().set_index('sample_id')
print(f'{len(meta)} samples:'); print(meta.groupby(['dataset', 'region', 'compartment']).size().to_string())


def paired(df, value):
    """Pivot a sample-level column into (dataset, region, unit) rows with dendrite / soma columns."""
    p = df.pivot_table(index=['dataset', 'region', 'unit'], columns='compartment', values=value)
    return p.dropna(subset=['dendrite', 'soma'])


def ptest(a, b):
    """Paired test: Wilcoxon when n>=6 (otherwise paired t), returns (stat_name, p)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if stats is None or len(a) < 3:
        return ('none', np.nan)
    if len(a) >= 6:
        try:
            return ('wilcoxon', stats.wilcoxon(a, b).pvalue)
        except ValueError:
            return ('wilcoxon', np.nan)
    return ('paired_t', stats.ttest_rel(a, b).pvalue)


# ---------------------------------------------------------------- 1. marker QC
rows = []
for g in DEND_MARKERS + SOMA_MARKERS:
    if g not in G.index:
        continue
    s = G.loc[g].rename('TPM').to_frame().join(meta)
    p = paired(s.reset_index(), 'TPM')
    lr = np.log2((p.dendrite + PSEUDO) / (p.soma + PSEUDO))
    for ds, grp in lr.groupby(level='dataset'):
        rows.append((ds, g, 'dendritic' if g in DEND_MARKERS else 'somatic/nuclear', grp.mean(), grp.median(), len(grp)))
M = pd.DataFrame(rows, columns=['dataset', 'gene', 'expected', 'mean_log2_dend_soma', 'median_log2', 'n_pairs'])
M.to_csv(f'{OUT}/tables/dendrite_marker_qc.csv', index=False)
print('\nmarker QC (mean log2 dendrite/soma):'); print(M.pivot_table(index=['expected', 'gene'], columns='dataset', values='mean_log2_dend_soma').round(2).to_string())

# ---------------------------------------------------------------- 2. gene level
gene_tpm = D.groupby(['dataset', 'sample_id', 'region', 'compartment', 'unit', 'gene'], as_index=False).agg(TPM=('TPM', 'sum'), reads=('reads', 'sum'))
# alternative normalisation: TPM relative to the median of the dendritic-marker set, to guard against the
# compartment's bulk composition (e.g. glial/axonal RNA in neuropil) shifting all neuronal TPMs together
ref = G.loc[[g for g in DEND_MARKERS if g in G.index]].median()
gene_tpm['TPM_rel_dend_markers'] = gene_tpm.TPM / gene_tpm.sample_id.map(ref)
res = []
for val in ['TPM', 'TPM_rel_dend_markers']:
    for g in GENES:
        p = paired(gene_tpm[gene_tpm.gene == g], val)
        for (ds, rg), grp in p.groupby(level=['dataset', 'region']):
            lr = np.log2((grp.dendrite + PSEUDO) / (grp.soma + PSEUDO))
            name, pv = ptest(np.log2(grp.dendrite + PSEUDO), np.log2(grp.soma + PSEUDO))
            res.append((val, ds, rg, g, PANEL_OF[g], grp.soma.mean(), grp.dendrite.mean(), lr.mean(), lr.std(ddof=1) if len(lr) > 1 else np.nan, len(lr), name, pv))
        # Farris: all subfields pooled (paired by animal x subfield)
        if val == 'TPM' or True:
            pf = p.loc['farris'] if 'farris' in p.index.get_level_values(0) else None
            if pf is not None and len(pf) > 1:
                lr = np.log2((pf.dendrite + PSEUDO) / (pf.soma + PSEUDO))
                name, pv = ptest(np.log2(pf.dendrite + PSEUDO), np.log2(pf.soma + PSEUDO))
                res.append((val, 'farris', 'all', g, PANEL_OF[g], pf.soma.mean(), pf.dendrite.mean(), lr.mean(), lr.std(ddof=1), len(lr), name, pv))
R = pd.DataFrame(res, columns=['normalisation', 'dataset', 'region', 'gene', 'panel', 'mean_soma', 'mean_dendrite', 'mean_log2_dend_soma', 'sd_log2', 'n_pairs', 'test', 'p'])
R.to_csv(f'{OUT}/tables/dendrite_gene_level.csv', index=False)
print('\ngene level, log2(dendrite/soma) TPM:')
print(R[R.normalisation == 'TPM'].pivot_table(index='gene', columns=['dataset', 'region'], values='mean_log2_dend_soma').round(2).reindex(GENES).to_string())

# figure: one panel per dataset/region block, bars = mean log2 ratio with per-pair points
blocks = [(ds, rg) for ds, rg in R[(R.normalisation == 'TPM')][['dataset', 'region']].drop_duplicates().values if rg != 'all']
fig, axes = plt.subplots(1, len(blocks), figsize=(2.6 * len(blocks) + 1, 4.2), sharey=True)
axes = np.atleast_1d(axes)
for ax, (ds, rg) in zip(axes, blocks):
    for i, g in enumerate(GENES):
        p = paired(gene_tpm[(gene_tpm.gene == g) & (gene_tpm.dataset == ds) & (gene_tpm.region == rg)], 'TPM')
        lr = np.log2((p.dendrite + PSEUDO) / (p.soma + PSEUDO))
        col = {'Kv4': '#2E6FBA', 'KChIP': '#4CAF50', 'DPP': '#7B4FA3'}[PANEL_OF[g]]
        ax.bar(i, lr.mean(), color=col, alpha=0.8)
        ax.scatter(np.full(len(lr), i) + np.random.uniform(-0.15, 0.15, len(lr)), lr, s=10, color='k', zorder=3)
    ax.axhline(0, color='k', lw=0.8); ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=90, fontsize=8)
    ax.set_title(f'{ds} {rg}', fontsize=10)
axes[0].set_ylabel('log2 (dendrite / soma) TPM')
fig.suptitle('Kv4-complex genes: dendritic vs somatic compartment (points = animals / neurons)', fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/figures/dendrite_gene_level.png', dpi=200); plt.close(fig)

# ---------------------------------------------------------------- 3. variant composition by compartment
V = D.dropna(subset=['vclass']).copy()
is_other = V.vclass.str.contains('other|fragment', case=False, regex=True)
V = V[~((V.biotype == 'retained_intron') & is_other)]
# pooled-read composition per (dataset, region, compartment) for the figure
agg = V.groupby(['dataset', 'region', 'compartment', 'panel', 'vclass'], as_index=False).reads.sum()
agg['panel_total'] = agg.groupby(['dataset', 'region', 'compartment', 'panel']).reads.transform('sum')
agg['frac'] = agg.reads / agg.panel_total.replace(0, np.nan)
agg.to_csv(f'{OUT}/tables/dendrite_variant_fracs_by_compartment.csv', index=False)
# per-sample variant fraction of gene, then paired difference dendrite - soma
vs = V.groupby(['dataset', 'sample_id', 'region', 'compartment', 'unit', 'gene', 'vclass'], as_index=False).agg(reads=('reads', 'sum'), gene_reads=('gene_reads', 'first'))
vs = vs[vs.gene_reads >= MIN_GENE_READS]
vs['frac'] = vs.reads / vs.gene_reads
vres = []
for (g, v), sub in vs.groupby(['gene', 'vclass']):
    p = paired(sub, 'frac')
    for (ds, rg), grp in p.groupby(level=['dataset', 'region']):
        name, pv = ptest(grp.dendrite, grp.soma)
        vres.append((ds, rg, g, v, grp.soma.mean(), grp.dendrite.mean(), (grp.dendrite - grp.soma).mean(), len(grp), name, pv))
    if 'farris' in p.index.get_level_values(0):
        pf = p.loc['farris']
        if len(pf) > 1:
            name, pv = ptest(pf.dendrite, pf.soma)
            vres.append(('farris', 'all', g, v, pf.soma.mean(), pf.dendrite.mean(), (pf.dendrite - pf.soma).mean(), len(pf), name, pv))
VR = pd.DataFrame(vres, columns=['dataset', 'region', 'gene', 'vclass', 'mean_frac_soma', 'mean_frac_dendrite', 'mean_diff', 'n_pairs', 'test', 'p'])
VR.to_csv(f'{OUT}/tables/dendrite_variant_frac_paired.csv', index=False)
print('\nvariant fraction of gene, dendrite - soma (mean paired difference), |diff| > 0.05:')
print(VR[VR.mean_diff.abs() > 0.05].sort_values(['gene', 'dataset', 'region']).round(3).to_string(index=False))

# figure: stacked bars, rows = region x compartment, columns = panels (same scheme as the per-type figures)
rows_ = [(ds, rg, c) for ds in ['farris', 'middleton'] for rg in (['CA1', 'CA2', 'CA3', 'DG'] if ds == 'farris' else ['culture']) for c in ['soma', 'dendrite']]
rows_ = [r for r in rows_ if ((agg.dataset == r[0]) & (agg.region == r[1]) & (agg.compartment == r[2])).any()]
ypos = {}
y = 0.0
for i, r in enumerate(rows_):
    if i > 0 and (r[0], r[1]) != (rows_[i - 1][0], rows_[i - 1][1]):
        y += 0.8
    ypos[r] = y; y += 1.0
ymax = max(ypos.values())
fig, axes = plt.subplots(1, 3, figsize=(16, 0.42 * (ymax + 2) + 4.5), sharey=True)
for ax, panel in zip(axes, PANELS):
    sub = agg[agg.panel == panel]
    for r in rows_:
        row = sub[(sub.dataset == r[0]) & (sub.region == r[1]) & (sub.compartment == r[2])]
        yy = ymax - ypos[r]
        if row.reads.sum() < MIN_GENE_READS or row.empty:
            ax.barh(yy, 1.0, color='#EEEEEE', height=0.75); continue
        left = 0.0
        for _, q in row[row.frac > 0].sort_values('frac', ascending=False).iterrows():
            base, kind, pcolor = VSPEC[q.vclass]
            ax.barh(yy, q.frac, left=left, color=base, edgecolor=pcolor, hatch=HATCH[kind], height=0.75, linewidth=0.3)
            left += q.frac
        ax.text(1.01, yy, f'{int(row.reads.sum()):,}', va='center', fontsize=7, color='#555')
    ax.set_xlim(0, 1.15); ax.set_ylim(-0.75, ymax + 0.75); ax.set_xticks([]); ax.set_title(panel, fontsize=13, fontweight='bold')
    for sp in ax.spines.values(): sp.set_visible(False)
    groups = LEGEND_GROUPS[panel]; xs = [0.5] if len(groups) == 1 else [0.28, 0.76]
    legs = []
    for x, vl in zip(xs, groups):
        pats = [mpatches.Patch(facecolor='none', edgecolor='none', label=' ') if v == BLANK else
                mpatches.Patch(facecolor=VSPEC[v][0], edgecolor=VSPEC[v][2], hatch=HATCH[VSPEC[v][1]], label=v, linewidth=0.3) for v in vl]
        legs.append(ax.legend(handles=pats, loc='upper center', bbox_to_anchor=(x, -0.02), fontsize=8, frameon=False, handlelength=3.0, handleheight=2.4))
    for leg in legs[:-1]: ax.add_artist(leg)
axes[0].set_yticks([ymax - ypos[r] for r in rows_])
axes[0].set_yticklabels([f'{r[0][:3]} {r[1]} {r[2]}' for r in rows_], fontsize=9)
fig.suptitle('Variant composition by compartment — pooled reads per dataset / subfield / compartment\n'
             'Farris 2019 laminar LCM (cell-body vs dendritic layers), Middleton 2019 single-neuron soma vs dendrite; number = family reads',
             fontsize=11)
fig.subplots_adjust(left=0.14, right=0.98, top=0.88, bottom=0.42, wspace=0.08)
fig.savefig(f'{OUT}/figures/dendrite_variant_scheme_by_compartment.png', dpi=200); plt.close(fig)

# Kv4 pool composition explicitly (the headline question)
kv4 = agg[agg.panel == 'Kv4'].pivot_table(index=['dataset', 'region', 'compartment'], columns='vclass', values='frac').fillna(0)
kv4.round(3).to_csv(f'{OUT}/tables/dendrite_kv4_pool_composition.csv')
print('\nKv4 pool composition (fraction of Kv4-family reads):'); print(kv4.round(3).to_string())

# ---------------------------------------------------------------- 4. STAR exon / intron (Farris)
ei_path = f'{TABLES}/star_exon_intron_counts.tsv'
if os.path.exists(ei_path):
    E = pd.read_csv(ei_path, sep='\t')
    # pick the strand setting that assigns the most exonic reads overall (TruSeq stranded -> expect s=2)
    s_best = E[E.feature == 'exon'].groupby('strand_setting')['count'].sum().idxmax()
    print(f'\nSTAR/featureCounts: using strand setting -s {s_best}')
    E = E[E.strand_setting == s_best].pivot_table(index=['run', 'gene'], columns='feature', values='count').reset_index()
    run2sample = pd.read_csv(f'{TABLES}/dendrite_samples.tsv', sep='\t', dtype=str)
    if 'extra' in run2sample.columns:
        run2sample = run2sample[run2sample.extra != '1']
    run2sample['sample_id'] = run2sample.dataset + '_' + run2sample.region + '_' + run2sample.compartment + '_' + run2sample.unit
    E = E.merge(run2sample[['run_accession', 'sample_id', 'region', 'compartment', 'unit']], left_on='run', right_on='run_accession')
    E['exon_intron'] = (E.exon + 1) / (E.intron + 1)
    E['intron_frac'] = E.intron / (E.exon + E.intron).replace(0, np.nan)
    E['dataset'] = 'farris'
    E.to_csv(f'{OUT}/tables/dendrite_exon_intron_per_sample.csv', index=False)
    eres = []
    for g in GENES:
        p = paired(E[E.gene == g], 'intron_frac')
        for (ds, rg), grp in p.groupby(level=['dataset', 'region']):
            name, pv = ptest(grp.dendrite, grp.soma)
            eres.append((rg, g, grp.soma.mean(), grp.dendrite.mean(), len(grp), name, pv))
        if len(p):
            name, pv = ptest(p.dendrite, p.soma)
            eres.append(('all', g, p.soma.mean(), p.dendrite.mean(), len(p), name, pv))
    ER = pd.DataFrame(eres, columns=['region', 'gene', 'intron_frac_soma_layer', 'intron_frac_dendrite_layer', 'n_pairs', 'test', 'p'])
    ER.to_csv(f'{OUT}/tables/dendrite_exon_intron_paired.csv', index=False)
    print('\nintronic fraction of gene reads (Farris), cell-body layer vs dendritic layer:')
    print(ER.pivot_table(index='gene', columns='region', values=['intron_frac_soma_layer', 'intron_frac_dendrite_layer']).round(3).reindex(GENES).to_string())

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.8), sharey=True)
    for ax, rg in zip(axes, ['CA1', 'CA2', 'CA3', 'DG']):
        sub = E[E.region == rg]
        for i, g in enumerate(GENES):
            for c, dx, col in [('soma', -0.18, '#C1593C'), ('dendrite', 0.18, '#2E6FBA')]:
                v = sub[(sub.gene == g) & (sub.compartment == c)].intron_frac
                ax.bar(i + dx, v.mean(), width=0.34, color=col, alpha=0.85)
                ax.scatter(np.full(len(v), i + dx), v, s=8, color='k', zorder=3)
        ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=90, fontsize=8); ax.set_title(rg)
    axes[0].set_ylabel('intronic fraction of gene reads')
    axes[0].legend(handles=[mpatches.Patch(color='#C1593C', label='cell-body layer'), mpatches.Patch(color='#2E6FBA', label='dendritic layer')], fontsize=8, frameon=False)
    fig.suptitle('Farris 2019: intronic read fraction per gene, cell-body vs dendritic layer (STAR genome alignment)', fontsize=11)
    fig.tight_layout(); fig.savefig(f'{OUT}/figures/dendrite_exon_intron.png', dpi=200); plt.close(fig)

    # per-intron Kcnd2 depth, normalised to that sample's Kcnd2 exonic read count
    dp = f'{TABLES}/star_intron_depth.tsv'
    if os.path.exists(dp):
        I = pd.read_csv(dp, sep='\t')
        I = I[I.intron.str.startswith('Kcnd2_')].merge(run2sample[['run_accession', 'region', 'compartment', 'unit']], left_on='run', right_on='run_accession')
        ex = E[E.gene == 'Kcnd2'][['run', 'exon']]
        I = I.merge(ex, on='run'); I['depth_per_1k_exon_reads'] = I.mean_depth / (I.exon / 1000).replace(0, np.nan)
        I['k'] = I.intron.str.extract(r'_intron_(\d+)').astype(int)
        I.to_csv(f'{OUT}/tables/dendrite_kcnd2_intron_depth.csv', index=False)
        piv = I.groupby(['region', 'compartment', 'k']).depth_per_1k_exon_reads.mean().unstack('k')
        print('\nKcnd2 per-intron mean depth per 1,000 exonic reads (Farris):'); print(piv.round(3).to_string())
        fig, axes = plt.subplots(1, 4, figsize=(14, 3.6), sharey=True)
        for ax, rg in zip(axes, ['CA1', 'CA2', 'CA3', 'DG']):
            for c, col in [('soma', '#C1593C'), ('dendrite', '#2E6FBA')]:
                sub = I[(I.region == rg) & (I.compartment == c)].groupby('k').depth_per_1k_exon_reads.agg(['mean', 'sem'])
                ax.errorbar(sub.index, sub['mean'], yerr=sub['sem'], color=col, marker='o', ms=3, lw=1, label=f'{c} layer')
            ax.set_title(rg); ax.set_xlabel('Kcnd2 intron (5\' → 3\')')
        axes[0].set_ylabel('mean depth / 1,000 exonic reads'); axes[0].legend(fontsize=8, frameon=False)
        fig.suptitle('Kcnd2 intron coverage by layer — which introns carry the intronic signal, and is it compartment-specific?', fontsize=11)
        fig.tight_layout(); fig.savefig(f'{OUT}/figures/dendrite_kcnd2_intron_depth.png', dpi=200); plt.close(fig)
else:
    print('\n(no STAR tables found — exon/intron section skipped)')
print('\ndone')
