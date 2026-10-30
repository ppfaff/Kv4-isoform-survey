#!/usr/bin/env python3
"""Read-level anatomy of Kcnd2 in Kv4-null vs Kv4-low vs Kv4-high DG granule cells (Allen SSv4).

Inputs (from pipeline/run_dg_star.sh): tables/dg_readlevel_coverage.tsv, tables/dg_readlevel_junctions.tsv,
tables/dg_readlevel_runlist.tsv. Three questions, each a distinct read-level signature:
  A. Where along the gene are the intronic reads?  Nascent transcription caught mid-gene -> 5'-biased; complete but
     unspliced pre-mRNA -> flat; a detained intron -> one intron carries most of the signal.
  B. Is there residual spliced mRNA in the null cells?  Junction reads per cell; exon-block depth 5'->3' (CDS vs 3'UTR).
  C. Are Kcnip4 and Dpp6 read profiles normal in the same cells? (gene-specific vs general processing defect)
Outputs: tables/dg_readlevel_summary.csv, tables/dg_readlevel_per_cell.csv, figures/dg_readlevel.png
Usage:   python3 scripts/dg_readlevel.py [TABLES_DIR] [OUT_DIR]
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)
GRPS = ['Kv4-null', 'Kv4-low', 'Kv4-high']; COL = {'Kv4-null': '#D6372E', 'Kv4-low': '#E8822D', 'Kv4-high': '#2E6FBA'}

C = pd.read_csv(f'{TABLES}/dg_readlevel_coverage.tsv', sep='\t')
C['gene'] = C.name.str.replace(r'_(intron|exon)_\d+$', '', regex=True)
C['k'] = C.name.str.extract(r'_(\d+)$').astype(int)
C['len'] = C.end - C.start
C['depth'] = C.depth_sum / C.len
strand = C.groupby('gene').strand.first()
# order features 5'->3' in transcript orientation: introns were numbered 5'->3' already; exon blocks were numbered by
# genomic position, so reverse them for minus-strand genes
C['k5'] = np.where((C.feature == 'exon') & (C.gene.map(strand) == '-'),
                   C.groupby(['gene', 'feature']).k.transform('max') - C.k + 1, C.k)
R = pd.read_csv(f'{TABLES}/dg_readlevel_runlist.tsv', sep='\t')
ncell = R.grp.value_counts()
print('cells per group:', dict(ncell))

# ---------------------------------------------------------------- A/B: per-cell profiles for the three genes
rows = []
for (s, g), d in C.groupby(['sample_name', 'gene']):
    grp = d.grp.iloc[0]
    intr = d[d.feature == 'intron'].sort_values('k5'); ex = d[d.feature == 'exon'].sort_values('k5')
    it = intr.depth_sum.sum(); et = ex.depth_sum.sum()
    # position of intronic signal: depth-weighted mean of relative position (0 = 5', 1 = 3'); detained-intron score = max share
    pos = (intr.depth * (intr.k5 - 1) / max(len(intr) - 1, 1)).sum() / intr.depth.sum() if intr.depth.sum() > 0 else np.nan
    top = intr.depth_sum.max() / it if it > 0 else np.nan
    n_ex = len(ex); half = max(n_ex // 2, 1)
    cds5 = ex.iloc[:half].depth_sum.sum(); utr3 = ex.iloc[half:].depth_sum.sum()
    rows.append((s, grp, g, it, et, et / (et + it) if et + it > 0 else np.nan, pos, top, intr.depth_sum.idxmax() if it > 0 else None,
                 intr.loc[intr.depth_sum.idxmax(), 'k5'] if it > 0 else np.nan, cds5, utr3))
P = pd.DataFrame(rows, columns=['sample_name', 'grp', 'gene', 'intron_bases', 'exon_bases', 'exon_frac', 'intron_pos_5to3',
                                'top_intron_share', '_i', 'top_intron_k', 'exon_5half_bases', 'exon_3half_bases']).drop(columns='_i')
J = None
jp = f'{TABLES}/dg_readlevel_junctions.tsv'
if os.path.exists(jp):
    J = pd.read_csv(jp, sep='\t')
    genes_bed = C[['gene', 'chr']].drop_duplicates()
    # assign each junction to a gene by its intron intervals
    iv = C[C.feature == 'intron'][['gene', 'chr', 'start', 'end']].drop_duplicates()
    def gene_of(r):
        m = iv[(iv.chr == r.chr) & (iv.start <= r.start) & (iv.end >= r.end)]
        return m.gene.iloc[0] if len(m) else None
    J['gene'] = J.apply(gene_of, axis=1)
    jc = J.dropna(subset=['gene']).groupby(['sample_name', 'gene']).uniq.sum().rename('junction_reads').reset_index()
    P = P.merge(jc, on=['sample_name', 'gene'], how='left'); P['junction_reads'] = P.junction_reads.fillna(0)
P.to_csv(f'{OUT}/tables/dg_readlevel_per_cell.csv', index=False)

summ = P.groupby(['gene', 'grp']).agg(n=('sample_name', 'nunique'),
                                      intron_bases_med=('intron_bases', 'median'), exon_bases_med=('exon_bases', 'median'),
                                      exon_frac_med=('exon_frac', 'median'),
                                      intron_pos_mean=('intron_pos_5to3', 'mean'), top_intron_share_med=('top_intron_share', 'median'),
                                      pct_cells_detained=('top_intron_share', lambda x: (x > 0.8).mean() * 100),
                                      junction_reads_med=('junction_reads', 'median') if 'junction_reads' in P else ('exon_bases', 'size'),
                                      pct_cells_any_junction=('junction_reads', lambda x: (x > 0).mean() * 100) if 'junction_reads' in P else ('exon_bases', 'size'),
                                      utr3_over_cds5=('exon_3half_bases', lambda x: x.sum()) ).reset_index()
summ['utr3_over_cds5'] = P.groupby(['gene', 'grp']).apply(lambda d: d.exon_3half_bases.sum() / max(d.exon_5half_bases.sum(), 1)).values
summ.to_csv(f'{OUT}/tables/dg_readlevel_summary.csv', index=False)
pd.set_option('display.width', 220)
print(summ.round(3).to_string(index=False))

# group-mean intron profiles for Kcnd2 (each cell normalised to its own intronic total)
prof = {}
for g in ['Kcnd2', 'Kcnip4', 'Dpp6']:
    d = C[(C.gene == g) & (C.feature == 'intron')].copy()
    d['tot'] = d.groupby('sample_name').depth_sum.transform('sum'); d = d[d.tot > 0]
    d['share'] = d.depth / d.groupby('sample_name').depth.transform('sum')
    prof[g] = d.groupby(['grp', 'k5']).share.agg(['mean', 'sem']).reset_index()
    print(f'\n{g} intron share by position (mean over cells with intronic reads):')
    print(prof[g].pivot(index='grp', columns='k5', values='mean').round(3).to_string())

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 4, figsize=(17, 4))
ax = axes[0]
for grp in GRPS:
    p = prof['Kcnd2'][prof['Kcnd2'].grp == grp]
    ax.errorbar(p.k5, p['mean'], yerr=p['sem'], color=COL[grp], marker='o', ms=3, label=grp)
ax.set_xlabel("Kcnd2 intron (5' → 3')"); ax.set_ylabel('share of cell\'s intronic depth'); ax.set_title('A. Where the intronic reads are'); ax.legend(fontsize=8, frameon=False)
ax = axes[1]
for i, grp in enumerate(GRPS):
    v = P[(P.gene == 'Kcnd2') & (P.grp == grp)].top_intron_share.dropna()
    ax.scatter(np.full(len(v), i) + np.random.uniform(-0.2, 0.2, len(v)), v, s=6, color=COL[grp], alpha=0.6)
ax.axhline(0.8, color='k', lw=0.6, ls='--'); ax.set_xticks(range(3)); ax.set_xticklabels(GRPS); ax.set_ylabel('largest single-intron share'); ax.set_title('B. Detained-intron test')
ax = axes[2]
if 'junction_reads' in P:
    for i, grp in enumerate(GRPS):
        v = P[(P.gene == 'Kcnd2') & (P.grp == grp)].junction_reads
        ax.scatter(np.full(len(v), i) + np.random.uniform(-0.2, 0.2, len(v)), v + 0.5, s=6, color=COL[grp], alpha=0.6)
    ax.set_yscale('log'); ax.set_xticks(range(3)); ax.set_xticklabels(GRPS); ax.set_ylabel('Kcnd2 splice-junction reads (+0.5)'); ax.set_title('C. Residual spliced mRNA')
ax = axes[3]
for g, m in [('Kcnd2', 'o'), ('Kcnip4', 's'), ('Dpp6', '^')]:
    v = P[P.gene == g].groupby('grp').exon_frac.median().reindex(GRPS)
    ax.plot(range(3), v, marker=m, label=g)
ax.set_xticks(range(3)); ax.set_xticklabels(GRPS); ax.set_ylabel('median exonic fraction of gene reads'); ax.set_title('D. Gene-specific or general?'); ax.legend(fontsize=8, frameon=False)
fig.suptitle(f'Read-level Kcnd2 in DG granule cells (Allen SSv4): null n={ncell.get("Kv4-null", 0)}, low n={ncell.get("Kv4-low", 0)}, high n={ncell.get("Kv4-high", 0)}', fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/figures/dg_readlevel.png', dpi=200)
print('\ndone')
