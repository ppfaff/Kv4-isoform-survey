#!/usr/bin/env python3
"""Kv4-specific readouts of the V404M bulk RNA-seq.

  1. allele-specific Kcnd2 expression at codon 404 — zygosity and whether the mutant mRNA is under-represented
  2. Kv4-complex gene levels and auxiliary-subunit variant composition, mutant vs wild type, per tissue
  3. Kcnd2 exonic:intronic ratio — a population readout of the Kv4-null fraction seen in single granule cells
  4. junction readouts: XBP1 splicing (the IRE1 arm of the UPR), the dFosB splice, Kv4.3L and the Kcnip3 acceptor

Allele counting uses the anchor counts taken straight from the FASTQs (tables/v404m_allele_counts.tsv) in
preference to the alignment pileup, because a read carrying the mutation has one mismatch and aligns marginally
less well, which biases alignment-based allele ratios against the mutant; anchor matching has no such bias.

Usage:  python3 scripts/v404m_kv4.py [TABLES_DIR] [OUT_DIR]
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from v404m_panels import load_counts, load_meta, cpm_log  # noqa: E402
try:
    from hip_common import VSPEC, HATCH  # noqa: F401
except Exception:
    VSPEC = HATCH = None

TABLES = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(f'{OUT}/tables', exist_ok=True); os.makedirs(f'{OUT}/figures', exist_ok=True)
GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']
TISSUE_ORDER = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']
CODON_TABLE = {}
for i, b1 in enumerate('TCAG'):
    for j, b2 in enumerate('TCAG'):
        for k, b3 in enumerate('TCAG'):
            CODON_TABLE[b1 + b2 + b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[i * 16 + j * 4 + k]

meta = load_meta(TABLES)
tissues = [t for t in TISSUE_ORDER if t in set(meta.tissue)]
C, source = load_counts(TABLES)
common = [s for s in C.columns if s in meta.index]
C, m2 = C[common], meta.loc[common]
kcnd2_cpm = (C.loc['Kcnd2'] / C.sum() * 1e6) if 'Kcnd2' in C.index else None


def pp(df, title):
    print(f'\n--- {title}')
    if isinstance(df, pd.DataFrame):
        print(df.to_string(index=not isinstance(df.index, pd.RangeIndex)))
    else:
        print(df.to_string())


def mw(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or len(b) < 3:
        return np.nan
    return stats.mannwhitneyu(a, b).pvalue


# ---------------------------------------------------------------- 1. allele-specific expression
allele_path = f'{TABLES}/v404m_allele_counts.tsv'
codon_path = f'{TABLES}/v404m_codon.tsv'
if os.path.exists(allele_path) and os.path.exists(codon_path):
    cfg = dict(l.rstrip('\n').split('\t', 1) for l in open(codon_path) if '\t' in l)
    wt_codon = cfg['wt_codon']
    A = pd.read_csv(allele_path, sep='\t').merge(meta, left_on='sample', right_index=True, how='left')
    A['aa'] = A.triplet.map(lambda t: CODON_TABLE.get(t, 'X'))
    A['is_wt'] = A.triplet == wt_codon
    A['wt_count'] = A['count'].where(A.is_wt, 0)
    tot = A.groupby(['sample', 'tissue', 'genotype', 'anchor'], as_index=False).agg(
        reads=('count', 'sum'), wt=('wt_count', 'sum'))
    tot['alt'] = tot.reads - tot.wt
    tot['alt_frac'] = tot.alt / tot.reads.replace(0, np.nan)
    tot.to_csv(f'{OUT}/tables/v404m_allele_per_sample.csv', index=False)
    pp(A.groupby(['genotype', 'triplet', 'aa'], as_index=False)['count'].sum().sort_values('count', ascending=False).head(12),
       f'triplets observed at codon {cfg.get("residue", "404")} (wild-type codon {wt_codon} = {CODON_TABLE.get(wt_codon, "?")})')
    pp(tot.groupby(['genotype', 'anchor']).agg(n=('sample', 'nunique'), median_reads=('reads', 'median'),
                                               median_alt_frac=('alt_frac', 'median')).round(3),
       'non-wild-type fraction of Kcnd2 reads at the codon (0.5 = heterozygous and balanced)')
    pp(tot[tot.genotype == 'mutant'].groupby(['tissue', 'anchor']).alt_frac.agg(['median', 'min', 'max', 'size']).round(3),
       'mutant animals, by tissue and anchor')
    if 'haplotype' in set(tot.anchor):
        h = tot[tot.anchor == 'haplotype']
        print('\nPRIMARY READOUT — whole-allele haplotype counts (immune to knock-in changes inside an anchor):')
        print(h.groupby(['tissue', 'genotype']).alt_frac.agg(['median', 'min', 'max', 'size']).round(3).to_string())
        hm = h[h.genotype == 'mutant'].alt_frac.median()
        print(f'\nmutant animals: {hm:.3f} of Kcnd2 reads carry the knock-in allele -> '
              + ('heterozygous, both alleles expressed about equally' if 0.4 <= hm <= 0.6 else
                 'homozygous' if hm > 0.85 else
                 'the knock-in allele is UNDER-represented relative to wild type — worth checking before '
                 'interpreting any physiology'))
    agree = tot.pivot_table(index=['sample', 'genotype'], columns='anchor', values='alt_frac')
    if {'left', 'right'} <= set(agree.columns):
        d = (agree.left - agree.right).abs()
        print(f'\nleft vs right anchor agreement: median |difference| {d.median():.3f}, max {d.max():.3f} '
              f'({"anchors agree — neither overlaps an engineered change" if d.max() < 0.1 else "ANCHORS DISAGREE — one overlaps an edit in the knock-in allele; trust the one that still recovers both alleles"})')

    # The knock-in carries at least one silent change besides V404M (a diagnostic restriction site). If it falls
    # inside an anchor, that anchor stops matching the mutant allele. Two signatures reveal it: the captured
    # windows (the mutant haplotype is simply a second window sequence), and a drop in anchored reads in mutants
    # relative to their Kcnd2 expression — if BOTH anchors were blind to the mutant allele, heterozygotes would
    # show ~0% alt while recovering only ~half as many anchored reads as wild types.
    if 'window' in A.columns:
        W = A.groupby(['genotype', 'anchor', 'window'], as_index=False)['count'].sum()
        W['frac_within'] = W['count'] / W.groupby(['genotype', 'anchor'])['count'].transform('sum')
        pp(W[W.frac_within > 0.005].sort_values(['anchor', 'genotype', 'count'], ascending=[True, True, False]).round(4),
           'sequence windows captured around the codon (left anchor = codon + 12 nt downstream; '
           'right anchor = 12 nt upstream + codon). A second window in mutants only is the knock-in haplotype.')
    if kcnd2_cpm is not None:
        depth = tot.groupby(['sample', 'genotype', 'tissue'], as_index=False).reads.sum()
        depth['kcnd2_cpm'] = depth['sample'].map(kcnd2_cpm)
        depth['anchored_per_cpm'] = depth.reads / depth.kcnd2_cpm.replace(0, np.nan)
        s = depth.groupby('genotype').anchored_per_cpm.median()
        pp(depth.groupby(['tissue', 'genotype']).anchored_per_cpm.agg(['median', 'size']).round(2),
           'anchored reads per unit of Kcnd2 expression (a mutant value near half the wild-type one means the '
           'anchors are blind to the knock-in allele, NOT that the mutant mRNA is absent)')
        if {'mutant', 'wildtype'} <= set(s.index) and s['wildtype'] > 0:
            r = s['mutant'] / s['wildtype']
            print(f'\nmutant / wild-type anchored-read recovery: {r:.2f} '
                  + ('(both alleles are being counted)' if r > 0.8 else
                     'WARNING: mutant reads are being lost — redefine the anchors from the windows above and re-run '
                     'count_v404m_allele.sh'))
else:
    tot = None
    print('(no allele counts — run find_v404m_codon.py then count_v404m_allele.sh, and copy '
          'tables/v404m_allele_counts.tsv and tables/v404m_codon.tsv here)')

# ---------------------------------------------------------------- 2. Kv4-complex genes and variants
L = cpm_log(C.round().astype(int))
rows = []
for g in GENES:
    if g not in L.index:
        continue
    for t in tissues:
        idx = m2.index[m2.tissue == t]
        mut = L.loc[g, [i for i in idx if m2.at[i, 'genotype'] == 'mutant']]
        wt = L.loc[g, [i for i in idx if m2.at[i, 'genotype'] == 'wildtype']]
        rows.append((t, g, wt.mean(), mut.mean(), mut.mean() - wt.mean(), len(mut), len(wt), mw(mut, wt)))
G = pd.DataFrame(rows, columns=['tissue', 'gene', 'mean_logCPM_wt', 'mean_logCPM_mut', 'log2FC', 'n_mut', 'n_wt', 'p'])
G.to_csv(f'{OUT}/tables/v404m_kv4_genes.csv', index=False)
pp(G.pivot_table(index='gene', columns='tissue', values='log2FC').reindex(GENES)[tissues].round(3),
   f'Kv4-complex genes, log2FC mutant/WT ({source})')

iso_path = f'{TABLES}/v404m_isoform_per_sample.csv'
if os.path.exists(iso_path):
    I = pd.read_csv(iso_path).dropna(subset=['vclass'])
    I = I[~((I.biotype == 'retained_intron') & I.vclass.str.contains('other|fragment', case=False, na=False))]
    vs = I.groupby(['sample', 'tissue', 'genotype', 'gene', 'vclass'], as_index=False).agg(
        reads=('reads', 'sum'), gene_reads=('gene_reads', 'first'))
    vs = vs[vs.gene_reads >= 50]
    vs['frac'] = vs.reads / vs.gene_reads
    vr = []
    for (t, g, v), sub in vs.groupby(['tissue', 'gene', 'vclass']):
        a = sub[sub.genotype == 'mutant'].frac; b = sub[sub.genotype == 'wildtype'].frac
        if len(a) >= 3 and len(b) >= 3:
            vr.append((t, g, v, b.mean(), a.mean(), a.mean() - b.mean(), mw(a, b)))
    V = pd.DataFrame(vr, columns=['tissue', 'gene', 'vclass', 'frac_wt', 'frac_mut', 'diff', 'p'])
    if len(V):
        V['fdr'] = stats.false_discovery_control(V.p.fillna(1))
    V.to_csv(f'{OUT}/tables/v404m_variant_fracs.csv', index=False)
    pp(V[(V['diff'].abs() > 0.02) & (V.p < 0.05)].sort_values(['tissue', 'gene']).round(4) if len(V) else V,
       'auxiliary-subunit variant fractions differing by genotype (|diff| > 2%, p < 0.05)')
else:
    V = None

# ---------------------------------------------------------------- 3. Kcnd2 exon:intron
ei_path = f'{TABLES}/v404m_exon_intron.tsv'
if os.path.exists(ei_path):
    E = pd.read_csv(ei_path, sep='\t')
    s_best = E[(E.feature == 'exon') & (E.strand_setting > 0)].groupby('strand_setting')['count'].sum().idxmax()
    E = E[E.strand_setting == s_best].pivot_table(index=['sample', 'gene'], columns='feature', values='count').reset_index()
    E = E.merge(meta, left_on='sample', right_index=True)
    E['intron_frac'] = E.intron / (E.exon + E.intron).replace(0, np.nan)
    E.to_csv(f'{OUT}/tables/v404m_exon_intron_per_sample.csv', index=False)
    er = []
    for (t, g), sub in E.groupby(['tissue', 'gene']):
        a = sub[sub.genotype == 'mutant'].intron_frac; b = sub[sub.genotype == 'wildtype'].intron_frac
        er.append((t, g, b.mean(), a.mean(), a.mean() - b.mean(), mw(a, b)))
    EI = pd.DataFrame(er, columns=['tissue', 'gene', 'intron_frac_wt', 'intron_frac_mut', 'diff', 'p'])
    EI.to_csv(f'{OUT}/tables/v404m_exon_intron_tests.csv', index=False)
    pp(EI[EI.gene.isin(['Kcnd2', 'Kcnip4', 'Dpp6'])].round(4), f'intronic fraction of gene reads (strand setting -s {s_best})')
else:
    EI = None

# ---------------------------------------------------------------- 4. junction readouts
jp = f'{TABLES}/v404m_junctions_all.tsv.gz'
JR = None
if os.path.exists(jp):
    J = pd.read_csv(jp, sep='\t')
    J = J.merge(meta, left_on='sample', right_index=True, how='left')
    J['span'] = J.end - J.start + 1
    regions = {}
    rp = f'{TABLES}/v404m_regions.bed'
    if os.path.exists(rp):
        for line in open(rp):
            p = line.split('\t')
            if len(p) >= 4:
                regions.setdefault(p[3].strip(), []).append((p[0], int(p[1]), int(p[2])))

    def in_gene(df, gene):
        iv = regions.get(gene, [])
        if not iv:
            return df.iloc[0:0]
        m = np.zeros(len(df), bool)
        for c, s, e in iv:
            m |= (df.chr.values == c) & (df.start.values >= s) & (df.end.values <= e)
        return df[m]

    gene_tot = C.reindex(['Xbp1', 'Fosb']).fillna(0)
    out = []
    for gene, lo, hi, label in [('Xbp1', 20, 35, 'XBP1 spliced (IRE1, 26-nt intron)'),
                                ('Fosb', 120, 170, 'dFosB splice')]:
        g = in_gene(J, gene)
        if not len(g):
            continue
        cand = g[g.span.between(lo, hi)].groupby(['chr', 'start', 'end']).uniq.sum().sort_values(ascending=False)
        if not len(cand) or cand.iloc[0] < 20:
            print(f'\n{label}: no candidate junction with enough reads (max {int(cand.iloc[0]) if len(cand) else 0}) — '
                  f'this event cannot be measured in these libraries')
            continue
        (ch, st, en) = cand.index[0]
        ev = g[(g.chr == ch) & (g.start == st) & (g.end == en)].groupby('sample').uniq.sum()
        const = g[~((g.start == st) & (g.end == en))].groupby('sample').uniq.sum()
        d = pd.DataFrame({'event_reads': ev, 'other_junction_reads': const}).reindex(meta.index).fillna(0).join(meta)
        d['event_per_1k_gene_reads'] = d.event_reads / (gene_tot.loc[gene].reindex(d.index) / 1000).replace(0, np.nan)
        d['event_over_other_junctions'] = d.event_reads / d.other_junction_reads.replace(0, np.nan)
        d['event'] = label; d['locus'] = f'{ch}:{st}-{en} ({en - st + 1} nt)'
        out.append(d.reset_index().rename(columns={'index': 'sample'}))
        print(f'\n{label}: {ch}:{st}-{en} ({en - st + 1} nt), {int(cand.iloc[0]):,} reads pooled')
        summ = d.groupby(['tissue', 'genotype']).agg(n=('event_reads', 'size'), median_event_reads=('event_reads', 'median'),
                                                     median_per_1k=('event_per_1k_gene_reads', 'median'),
                                                     median_ratio=('event_over_other_junctions', 'median')).round(4)
        print(summ.to_string())
        for t in tissues:
            a = d[(d.tissue == t) & (d.genotype == 'mutant')].event_over_other_junctions
            b = d[(d.tissue == t) & (d.genotype == 'wildtype')].event_over_other_junctions
            p = mw(a.dropna(), b.dropna())
            if p == p:
                print(f'   {t}: mutant {a.median():.4f} vs WT {b.median():.4f}, p = {p:.3f}')

    # known cassette events from the dendrite analysis (GRCm39 / vM36 coordinates)
    for name, inc, skip in [('Kv4.3L exon inclusion', [('chr3', 105572863, 105574277), ('chr3', 105574335, 105575395)],
                             [('chr3', 105572863, 105575395)]),
                            ('Kcnip3 cys-rich acceptor', [('chr2', 127307895, 127352748)],
                             [('chr2', 127307899, 127352748)])]:
        def jsum(spec):
            tot = pd.Series(0.0, index=meta.index)
            for c, s, e in spec:
                v = J[(J.chr == c) & (J.start == s) & (J.end == e)].groupby('sample').uniq.sum()
                tot = tot.add(v.reindex(meta.index).fillna(0), fill_value=0)
            return tot / len(spec)
        i_, s_ = jsum(inc), jsum(skip)
        psi = (i_ / (i_ + s_).replace(0, np.nan)).rename('psi').to_frame().join(meta)
        psi['event'] = name
        if (i_ + s_).sum() < 50:
            print(f'\n{name}: too few junction reads in these libraries to measure'); continue
        out.append(psi.reset_index().rename(columns={'index': 'sample'}))
        print(f'\n{name}: PSI by tissue and genotype')
        print(psi.groupby(['tissue', 'genotype']).psi.agg(['median', 'size']).round(3).to_string())
        for t in tissues:
            p = mw(psi[(psi.tissue == t) & (psi.genotype == 'mutant')].psi.dropna(),
                   psi[(psi.tissue == t) & (psi.genotype == 'wildtype')].psi.dropna())
            if p == p:
                print(f'   {t}: p = {p:.3f}')
    if out:
        JR = pd.concat(out, ignore_index=True)
        JR.to_csv(f'{OUT}/tables/v404m_junction_events.csv', index=False)

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
ax = axes[0]
if tot is not None and len(tot):
    pref = 'haplotype' if 'haplotype' in set(tot.anchor) else ('right' if 'right' in set(tot.anchor) else 'left')
    d = tot[tot.anchor == pref]
    for j, t in enumerate(tissues):
        for gt, dx, col in [('wildtype', -0.17, '#4C72B0'), ('mutant', 0.17, '#C44E52')]:
            v = d[(d.tissue == t) & (d.genotype == gt)].alt_frac.dropna()
            ax.scatter(np.full(len(v), j + dx) + np.random.uniform(-0.04, 0.04, len(v)), v, s=14, color=col,
                       label=gt if j == 0 else None)
    ax.axhline(0.5, color='#999', lw=0.7, ls='--'); ax.set_ylim(-0.02, 1.02)
    ax.set_xticks(range(len(tissues))); ax.set_xticklabels([t[:4] for t in tissues])
    ax.set_ylabel('non-wild-type fraction of Kcnd2 codon-404 reads'); ax.legend(fontsize=7, frameon=False)
ax.set_title('A. Allele-specific Kcnd2', fontsize=10)
ax = axes[1]
sub = G[G.gene.isin(GENES)]
for j, g in enumerate(GENES):
    for t, mk in zip(tissues, 'os^D'):
        r = sub[(sub.gene == g) & (sub.tissue == t)]
        if len(r):
            ax.scatter(j, r.log2FC.iloc[0], marker=mk, s=22, label=t if j == 0 else None,
                       color='k' if r.p.iloc[0] >= 0.05 else '#C44E52')
ax.axhline(0, color='#999', lw=0.7); ax.set_xticks(range(len(GENES))); ax.set_xticklabels(GENES, rotation=90, fontsize=8)
ax.set_ylabel('log2FC mutant / WT'); ax.legend(fontsize=7, frameon=False); ax.set_title('B. Kv4-complex genes (red p<0.05)', fontsize=10)
ax = axes[2]
if EI is not None:
    e = pd.read_csv(f'{OUT}/tables/v404m_exon_intron_per_sample.csv')
    e = e[e.gene == 'Kcnd2']
    for j, t in enumerate(tissues):
        for gt, dx, col in [('wildtype', -0.17, '#4C72B0'), ('mutant', 0.17, '#C44E52')]:
            v = e[(e.tissue == t) & (e.genotype == gt)].intron_frac.dropna()
            ax.scatter(np.full(len(v), j + dx) + np.random.uniform(-0.04, 0.04, len(v)), v, s=14, color=col)
    ax.set_xticks(range(len(tissues))); ax.set_xticklabels([t[:4] for t in tissues])
    ax.set_ylabel('Kcnd2 intronic fraction')
ax.set_title('C. Kcnd2 exon:intron', fontsize=10)
fig.suptitle('Kv4.2-V404M vs wild-type littermates, 3 months', fontsize=11)
fig.tight_layout(); fig.savefig(f'{OUT}/figures/v404m_kv4.png', dpi=200)
print('\ndone')
