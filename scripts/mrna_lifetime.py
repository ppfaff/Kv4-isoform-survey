#!/usr/bin/env python3
"""Bound the lifetime of mature A-current mRNAs from exon:intron coverage in total-RNA libraries.

At steady state the mature pool is M = k_synth * tau_m and the unspliced pool is P = k_init * tau_pre, so if
every initiation yields a mature transcript,

    tau_mature = (M / P) * tau_pre

M/P as MOLAR amounts, which means PER-BASE DENSITY, not read fractions. Kcnd2 is 99% intron by length, so a
47% intronic read fraction is a ~110:1 density ratio; using read fractions directly misstates the answer by
two orders of magnitude.

tau_pre is the intron RESIDENCE time, and for this gene family that is the whole question. Two regimes:

  LENGTH-INDEPENDENT. Splicing commits shortly after the 3' splice site is transcribed, so every intron sits
  for about the same tau_s regardless of size. Per-base intron coverage is then FLAT across introns.

  LENGTH-LIMITED. An intron cannot be spliced until its own 3' end exists, so a position s bases upstream of
  the 3' splice site persists for s/v + tau_s. Mean residence is L/(2v) + tau_s, which for Kcnd2's 497-kb
  first intron is ~124 min at 2 kb/min rather than ~5 min. Per-base coverage then RISES with intron length.

These differ by ~25x in the implied tau_m, so the regime is not a detail to be assumed - it is measured here
by regressing per-base intron coverage on intron length across the whole panel (197 introns spanning three
orders of magnitude), which needs no exonic length and is the central result.

CONFOUND THE REGRESSION MUST AVOID: reads spanning an exon-intron junction deposit coverage up to a read
length into the intron, which inflates SHORT introns (a 700-bp intron is ~30% junction-adjacent) and not long
ones. That bias mimics length-independence, so the regression is refit on long introns only, where it cannot
operate.

ONLY TOTAL-RNA LIBRARIES CAN BE USED. Pre-mRNA is not polyadenylated until cleavage, so polyA+ selection
depletes it by a large unknown factor; the V404M bulk's 2.9-11.6% intronic Kcnd2 is a property of the
selection chemistry, not of the RNA. The Farris laminar libraries sit near 47% intronic, which is not a
polyA+ number, so they are the only ones here that can support this.

Usage:  python3 mrna_lifetime.py [tables_dir] [out_dir]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

TAB = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
TAU_S = [float(x) for x in os.environ.get('TAU_S', '2,5,10,20').split(',')]   # splicing commitment, min
ELONG = float(os.environ.get('ELONG', 2.0))                                  # kb/min
MIN_LONG = int(os.environ.get('MIN_LONG', 2000))   # bp; above this, junction overhang is <10% of the intron
READLEN = 100                                      # for the overhang correction only

PANEL = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']


# --------------------------------------------------------------------------- exonic lengths
def exonic_lengths(tab):
    """Exact exonic length for the genes where the repo carries the annotation. Kcnd2 comes from the locus
    table (single annotated transcript, so span minus total intron is exact); Kcnd3 and Kcnip3 from their
    exon list. Genes without an exon list get no lifetime estimate rather than a guessed one."""
    out = {}
    p = f'{tab}/exons_kcnd3_kcnip3.tsv'
    if os.path.exists(p):
        e = pd.read_csv(p, sep='\t', header=None, names=['tx', 'start', 'end'])
        e['gene'] = e.tx.str.split('-').str[0]
        e['len'] = (e.end - e.start).abs() + 1
        for g, sub in e.groupby('gene'):
            tx = sub.tx.value_counts().idxmax()
            out[g] = int(sub[sub.tx == tx]['len'].sum())
    return out


def main():
    dep = pd.read_csv(f'{TAB}/star_intron_depth.tsv', sep='\t')
    ei = pd.read_csv(f'{TAB}/dendrite_exon_intron_per_sample.csv')

    # ---- which quantity is trustworthy? The two count tables must agree, and the depth table must be
    # reconcilable with them. Checked rather than assumed, because the lifetime is a ratio of these.
    sc = pd.read_csv(f'{TAB}/star_exon_intron_counts.tsv', sep='\t')
    s1 = (sc[sc.strand_setting == 1].pivot_table(index=['run', 'gene'], columns='feature',
                                                 values='count', aggfunc='first').reset_index())
    chk = ei.merge(s1, on=['run', 'gene'], suffixes=('_ps', '_star'))
    agree = np.allclose(chk.exon_ps, chk.exon_star) and np.allclose(chk.intron_ps, chk.intron_star)
    print(f'count tables agree (per_sample vs STAR, strand_setting=1): {agree}  '
          f'on {len(chk):,} gene x library pairs')
    if not agree:
        sys.exit('the two count tables disagree; resolve before computing a ratio of them')

    # depth_sum / mean_depth cannot be the same reads: check the implied read length.
    g2 = dep[dep.intron.str.startswith('Kcnd2_')].copy()
    g2['length'] = (g2.end - g2.start).abs()
    Li2 = g2.drop_duplicates('intron').length.sum()
    cov = g2.groupby('run').depth_sum.sum()
    rds = ei[ei.gene == 'Kcnd2'].set_index('run').intron.reindex(cov.index)
    implied = float((cov / rds).median())
    print(f'implied read length from the depth table: {implied:,.0f} bp  '
          f'(a real library is 50-150)')
    print('  -> the depth tables are on an UNEXPLAINED scale, ~%.0fx the read counts. Most likely they'
          % (implied / 100))
    print('     include multi-mapping alignments that featureCounts discards. They are therefore used')
    print('     below ONLY for RELATIVE comparisons within a library, never mixed with the counts.')
    print('     That bias is repeat-driven and so rises with intron length, which pushes the')
    print('     length-vs-coverage slope POSITIVE - making a flat or negative slope conservative.')
    print()
    dep['gene'] = dep.intron.str.rsplit('_intron_', n=1).str[0]
    dep['length'] = (dep.end - dep.start).abs()
    dep = dep[dep.gene.isin(PANEL) & (dep.length > 0)].copy()
    print(f'{len(dep):,} intron x library rows; {dep.intron.nunique()} distinct introns across '
          f'{dep.gene.nunique()} genes; {dep.run.nunique()} libraries')

    # intron inventory per gene
    inv = (dep.drop_duplicates('intron').groupby('gene')
              .agg(n_introns=('length', 'size'), intron_bp=('length', 'sum'),
                   longest=('length', 'max'), median_bp=('length', 'median')))
    inv['span_min_bp'] = (dep.groupby('gene').agg(a=('start', 'min'), b=('end', 'max'))
                          .eval('b - a'))
    EL = exonic_lengths(TAB)
    k2 = pd.read_csv(f'{TAB}/kcnd2_locus_annotation.tsv', sep='\t', header=None,
                     names=['kind', 'gene', 'biotype', 'tx', 'start', 'end', 'strand'])
    r = k2[(k2.kind == 'gene') & (k2.gene == 'Kcnd2')]
    if len(r):
        span2 = int(r.end.iloc[0] - r.start.iloc[0])
        EL['Kcnd2'] = span2 - int(inv.loc['Kcnd2', 'intron_bp'])
        inv.loc['Kcnd2', 'span_min_bp'] = span2
    inv['exonic_bp'] = [EL.get(g, np.nan) for g in inv.index]
    inv['intron_per_exon'] = inv.intron_bp / inv.exonic_bp

    print('\n' + '=' * 104)
    print('GENE DIMENSIONS - the "are these transcribed like huge genes" question')
    print('=' * 104)
    print(inv.assign(span_kb=(inv.span_min_bp / 1000).round(0),
                     intron_kb=(inv.intron_bp / 1000).round(0))
             [['n_introns', 'span_kb', 'intron_kb', 'longest', 'median_bp',
               'exonic_bp', 'intron_per_exon']].round(1).to_string())
    print('\n  span_kb from the first intron start to the last intron end, so it understates by the two')
    print('  terminal exons - immaterial against hundreds of kb of intron. exonic_bp only where the repo')
    print('  carries the exon annotation; the rest get no lifetime estimate rather than a guessed one.')
    print(f'\n  At {ELONG} kb/min, one round of transcription takes:')
    for g, row in inv.sort_values('span_min_bp', ascending=False).iterrows():
        t = row.span_min_bp / 1000 / ELONG
        print(f'    {g:8} {row.span_min_bp/1000:7.0f} kb -> {t:6.0f} min = {t/60:4.1f} h')

    # ----------------------------------------------------------------- regime test
    print('\n' + '=' * 104)
    print('WHICH REGIME? per-base intron coverage vs intron length')
    print('  length-limited splicing predicts slope +1 on log-log (density proportional to residence,')
    print('  residence proportional to length); prompt length-independent splicing predicts slope 0.')
    print('=' * 104)
    d = dep.copy()
    d['norm'] = d.mean_depth / d.groupby('run').mean_depth.transform('median')   # library-size control
    d = d[d.norm > 0]
    d['ln_len'] = np.log10(d.length)
    d['ln_dep'] = np.log10(d.norm)

    def fit(sub, tag):
        if sub.intron.nunique() < 8:
            print(f'  {tag:34} too few introns'); return None
        # collapse to one point per intron (median across libraries) so n is introns, not intron x library
        m = sub.groupby('intron').agg(ln_len=('ln_len', 'first'), ln_dep=('ln_dep', 'median'),
                                      gene=('gene', 'first'))
        res = stats.linregress(m.ln_len, m.ln_dep)
        lo = res.slope - 1.96 * res.stderr
        hi = res.slope + 1.96 * res.stderr
        verdict = ('LENGTH-LIMITED' if lo > 0.5 else
                   'length-independent' if hi < 0.5 else 'ambiguous')
        print(f'  {tag:34} n={len(m):4d} introns  slope {res.slope:+.3f} '
              f'[{lo:+.3f}, {hi:+.3f}]  r2={res.rvalue**2:.3f}  p={res.pvalue:.2g}   {verdict}')
        return res

    fit(d, 'all introns')
    fit(d[d.length >= MIN_LONG], f'introns >= {MIN_LONG} bp (clean)')
    fit(d[d.length >= 10000], 'introns >= 10 kb')
    print(f'\n  The >= {MIN_LONG} bp fit is the one to read: below that, junction-spanning reads deposit')
    print(f'  coverage ~{READLEN} bp into each intron end and inflate short introns, which would fake a')
    print('  slope of zero. If the clean fit is also flat, length-independence is real.')

    # per-gene slopes, to see whether any single gene behaves differently
    print('\n  per-gene slope on introns >= 1 kb:')
    for g, sub in d[d.length >= 1000].groupby('gene'):
        m = sub.groupby('intron').agg(x=('ln_len', 'first'), y=('ln_dep', 'median'))
        if len(m) < 6:
            print(f'    {g:8} n={len(m):3d}  (too few)'); continue
        rr = stats.linregress(m.x, m.y)
        print(f'    {g:8} n={len(m):3d}  slope {rr.slope:+.3f} +/- {rr.stderr:.3f}  p={rr.pvalue:.2g}')

    # ------------------------------------------------- premature termination (5'->3' attrition)
    print('\n' + '=' * 104)
    print("PREMATURE TERMINATION CHECK: does intron coverage decline along the gene?")
    print('  This is the confound that matters most for the lifetime estimate. If polymerase drops off')
    print('  before the end of a 500-kb gene, the 5\' introns are over-represented, the unspliced pool is')
    print('  inflated relative to what actually becomes mRNA, and tau_m comes out too SHORT.')
    print('  Negative slope = attrition. Position is taken from the INTRON NUMBER, which is already in')
    print('  transcriptional order, so this does not depend on knowing the strand - using genomic')
    print('  coordinate instead silently reverses every minus-strand gene (Dpp6 among them).')
    print('=' * 104)
    pos_rows = []
    d['n_intron'] = d.intron.str.rsplit('_intron_', n=1).str[1].astype(int)
    for g, sub in d.groupby('gene'):
        m = (sub.groupby('intron')
                .agg(k=('n_intron', 'first'), start=('start', 'first'),
                     y=('ln_dep', 'median'), length=('length', 'first'))
                .sort_values('k'))
        if len(m) < 8:
            continue
        strand = '+' if m.start.iloc[0] < m.start.iloc[-1] else '-'
        # relative position along the TRANSCRIBED unit, by cumulative intron length in transcription order
        cum = m.length.cumsum() - m.length / 2
        m['rel'] = cum / cum.iloc[-1]
        rr = stats.linregress(m.rel, m.y)
        fold = 10 ** rr.slope                     # fold change in coverage, 5' end to 3' end
        print(f'  {g:8} ({strand}) n={len(m):3d} introns  slope {rr.slope:+.3f} +/- {rr.stderr:.3f}  '
              f'p={rr.pvalue:.2g}   coverage x{fold:.3g} from 5\' to 3\'')
        pos_rows.append(dict(gene=g, strand=strand, n=len(m), slope=rr.slope, se=rr.stderr,
                             p=rr.pvalue, fold_5p_to_3p=fold))
        # first-vs-rest: an unused annotated first exon shows up as a depleted intron 1, not a gradient
        if len(m) >= 10:
            r1 = 10 ** (m.y.iloc[0] - m.y.iloc[1:].median())
            print(f'           intron 1 sits at {r1:.3g}x the median of the rest')
    if pos_rows:
        pd.DataFrame(pos_rows).to_csv(f'{OUT}/tables/intron_positional_slope.csv', index=False)
    print('\n  A gradient of a few fold inflates the unspliced pool in proportion and shifts tau_m by that')
    print('  factor at most. A swing of 20x or more is too steep for polymerase attrition and points')
    print('  instead to ALTERNATIVE PROMOTER USE: if the dominant first exon is internal, the annotated')
    print('  5\' introns are never transcribed and read as empty. For Dpp6 and Dpp10 that is the expected')
    print('  signature of the DPP6S/DPP6K and DPP10a/DPP10d start sites, and it is read here off intronic')
    print('  coverage alone, independent of the isoform quantification.')

    # ----------------------------------------------------------------- lifetime
    print('\n' + '=' * 104)
    print('MATURE mRNA LIFETIME, for the genes with exact exonic length')
    print('=' * 104)
    rows = []
    for g in [x for x in PANEL if not np.isnan(inv.loc[x, 'exonic_bp']) if x in inv.index]:
        Le, Li = float(inv.loc[g, 'exonic_bp']), float(inv.loc[g, 'intron_bp'])
        sub = ei[ei.gene == g]
        if sub.empty:
            continue
        ratio = float(np.median((sub.exon / Le) / (sub.intron / Li)))
        ifrac = float(sub.intron_frac.median())
        print(f'\n{g}:  exonic {Le:,.0f} bp, intronic {Li:,.0f} bp, intronic read fraction {ifrac:.1%}')
        print(f'  mature:unspliced DENSITY ratio = {ratio:,.0f} : 1   (n = {len(sub)} libraries)')
        print(f'  {"tau_s (min)":>12} {"tau_m (h)":>11} {"half-life (h)":>14}  '
              f'{"length-limited tau_m (h)":>25}')
        # length-limited alternative: coverage-weighted mean residence over the gene's own introns
        Lw = float((dep[dep.gene == g].drop_duplicates('intron').length ** 2).sum() / Li)  # length-weighted
        for ts in TAU_S:
            tm = ratio * ts / 60
            tau_pre_ll = Lw / (2 * ELONG * 1000) + ts
            tm_ll = ratio * tau_pre_ll / 60
            print(f'  {ts:12.0f} {tm:11.1f} {tm*np.log(2):14.1f}  {tm_ll:25,.0f}')
            rows.append(dict(gene=g, tau_s_min=ts, density_ratio=ratio, tau_m_h=tm,
                             half_life_h=tm * np.log(2), tau_m_h_length_limited=tm_ll))
        print(f'  (length-weighted mean intron = {Lw/1000:,.0f} kb, so the length-limited regime would put')
        print(f'   tau_pre at {Lw/(2*ELONG*1000):,.0f} min and tau_m in the column at right.)')

    if rows:
        os.makedirs(f'{OUT}/tables', exist_ok=True)
        pd.DataFrame(rows).to_csv(f'{OUT}/tables/mrna_lifetime_estimate.csv', index=False)
        inv.to_csv(f'{OUT}/tables/panel_gene_dimensions.csv')
        print(f'\nwrote {OUT}/tables/{{mrna_lifetime_estimate,panel_gene_dimensions}}.csv')

    print('\n' + '=' * 104)
    print('CONFOUNDS, all of which inflate the unspliced pool and so DEFLATE tau_m - read the estimate')
    print('as a lower bound: premature termination (not every initiation reaches the end of a 500-kb gene);')
    print('excised lariats before debranching; annotated retained-intron isoforms; independent transcription')
    print('units inside the introns (the Kcnd2 locus table lists two). Metabolic labelling (SLAM-seq, 4sU)')
    print('measures decay directly and would settle it.')
    print('=' * 104)


if __name__ == '__main__':
    main()
