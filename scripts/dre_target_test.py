#!/usr/bin/env python3
"""Are DREAM/KChIP transcriptional targets changed in V404M?

THE HYPOTHESIS BEING TESTED. KChIPs are not only channel subunits; all four bind downstream
regulatory element (DRE) sites and repress transcription (KChIP3/DREAM/calsenilin is the best
characterised, but KChIP1-4 all repress DRE-containing promoters). DNA binding requires the tetramer
and is divalent-cation sensitive: Mg2+ promotes it, and RISING Ca2+ PREVENTS it. So the direction is
specific and counterintuitive - DREAM represses at rest, and activity DEREPRESSES its targets. After
fear conditioning KChIP3 relocates from membrane to nucleus in hippocampal neurons, and sumoylation
favours the nuclear pool while palmitoylation and phosphorylation hold it at the membrane.

If a hypofunctional Kv4.2 raises excitability and Ca2+ entry, DREAM should come off DNA and its
targets should rise. Named brain DRE targets: Pdyn, Bdnf, Fos, Crem, Fosl2 (Fra-2); Cacnb2 is a
KChIP2 target. Aanat is a pineal target and is excluded.

WHY THIS TEST IS WEAK ON ITS OWN, AND HOW IT IS SHARPENED. Fos and Bdnf are activity-regulated by
many routes, so a mutant that is simply more excitable moves them whether or not DREAM is involved.
A bare "DRE targets are up" result is therefore uninformative. The test implemented here is
COMPETITIVE against a matched set of immediate-early genes that are activity-regulated but are NOT
known DRE targets (Egr1, Junb, Npas4, Arc, Nr4a1, Nr4a2, Dusp1, Btg2, Ier2). Only a DRE-specific
shift ON TOP of the general activity shift would point at this mechanism. Expect a negative result
to be the likely one, and read it as such.

A second, independent prediction is also checked. Deleting Kv4.2 reduces KChIP1-3 PROTEIN in mouse
brain, and loss of Kv4.2 in cortical pyramidal neurons causes targeted degradation of KChIP3 - the
channel protects its auxiliary subunit from degradation. V404M is a point mutation rather than a
null, so if it destabilises the complex the same thing could happen. That is a protein-level
prediction, so Kcnip mRNA being flat is the EXPECTED result and would not refute it; it is reported
to make the distinction explicit, because an unchanged transcript is exactly what a
post-translational mechanism looks like.

Usage:  python3 dre_target_test.py [tables_dir] [out_dir]
"""
import os, sys, glob
import numpy as np, pandas as pd
from scipy import stats

TAB = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'

DRE = ['Pdyn', 'Bdnf', 'Fos', 'Crem', 'Fosl2', 'Cacnb2']
IEG_CTRL = ['Egr1', 'Egr2', 'Junb', 'Jun', 'Npas4', 'Arc', 'Nr4a1', 'Nr4a2', 'Dusp1', 'Btg2', 'Ier2']
KCHIP = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']
KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']


def competitive(d, genes, label):
    """Wilcoxon rank-sum of the set's log2FC against all other genes. Competitive, not
    self-contained: the question is whether these genes moved MORE than the transcriptome did,
    which is the only version that is not answered by a global activity shift."""
    s = d[d.gene.isin(genes)].dropna(subset=['log2FoldChange'])
    r = d[~d.gene.isin(genes)].dropna(subset=['log2FoldChange'])
    if len(s) < 3:
        return None
    u, p = stats.mannwhitneyu(s.log2FoldChange, r.log2FoldChange, alternative='two-sided')
    return dict(set=label, n=len(s), median_log2FC=float(s.log2FoldChange.median()),
                median_rest=float(r.log2FoldChange.median()), p_competitive=float(p))


def main():
    files = sorted(glob.glob(f'{TAB}/v404m_de_*.csv'))
    if not files:
        sys.exit(f'no v404m_de_*.csv in {TAB}')
    rows, detail = [], []
    for f in files:
        tis = os.path.basename(f).replace('v404m_de_', '').replace('.csv', '')
        d = pd.read_csv(f)
        d['gene'] = d.gene.astype(str)
        print('=' * 100)
        print(f'{tis.upper()}   ({len(d):,} genes tested)')
        print('=' * 100)

        for genes, label in ((DRE, 'DRE targets'), (IEG_CTRL, 'IEG controls (not DRE)')):
            r = competitive(d, genes, label)
            if r:
                r['tissue'] = tis
                rows.append(r)
                print(f'  {label:24} n={r["n"]:2d}  median log2FC {r["median_log2FC"]:+.3f}  '
                      f'(rest {r["median_rest"]:+.3f})  competitive p = {r["p_competitive"]:.3f}')

        print(f'\n  {"gene":10} {"baseMean":>10} {"log2FC":>9} {"pvalue":>10} {"padj":>10}  set')
        for genes, label in ((DRE, 'DRE'), (IEG_CTRL, 'IEG'), (KCHIP, 'KChIP'), (KV4, 'Kv4')):
            for g in genes:
                s = d[d.gene == g]
                if s.empty:
                    print(f'  {g:10} {"":>10} {"":>9} {"":>10} {"":>10}  {label}  (not in table)')
                    continue
                r = s.iloc[0]
                star = '*' if (pd.notna(r.padj) and r.padj < 0.05) else ''
                print(f'  {g:10} {r.baseMean:10.1f} {r.log2FoldChange:+9.3f} {r.pvalue:10.3g} '
                      f'{r.padj if pd.notna(r.padj) else float("nan"):10.3g}  {label} {star}')
                detail.append(dict(tissue=tis, gene=g, set=label, baseMean=float(r.baseMean),
                                   log2FC=float(r.log2FoldChange), pvalue=float(r.pvalue),
                                   padj=float(r.padj) if pd.notna(r.padj) else np.nan))
        print()

    R = pd.DataFrame(rows)
    D = pd.DataFrame(detail)
    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/dre_target_competitive.csv', index=False)
    D.to_csv(f'{OUT}/tables/dre_target_detail.csv', index=False)

    print('=' * 100)
    print('READING THIS')
    print('=' * 100)
    if not R.empty:
        piv = R.pivot_table(index='tissue', columns='set', values='p_competitive')
        print('\ncompetitive p by tissue and set:')
        print(piv.round(3).to_string())
        sig = R[(R.set == 'DRE targets') & (R.p_competitive < 0.05)]
        if len(sig):
            print(f'\n  DRE targets shifted in: {", ".join(sig.tissue)}')
            print('  Check whether the IEG controls moved in the same tissues. If they did, this is a')
            print('  general activity shift and says nothing about DREAM. Only a DRE shift WITHOUT a')
            print('  matched IEG shift would support the mechanism.')
        else:
            print('\n  No tissue shows a DRE-specific shift. With 6 DRE genes the test has little')
            print('  power, so this is weak evidence against rather than a refutation - but there is')
            print('  no transcriptional footprint of KChIP derepression in this data.')
    print('\n  Kcnip mRNA: the protein-level prediction (channel loss destabilising KChIP) does NOT')
    print('  predict a transcript change, so flat Kcnip values here neither support nor refute it.')
    print('  Testing it needs a KChIP western or IP from V404M tissue, not RNA.')
    print(f'\nwrote {OUT}/tables/{{dre_target_competitive,dre_target_detail}}.csv')


if __name__ == '__main__':
    main()
