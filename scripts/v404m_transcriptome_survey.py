#!/usr/bin/env python3
"""What else changes transcriptionally in V404M, and is it region-specific?

The Kv4-complex work in this project has been mostly about occupancy and turnover, which are
post-translational. This script asks the separate question of what the transcriptome does, treating a
transcript change as a different mechanism rather than as evidence for the same one.

Three guards are built in, because the obvious readings here are the wrong ones.

  COMPOSITION. Cortical gene-level calls in this dataset are confounded by a cell-composition shift
  (established earlier: the apparent Kv4 decrease is fewer neurons, not less transcription). So any
  panel that moves in cortex is tested against a neuron-identity panel and a glia-identity panel. A
  panel that falls no faster than neuron identity itself is dilution; one that falls several-fold
  faster is not.

  COMPETITIVE, NOT SELF-CONTAINED. Panels are tested against all other genes in the same tissue, so a
  global shift in the libraries cannot produce a positive result.

  REGION SPECIFICITY IS THE POINT. A change in one tissue and not three others is far more
  interpretable than a change everywhere, and the four tissues here differ enormously in how much
  moves at all.

Usage:  python3 v404m_transcriptome_survey.py [tables_dir] [out_dir]
"""
import os, sys
import numpy as np, pandas as pd
from scipy import stats

TAB = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
TIS = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']

PANELS = {
    'immature/adult-born granule': ['Sox11', 'Sox4', 'Dcx', 'Igfbpl1', 'Calb2', 'Neurod1', 'Eomes',
                                    'Ascl1', 'Nes', 'Mki67', 'Cd24a', 'Stmn1', 'Tubb3', 'Bhlhe22'],
    'mature granule': ['Calb1', 'C1ql2', 'Prox1', 'Trpc6', 'Rasgrf1', 'Dsp', 'Pdzd2'],
    'K channels': ['Kcnq2', 'Kcnq3', 'Kcnq5', 'Kcna1', 'Kcna2', 'Kcnb1', 'Kcnc1', 'Kcnd2', 'Kcnd3',
                   'Kcnma1', 'Kcnn2', 'Kcnj10', 'Hcn1', 'Hcn2'],
    'Na/Ca channels': ['Scn1a', 'Scn2a', 'Scn8a', 'Cacna1a', 'Cacna1b', 'Cacna1e', 'Cacna1g', 'Cacnb4'],
    'IEG / ERK targets': ['Egr1', 'Egr2', 'Fos', 'Junb', 'Npas4', 'Arc', 'Nr4a1', 'Dusp1', 'Btg2'],
    'neuron identity': ['Rbfox3', 'Snap25', 'Syt1', 'Syn1', 'Nefl', 'Nefm', 'Map2', 'Tubb3', 'Eno2',
                        'Stmn2', 'Nrgn', 'Camk2a', 'Slc17a7', 'Gad1'],
    'glia identity': ['Gfap', 'Aqp4', 'Slc1a2', 'Slc1a3', 'Gja1', 'Aldoc', 'Cx3cr1', 'Mbp', 'Plp1',
                      'Sox10', 'Pdgfra', 'Cspg4'],
}


def competitive(d, genes):
    s = d[d.gene.isin(genes)].dropna(subset=['log2FoldChange'])
    r = d[~d.gene.isin(genes)].dropna(subset=['log2FoldChange'])
    if len(s) < 3:
        return None
    _, p = stats.mannwhitneyu(s.log2FoldChange, r.log2FoldChange, alternative='two-sided')
    return len(s), float(s.log2FoldChange.median()), float(p)


def main():
    D = {t: pd.read_csv(f'{TAB}/v404m_de_{t}.csv') for t in TIS}

    print('=' * 100)
    print('HOW MUCH MOVES AT ALL, BY TISSUE')
    print('=' * 100)
    print(f'  {"tissue":12} {"tested":>8} {"padj<0.05":>10} {"padj<0.01":>10} {"big & sig":>10}')
    for t, d in D.items():
        print(f'  {t:12} {len(d):8,} {int((d.padj < 0.05).sum()):10,} '
              f'{int((d.padj < 0.01).sum()):10,} '
              f'{int(((d.padj < 0.05) & (d.log2FoldChange.abs() > 1)).sum()):10,}')
    print('\n  The asymmetry is the first result: the forebrain regions with the astrogliosis are the')
    print('  ones with a transcriptome response, and thalamus and cerebellum have essentially none.')

    print('\n' + '=' * 100)
    print('PANELS, competitive against all other genes in the same tissue')
    print('=' * 100)
    rows = []
    print(f'  {"panel":30} ' + ' '.join(f'{t[:11]:>19}' for t in TIS))
    for lab, genes in PANELS.items():
        cells = []
        for t in TIS:
            r = competitive(D[t], genes)
            cells.append(f'{r[1]:+.3f} p={r[2]:.1e}' if r else '        --         ')
            if r:
                rows.append(dict(panel=lab, tissue=t, n=r[0], median_log2FC=round(r[1], 4),
                                 p_competitive=r[2]))
        print(f'  {lab:30} ' + ' '.join(f'{c:>19}' for c in cells))

    # ---- the composition guard, applied to the channel result
    print('\n' + '=' * 100)
    print('IS A CHANNEL DECREASE JUST NEURONAL DILUTION?')
    print('  dilution moves channels and neuron-identity genes by the SAME factor.')
    print('=' * 100)
    R = pd.DataFrame(rows)
    for t in TIS:
        g = R[R.tissue == t].set_index('panel').median_log2FC
        if 'neuron identity' not in g:
            continue
        ni = g['neuron identity']
        for pan in ('K channels', 'Na/Ca channels'):
            if pan in g:
                ratio = g[pan] / ni if abs(ni) > 1e-6 else np.inf
                verdict = ('consistent with dilution' if 0.5 < ratio < 2 else
                           f'{ratio:.1f}x the neuron-identity shift — NOT dilution alone'
                           if g[pan] < 0 and ni < 0 else 'opposite sign to dilution')
                print(f'  {t:12} {pan:16} {g[pan]:+.3f} vs neuron identity {ni:+.3f}   {verdict}')

    # ---- multi-tissue genes
    print('\n' + '=' * 100)
    print('GENES SIGNIFICANT IN MULTIPLE TISSUES')
    print('=' * 100)
    sig = {t: set(d[d.padj < 0.05].gene) for t, d in D.items()}
    multi = sorted(g for g in set().union(*sig.values())
                   if sum(g in sig[t] for t in TIS) >= 3)
    if not multi:
        print('  none in three or more')
    for g in multi:
        print(f'  {g}:')
        for t in TIS:
            r = D[t][D[t].gene == g]
            if len(r):
                x = r.iloc[0]
                print(f'    {t:12} log2FC {x.log2FoldChange:+.3f}  baseMean {x.baseMean:9.0f}  '
                      f'padj {x.padj:.3g}')
    for i, a in enumerate(TIS):
        for b in TIS[i + 1:]:
            print(f'  overlap {a} n {b}: {len(sig[a] & sig[b])}')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    R.to_csv(f'{OUT}/tables/v404m_panel_survey.csv', index=False)
    top = pd.concat([D[t].assign(tissue=t).nsmallest(40, 'padj') for t in TIS])
    top.to_csv(f'{OUT}/tables/v404m_top_genes_by_tissue.csv', index=False)
    print(f'\nwrote {OUT}/tables/{{v404m_panel_survey,v404m_top_genes_by_tissue}}.csv')


if __name__ == '__main__':
    main()
