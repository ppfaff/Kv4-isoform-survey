#!/usr/bin/env python3
"""Place the Kv4-complex genes on the measured RNA-turnover landscape of mouse cortex.

Source: Qiu, Zhang, Xia, Gao et al., Nat Neurosci 29(11):2685-2698 (2026), PMID 42675226.
Supplementary Table 2 of that paper assigns 12,483 genes to 100 RNA kinetics modules, clustered on
synthesis rate (alpha), degradation rate (gamma) and total abundance across 15 cortical cell types,
from in vivo 4sU labelling.

WHAT THIS CAN AND CANNOT DO. The paper computed per-gene half-lives (Supplementary Fig. 14a reports
their distribution across 16 cortical cell types; Supplementary Fig. 19 compares 13,591 transcripts
between excitatory neurons and astrocytes), but those VALUES are not in the supplementary tables -
Table 1 is the ECS-regulated gene list and Table 2 is module membership only. So this script places
our genes by module, and calibrates each module against gene sets whose turnover class is already
known, rather than reporting half-lives the table does not contain. The actual numbers need the
journal's Source Data files for Fig. 4 / Supplementary Fig. 19, or recomputation from the H5AD
objects in GSE249398.

A CALIBRATION THAT DOES NOT WORK, AND IS REPORTED RATHER THAN QUIETLY DROPPED. Modules are not
ordered by speed - the paper's own examples are module #3 at 15.2 h and module #13 at 2.34 h - so
the module number says nothing on its own. The obvious fix is to calibrate each module against gene
sets of known turnover: immediate-early genes, whose mRNAs survive minutes, against ribosomal-protein
and oxidative-phosphorylation genes, among the most stable transcripts in any cell. That fails here,
and the paper's two stated modules are exactly the test that shows it: module #3 (15.2 h) contains
0 IEGs and 0 stable anchors, module #13 (2.34 h) contains 1 IEG and 0 stable anchors. With 24 IEGs
and 234 stable genes spread over 12,483 genes in modules of up to 651, the anchors are too sparse to
separate the one pair of modules whose answer is known, so any enrichment score built from them is
driven by its own pseudocount, not by data. An earlier version of this script reported such a score
and placed module #3 on the FAST side of its axis - contradicting the paper's own 15.2 h. The score
is gone. What survives is co-membership, which needs no calibration, plus a qualitative read of what
else is in a module.

Because the modules cluster on three properties, not on degradation alone, two genes with similar
half-lives can land in different modules if their synthesis or abundance differ. Co-membership is
therefore strong evidence of shared kinetics; non-co-membership is weak evidence of different
kinetics. Read the asymmetry that way.

Usage:  python3 halflife_modules.py <MOESM3_ESM.xlsx> [out_dir]
"""
import os, re, sys
import numpy as np, pandas as pd
from scipy import stats

XLSX = sys.argv[1] if len(sys.argv) > 1 else 'MOESM3_ESM.xlsx'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'

PANEL = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']

# Anchors whose turnover class is established independently of this paper.
IEG = ['Fos', 'Fosb', 'Fosl2', 'Jun', 'Junb', 'Jund', 'Egr1', 'Egr2', 'Egr3', 'Egr4', 'Npas4',
       'Arc', 'Nr4a1', 'Nr4a2', 'Nr4a3', 'Ier2', 'Ier5', 'Dusp1', 'Btg2', 'Atf3', 'Gadd45b',
       'Myc', 'Klf4', 'Sgk1', 'Per1', 'Bdnf', 'Nptx2', 'Homer1', 'Rgs2', 'Trib1']
STABLE_RX = re.compile(r'^(Rp[ls]\d|Rplp\d|Rps\d|mt-|Cox\d|Ndufa\d|Ndufb\d|Atp5|Uqcr)', re.I)
NEURONAL = ['Snap25', 'Syt1', 'Syn1', 'Nefl', 'Nefm', 'Camk2a', 'Camk2b', 'Calm1', 'Calm2',
            'Gria1', 'Gria2', 'Grin1', 'Grin2b', 'Scn2a', 'Kcnq2', 'Cacna1a', 'Hcn1', 'Actb',
            'Tubb3', 'Eef1a1', 'Eef2', 'Map2', 'Dlg4', 'Gabra1', 'Gad1', 'Slc17a7']


def load_modules(path):
    d = pd.read_excel(path, sheet_name='Supplementary Table 2', header=1)
    d = d.iloc[:, :3]
    d.columns = ['module', 'n_genes', 'genes']
    d = d.dropna(subset=['genes']).copy()
    d['module'] = d.module.astype(str).str.extract(r'(\d+)').astype(int)
    rows = []
    for _, r in d.iterrows():
        gs = [g.strip() for g in re.split(r'[,;\s]+', str(r.genes)) if g.strip()]
        for g in gs:
            rows.append((g, int(r.module)))
    M = pd.DataFrame(rows, columns=['gene', 'module']).drop_duplicates('gene')
    return d, M


def main():
    d, M = load_modules(XLSX)
    g2m = dict(zip(M.gene, M.module))
    print(f'Supplementary Table 2: {d.module.nunique()} modules, {len(M):,} genes')

    M['is_ieg'] = M.gene.isin(IEG)
    M['is_stable'] = M.gene.str.match(STABLE_RX)
    S = (M.groupby('module')
          .agg(n=('gene', 'size'), n_ieg=('is_ieg', 'sum'), n_stable=('is_stable', 'sum')))
    tot_ieg, tot_stable, N = int(M.is_ieg.sum()), int(M.is_stable.sum()), len(M)

    # ---- the one available validation: can the anchors separate the paper's two stated modules?
    print('\n' + '=' * 100)
    print('CALIBRATION ATTEMPT AND ITS VALIDATION')
    print('  The paper states module 3 = 15.2 h (slow) and module 13 = 2.34 h (fast). Any score that')
    print('  claims to place modules on a speed axis must separate those two.')
    print('=' * 100)
    print(f'  anchor sets available: {tot_ieg} IEGs and {tot_stable} ribosomal/OxPhos genes '
          f'out of {N:,} ({tot_ieg/N*100:.2f}%, {tot_stable/N*100:.1f}%)')
    for m, label in ((3, 'stated 15.2 h, SLOW'), (13, 'stated 2.34 h, FAST')):
        if m in S.index:
            s = S.loc[m]
            print(f'  module {m:3d} ({label:20}) n={int(s.n):4d}  '
                  f'IEGs={int(s.n_ieg)}  stable anchors={int(s.n_stable)}')
    # Fisher's exact on the 2x2 of IEG membership in the slow vs fast module, rather than an
    # arbitrary count threshold: the question is whether the anchors could have detected the
    # difference, not whether the two counts happen to differ.
    a, b = int(S.loc[3, 'n_ieg']), int(S.loc[3, 'n']) - int(S.loc[3, 'n_ieg'])
    c, d_ = int(S.loc[13, 'n_ieg']), int(S.loc[13, 'n']) - int(S.loc[13, 'n_ieg'])
    _, pfish = stats.fisher_exact([[a, b], [c, d_]])
    sep = pfish < 0.05
    print(f'\n  Fisher exact on IEG membership, module 3 vs module 13: p = {pfish:.2f}')
    print(f'  VERDICT: the anchors {"separate" if sep else "DO NOT separate"} the two known modules, '
          f'so {"a speed score could be built" if sep else "no speed score is reported"}.')
    print('  Module membership here places genes relative to EACH OTHER, not on an absolute scale.')

    rows = []
    for g in PANEL:
        m = g2m.get(g)
        if m is None:
            print(f'  {g:8} not in the table'); continue
        s = S.loc[m]
        rows.append(dict(gene=g, module=m, module_size=int(s.n)))
    P = pd.DataFrame(rows)
    print('\n' + '=' * 100)
    print('WHERE THE PANEL SITS')
    print('=' * 100)
    print(P.to_string(index=False))

    # ---- qualitative read: what else is in each panel gene's module
    print('\n  what each module is made of (families present, as a qualitative read only):')
    FAM = {'K channel': r'^Kcn', 'Na channel': r'^Scn', 'Ca channel': r'^Cacn', 'HCN': r'^Hcn',
           'Glu receptor': r'^(Gria|Grin|Grik|Grm)', 'GABA-A': r'^Gabr',
           'synaptic': r'^(Syn|Syt|Snap|Stx|Nrxn|Nlgn|Dlg|Shank|Homer|Cplx)',
           'ribosomal': r'^Rp[ls]\d'}
    for m in sorted(P.module.unique()):
        g = M[M.module == m].gene.tolist()
        bits = []
        for k, p in FAM.items():
            h = [x for x in g if re.match(p, x)]
            if h:
                bits.append(f'{k} {len(h)}')
        print(f'    module {m:3d} (n={len(g):4d}): {"; ".join(bits) or "none of the listed families"}')

    # ---- the shared module is the headline
    shared = P.groupby('module').gene.apply(list)
    print('\n' + '=' * 100)
    print('CO-MEMBERSHIP: panel genes sharing a kinetics module')
    print('=' * 100)
    for m, gs in shared.items():
        if len(gs) > 1:
            co = M[(M.module == m)].gene.tolist()
            notable = [x for x in co if x in NEURONAL]
            print(f'\n  module {m} ({int(S.loc[m, "n"])} genes) holds {len(gs)} of the panel: '
                  f'{", ".join(gs)}')
            print(f'    other well-characterised members: {", ".join(notable[:18]) or "none in my list"}')
            print(f'    module size {int(S.loc[m, "n"])}, IEGs {int(S.loc[m, "n_ieg"])}, '
                  f'ribosomal/OxPhos {int(S.loc[m, "n_stable"])}')
    singles = [f'{r.gene} (module {r.module})' for r in P.itertuples() if len(shared[r.module]) == 1]
    if singles:
        print(f'\n  on their own: {", ".join(singles)}')

    # ---- is co-membership more than chance?
    sizes = S.n / len(M)
    p_same = float((sizes ** 2).sum())        # chance two random genes share a module
    print(f'\n  chance that two randomly chosen genes share a module: {p_same:.4f} '
          f'(1 in {1/p_same:.0f})')
    n_pan = len(P)
    exp_pairs = p_same * n_pan * (n_pan - 1) / 2
    obs_pairs = sum(len(g) * (len(g) - 1) / 2 for g in shared)
    print(f'  panel pairs in the same module: observed {obs_pairs:.0f}, expected {exp_pairs:.2f}')

    print('\n' + '=' * 100)
    print('WHAT IS STILL MISSING, and where to get it')
    print('=' * 100)
    print('  Table 2 carries module MEMBERSHIP only. The per-gene half-lives exist in the paper')
    print('  (Supplementary Fig. 14a, Supplementary Fig. 19: 13,591 transcripts, Ex neurons vs')
    print('  astrocytes) but are not tabulated, so the module -> hours mapping cannot be closed from')
    print('  this file. Two routes: the journal Source Data for Fig. 4 and Supplementary Fig. 19, or')
    print('  recomputation from GSE249398 (2.0 Gb tar of H5AD, per-cell new and total counts;')
    print('  tau = -t / ln(1 - NTR) with t = 4 h, half-life = tau * ln2).')
    print('  NOTE THE AGE: GSE249398 samples are E16.5, P7 and P12 cortex, not adult.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    P.to_csv(f'{OUT}/tables/panel_kinetics_modules.csv', index=False)
    S.to_csv(f'{OUT}/tables/kinetics_module_anchors.csv')
    M.to_csv(f'{OUT}/tables/gene_kinetics_module.csv', index=False)
    print(f'\nwrote {OUT}/tables/{{panel_kinetics_modules,kinetics_module_anchors,'
          f'gene_kinetics_module}}.csv')


if __name__ == '__main__':
    main()
