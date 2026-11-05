#!/usr/bin/env python3
"""Gene panels for the V404M stress-response analysis, and the shared loading helpers.

The panels are deliberately split by branch rather than lumped into one "stress" set, because the branches make
different claims: a PERK/ATF4 signature without an ATF6 chaperone signature is an integrated stress response that
need not involve the ER at all, whereas chaperone induction with XBP1 splicing is unfolded-protein stress. The
apoptosis panel is split from the gliosis panel for the same reason — in brain, ongoing neuronal death is usually
easier to see as microglial and astrocytic reaction than as caspase transcript induction, which is largely
post-translationally regulated and often invisible in bulk RNA.
"""
import os
import numpy as np
import pandas as pd

PANELS = {
    # --- integrated stress response: eIF2alpha kinases and the ATF4 transcriptional output
    'ISR (ATF4 output)': ['Atf4', 'Ddit3', 'Trib3', 'Chac1', 'Atf3', 'Sesn2', 'Asns', 'Slc7a11', 'Mthfd2',
                          'Ppp1r15a', 'Cebpb', 'Eif4ebp1', 'Nupr1', 'Psat1', 'Phgdh', 'Shmt2', 'Cth', 'Pycr1',
                          'Slc3a2', 'Slc7a5', 'Xpot', 'Aars1', 'Cars1', 'Gars1', 'Sars1', 'Wars1', 'Yars1'],
    'ISR kinases': ['Eif2ak1', 'Eif2ak2', 'Eif2ak3', 'Eif2ak4', 'Eif2s1', 'Atf5'],
    # --- unfolded protein response, by branch
    'UPR chaperones (ATF6)': ['Hspa5', 'Hsp90b1', 'Pdia3', 'Pdia4', 'Pdia6', 'Dnajb9', 'Dnajb11', 'Sel1l',
                              'Edem1', 'Herpud1', 'Manf', 'Creld2', 'Sdf2l1', 'Calr', 'Canx', 'Hyou1'],
    'UPR IRE1/XBP1': ['Ern1', 'Xbp1', 'Sec61a1', 'Dnajc3', 'Srpr', 'Ssr1'],
    # --- death pathways
    'Apoptosis (intrinsic)': ['Bax', 'Bak1', 'Bbc3', 'Pmaip1', 'Bid', 'Bcl2', 'Bcl2l1', 'Mcl1', 'Apaf1',
                              'Casp3', 'Casp6', 'Casp7', 'Casp9', 'Aifm1', 'Endog', 'Diablo'],
    'p53 / DNA damage': ['Trp53', 'Cdkn1a', 'Mdm2', 'Gadd45a', 'Gadd45b', 'Bbc3', 'Ccng1', 'Eda2r', 'Phlda3'],
    'Necroptosis / inflammasome': ['Ripk1', 'Ripk3', 'Mlkl', 'Zbp1', 'Casp8', 'Casp1', 'Nlrp3', 'Pycard', 'Gsdmd'],
    'Autophagy / lysosome': ['Sqstm1', 'Map1lc3b', 'Gabarapl1', 'Becn1', 'Atg5', 'Atg7', 'Lamp1', 'Lamp2',
                             'Ctsb', 'Ctsd', 'Ctsz', 'Tfeb'],
    # --- the readouts that actually detect neuronal injury in bulk brain RNA
    'Astrogliosis': ['Gfap', 'Vim', 'Aqp4', 'S100b', 'Serpina3n', 'Lcn2', 'Cd44', 'Timp1', 'Osmr', 'A2m'],
    'Microglia / neuroinflammation': ['Aif1', 'Cd68', 'Itgam', 'Tyrobp', 'Trem2', 'C1qa', 'C1qb', 'C1qc', 'C3',
                                      'C4b', 'B2m', 'H2-D1', 'H2-K1', 'Cxcl10', 'Il1b', 'Tnf', 'Il6', 'Ccl2', 'Fcgr3'],
    'Disease-associated microglia': ['Apoe', 'Cst7', 'Clec7a', 'Itgax', 'Lpl', 'Axl', 'Spp1', 'Gpnmb', 'Ank'],
    'Neuronal integrity': ['Rbfox3', 'Snap25', 'Syp', 'Syn1', 'Camk2a', 'Nefl', 'Nefm', 'Nefh', 'Map2', 'Eno2', 'Tubb3'],
    # --- controls and context
    'Heat shock (cytosolic)': ['Hspa1a', 'Hspa1b', 'Dnajb1', 'Hsph1', 'Hspb1', 'Hspa8', 'Hsp90aa1', 'Bag3'],
    'Oxidative stress (NRF2)': ['Nfe2l2', 'Hmox1', 'Nqo1', 'Gclc', 'Gclm', 'Txnrd1', 'Srxn1', 'Sod1', 'Sod2', 'Cat', 'Gpx1'],
    'Immediate-early / activity': ['Fos', 'Fosb', 'Arc', 'Egr1', 'Egr2', 'Egr3', 'Egr4', 'Npas4', 'Junb',
                                   'Nr4a1', 'Nr4a2', 'Bdnf', 'Homer1', 'Dusp1', 'Ier2', 'Sgk1', 'Per1'],
    'Kv4 complex': ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10'],
    'Other K channels (compensation)': ['Kcna1', 'Kcna2', 'Kcna4', 'Kcnb1', 'Kcnc1', 'Kcnc2', 'Kcnq2', 'Kcnq3',
                                        'Kcnma1', 'Kcnn2', 'Kcnj2', 'Hcn1', 'Hcn2', 'Kcnab1', 'Kcnab2'],
}
STRESS_PANELS = ['ISR (ATF4 output)', 'ISR kinases', 'UPR chaperones (ATF6)', 'UPR IRE1/XBP1',
                 'Apoptosis (intrinsic)', 'p53 / DNA damage', 'Necroptosis / inflammasome', 'Autophagy / lysosome',
                 'Astrogliosis', 'Microglia / neuroinflammation', 'Disease-associated microglia']
SEX_GENES = ['Xist', 'Ddx3y', 'Uty', 'Eif2s3y', 'Kdm5d']


def load_counts(tables='tables', prefer='star'):
    """Return (counts DataFrame genes x samples, source string). Prefers the STAR matrix when present."""
    star = os.path.join(tables, 'v404m_star_counts.csv.gz')
    salmon = os.path.join(tables, 'v404m_gene_counts.csv.gz')
    if prefer == 'star' and os.path.exists(star):
        C = pd.read_csv(star, index_col=0)
        m = os.path.join(tables, 'v404m_gene_id_map.csv')
        if os.path.exists(m):
            gm = pd.read_csv(m).drop_duplicates('gene_id').set_index('gene_id').symbol
            C.index = [gm.get(str(i).split('.')[0], str(i)) for i in C.index]
            C = C.groupby(level=0).sum()
        return C, 'STAR ReadsPerGene'
    C = pd.read_csv(salmon, index_col=0)
    return C.groupby(level=0).sum(), 'salmon (gene-summed NumReads)'


def load_meta(tables='tables'):
    p = os.path.join(tables, 'v404m_samples.tsv')
    if os.path.exists(p):
        S = pd.read_csv(p, sep='\t')
    else:
        S = pd.read_csv(os.path.join(tables, 'v404m_salmon_qc.csv'))
    S = S.rename(columns={'ID': 'sample'})
    return S.set_index('sample')[['tissue', 'genotype', 'sex']]


def cpm_log(C, lib=None):
    lib = C.sum() if lib is None else lib
    return np.log2(C / lib * 1e6 + 1)
