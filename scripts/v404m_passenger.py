#!/usr/bin/env python3
"""Is any of the V404M differential expression a genetic artefact rather than an effect of the mutation?

Two artefacts can masquerade as biology in a knock-in line, and they have different signatures:

 1. PASSENGER LOCUS. A knock-in targeted in 129-derived ES cells and backcrossed to B6 leaves a 129 interval
    flanking the edited gene. Genes inside that interval differ between genotypes because of strain
    polymorphism -- sequence divergence that costs reads at mapping, strain-specific Gm* loci, cis-eQTLs --
    not because of the mutation. Signature: changes CLUSTERED in a contiguous interval around Kcnd2
    (chr6:21.2-21.7 Mb) and CONSISTENT ACROSS ALL TISSUES, since a genetic difference is present in every cell.
    Consistency across tissues is the more specific half of that signature: real downstream effects of a
    channel mutation should vary with how much the channel matters in each region, whereas a strain variant
    does not care which tissue it is in.

 2. SEX IMBALANCE. With small groups, unequal sex composition puts every X- and Y-linked gene on the list.
    Signature: large fold-changes concentrated on chrX/chrY, tracking the sex split rather than the genotype.

The script tests both, and tests the clustering against a permutation null that slides the window to random
genomic positions -- necessary because genes are not uniformly distributed, so a fixed-width window around any
locus contains a variable number of genes and a raw count means nothing on its own.

Usage:  python3 scripts/v404m_passenger.py <tables_dir> [out_dir]
        tables_dir must hold v404m_de_<tissue>.csv, gene_positions.tsv and v404m_star_qc.csv
"""
import sys, os, glob
import numpy as np, pandas as pd
from scipy import stats

TAB = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
RNG = np.random.default_rng(0)
KCND2 = ('chr6', 21215502, 21729804)
WINDOWS = [1, 5, 10, 20, 40]          # Mb each side of the locus
N_PERM = 2000

# ---------------------------------------------------------------- inputs
pos = pd.read_csv(f'{TAB}/gene_positions.tsv', sep='\t', header=None,
                  names=['chr', 'start', 'end', 'gene'])
pos = pos.drop_duplicates('gene').set_index('gene')
de = {}
for f in sorted(glob.glob(f'{TAB}/v404m_de_*.csv')):
    t = os.path.basename(f)[len('v404m_de_'):-len('.csv')]
    d = pd.read_csv(f).drop_duplicates('gene').set_index('gene')
    de[t] = d
if not de:
    sys.exit(f'no v404m_de_*.csv in {TAB}')
tissues = list(de)
print(f'tissues: {", ".join(tissues)}')

L = pd.DataFrame({t: de[t]['log2FoldChange'] for t in tissues})
P = pd.DataFrame({t: de[t]['padj'] for t in tissues})
L = L.join(pos, how='inner')
P = P.reindex(L.index)
L['mean_lfc'] = L[tissues].mean(axis=1)
L['n_same_sign'] = (np.sign(L[tissues]).eq(np.sign(L['mean_lfc']), axis=0)
                    & L[tissues].abs().gt(0.5)).sum(axis=1)
L['mid'] = (L['start'] + L['end']) / 2
L['dist_kcnd2'] = np.where(L['chr'] == KCND2[0],
                           np.abs(L['mid'] - (KCND2[1] + KCND2[2]) / 2), np.nan)
print(f'{len(L):,} genes with positions and DE results')

# ---------------------------------------------------------------- 1. sex-chromosome check
print('\n' + '=' * 78)
print('1. SEX IMBALANCE')
print('=' * 78)
qc = pd.read_csv(f'{TAB}/v404m_star_qc.csv')
if 'sex' in qc.columns:
    tab = qc.pivot_table(index='tissue', columns=['genotype', 'sex'], values='sample', aggfunc='count')
    print(tab.fillna(0).astype(int).to_string())
auto = L[~L['chr'].isin(['chrX', 'chrY', 'chrM'])]
for name, sub in [('autosomes', auto), ('chrX', L[L['chr'] == 'chrX']), ('chrY', L[L['chr'] == 'chrY'])]:
    if len(sub) == 0:
        continue
    big = (sub[tissues].abs() > 1).any(axis=1).sum()
    print(f'{name:10} n={len(sub):6,}  |log2FC|>1 in any tissue: {big:4d} ({big/len(sub):.2%})  '
          f'median |mean log2FC| {sub["mean_lfc"].abs().median():.3f}')
ef = (L[L["chr"] == "chrX"][tissues].abs() > 1).any(axis=1).mean() / max(
     (auto[tissues].abs() > 1).any(axis=1).mean(), 1e-9)
print(f'chrX enrichment for large changes, relative to autosomes: {ef:.2f}x')

# ---------------------------------------------------------------- 2. clustering around Kcnd2
print('\n' + '=' * 78)
print('2. PASSENGER LOCUS AROUND Kcnd2 (chr6:21.2-21.7 Mb)')
print('=' * 78)
c6 = L[L['chr'] == KCND2[0]].copy()
print(f'chr6 genes with DE results: {len(c6):,}')
rows = []
for w in WINDOWS:
    near = c6[c6['dist_kcnd2'] <= w * 1e6]
    if len(near) == 0:
        continue
    obs = near['mean_lfc'].abs().median()
    nsig = int((P.loc[near.index, tissues] < 0.1).any(axis=1).sum())
    # permutation null: same number of consecutive genes elsewhere in the genome, preserving gene density
    order = L.sort_values(['chr', 'mid'])
    k = len(near)
    null_med, null_sig = np.empty(N_PERM), np.empty(N_PERM)
    for i in range(N_PERM):
        s = RNG.integers(0, len(order) - k)
        blk = order.iloc[s:s + k]
        null_med[i] = blk['mean_lfc'].abs().median()
        null_sig[i] = (P.loc[blk.index, tissues] < 0.1).any(axis=1).sum()
    p_med = (1 + np.sum(null_med >= obs)) / (1 + N_PERM)
    p_sig = (1 + np.sum(null_sig >= nsig)) / (1 + N_PERM)
    rows.append(dict(window_Mb=w, n_genes=k, median_abs_lfc=obs,
                     null_median=float(np.median(null_med)), p_lfc=p_med,
                     n_sig_padj10=nsig, null_n_sig=float(np.median(null_sig)), p_nsig=p_sig))
W = pd.DataFrame(rows)
print(W.round(4).to_string(index=False))

# cross-tissue-consistent strong changes: the specific passenger signature
cons = L[(L['mean_lfc'].abs() > 1) & (L['n_same_sign'] >= max(3, len(tissues) - 1))]
print(f'\ngenes with |mean log2FC| > 1 and same-direction change in >= {max(3, len(tissues)-1)} tissues: {len(cons)}')
if len(cons):
    byc = cons['chr'].value_counts()
    exp = L['chr'].value_counts().reindex(byc.index) / len(L) * len(cons)
    comp = pd.DataFrame({'observed': byc, 'expected': exp.round(1),
                         'ratio': (byc / exp).round(2)}).head(10)
    print(comp.to_string())
    n6 = cons[cons['chr'] == KCND2[0]]
    print(f'\non chr6: {len(n6)}')
    if len(n6):
        print(n6[['chr', 'mid', 'dist_kcnd2', 'mean_lfc', 'n_same_sign']]
              .assign(dist_Mb=lambda x: (x['dist_kcnd2'] / 1e6).round(2))
              .drop(columns='dist_kcnd2').sort_values('mid').to_string())
    print('\nthe 25 strongest cross-tissue-consistent changes and where they sit:')
    show = cons.reindex(cons['mean_lfc'].abs().sort_values(ascending=False).index).head(25)
    print(show[['chr', 'mid', 'mean_lfc', 'n_same_sign'] + tissues].round(2).to_string())

# ---------------------------------------------------------------- 3. the specific genes raised earlier
print('\n' + '=' * 78)
print('3. GENES FLAGGED EARLIER')
print('=' * 78)
for g in ['Gm14703', 'Onecut1', 'Kcnd2', 'Kcnd1', 'Gfap', 'Serpina3n']:
    if g in L.index:
        r = L.loc[g]
        d = '' if np.isnan(r['dist_kcnd2']) else f'  {r["dist_kcnd2"]/1e6:.1f} Mb from Kcnd2'
        print(f'{g:11} {r["chr"]:6} {int(r["mid"]):>11,}  mean log2FC {r["mean_lfc"]:+.2f}  '
              f'per tissue ' + ' '.join(f'{t[:4]}:{r[t]:+.2f}' for t in tissues) + d)

os.makedirs(f'{OUT}/tables', exist_ok=True)
L.to_csv(f'{OUT}/tables/v404m_gene_positions_lfc.csv')
W.to_csv(f'{OUT}/tables/v404m_passenger_windows.csv', index=False)
cons.to_csv(f'{OUT}/tables/v404m_cross_tissue_consistent.csv')
print(f'\nwrote {OUT}/tables/v404m_{{gene_positions_lfc,passenger_windows,cross_tissue_consistent}}.csv')
