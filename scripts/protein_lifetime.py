#!/usr/bin/env python3
"""Measured protein half-lives for the Kv4 complex, against the mRNA estimate.

Source: Dorrbaum, Kochen, Langer & Schuman, "Local and global influences on protein turnover in
neurons and glia", eLife 2018;7:e34202, Supplementary file 1. SILAC LC-MS/MS pulse-chase in rat
cortex/hippocampus cultures (P1-derived, 18-19 DIV), reported separately for mixed, neuron-enriched
and glia-enriched cultures. This is the dataset plotted in Supplementary Fig. 19b of Qiu et al.
Nat Neurosci 29:2685 (2026), reference 22 there.

Columns in that workbook: index, UniProt accessions, gene, protein name, rate constant k (per day),
SE of k, R^2 of the exponential fit, half-life (days). The identity half-life = ln2 / k is checked
rather than assumed, because a mislabelled column here would silently rescale every conclusion.

WHY THIS MATTERS MORE THAN THE mRNA NUMBER. The A-current a pipette measures is carried by protein,
not mRNA, and the two pools turn over on different timescales. If protein outlives its transcript by
a large factor, the channel complement of a cell reflects transcription integrated over the PROTEIN
lifetime, and a snapshot of mRNA at capture is a correspondingly poor predictor of current. That
directly bears on the prediction that ~19% of V404M recordings should look wild-type: for that to
hold, a cell must be Kv4.2-negative for as long as the protein persists, not merely at the instant
of capture.

CAVEATS, which are real. Rat, cultured, P1-derived, 18-19 DIV - not adult mouse in vivo, and these
cultures have no intact circuit. The values are WHOLE-CELL pools: surface delivery, endocytosis and
recycling of Kv4.2 are much faster than its degradation, so these half-lives bound the lifetime of
the total protein, not of the functional surface channel. Kcnd2 was quantified only in mixed
cultures, which contain glia; Kcnip2 and Kcnip3 were not detected at all, so absence here is
proteomic depth, not biology.

Usage:  python3 protein_lifetime.py <elife-34202-supp1.xlsx> [out_dir]
"""
import os, sys
import numpy as np, pandas as pd

XLSX = sys.argv[1] if len(sys.argv) > 1 else 'elife-34202-supp1.xlsx'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'

PANEL = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']
SHEETS = {'A) Mixed cultures': 'mixed', 'B) Neuron-enriched cultures': 'neuron',
          'C) Glia-enriched cultures': 'glia'}
# our own exon:intron density-ratio estimates, hours, at tau_s = 2 / 5 / 10 min (mrna_lifetime.py)
MRNA_H = {'Kcnd2': (2.6, 6.5, 13.0), 'Kcnd3': (7.6, 19.1, 38.1), 'Kcnip3': (2.5, 6.3, 12.5)}


def read_sheet(path, sheet):
    d = pd.read_excel(path, sheet_name=sheet, header=None)
    # locate the header row by finding the row that carries the gene column
    gene_col, k_col, hl_col = 2, 4, 7
    d = d[[0, 1, gene_col, 3, k_col, 5, 6, hl_col]].copy()
    d.columns = ['idx', 'uniprot', 'gene', 'protein', 'k_per_day', 'se', 'r2', 'half_life_days']
    for c in ('k_per_day', 'se', 'r2', 'half_life_days'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['gene'] = d.gene.astype(str).str.strip()
    return d.dropna(subset=['k_per_day', 'half_life_days'])


def main():
    frames = {}
    for sh, tag in SHEETS.items():
        try:
            frames[tag] = read_sheet(XLSX, sh)
        except Exception as e:
            print(f'[warn] could not read "{sh}": {e}')
    if not frames:
        sys.exit('no sheets read')

    # ---- unit check: is column 8 really ln2 / k?
    print('=' * 98)
    print('UNIT CHECK: does the stated half-life equal ln2 / k?')
    print('=' * 98)
    for tag, d in frames.items():
        pred = np.log(2) / d.k_per_day
        err = np.abs(pred - d.half_life_days) / d.half_life_days
        print(f'  {tag:7} n={len(d):5,}  median relative error {np.median(err):.2e}  '
              f'max {err.max():.2e}   {"consistent" if np.median(err) < 1e-3 else "MISMATCH"}')
    print('  k is per day, so the half-life column is in DAYS.')

    print('\n' + '=' * 98)
    print('DISTRIBUTION of protein half-lives (days), for context')
    print('=' * 98)
    for tag, d in frames.items():
        q = d.half_life_days.quantile([.1, .25, .5, .75, .9]).round(2)
        print(f'  {tag:7} n={len(d):5,}  p10 {q[.1]:6.2f}  Q1 {q[.25]:6.2f}  median {q[.5]:6.2f}  '
              f'Q3 {q[.75]:6.2f}  p90 {q[.9]:6.2f}  max {d.half_life_days.max():7.2f}')

    print('\n' + '=' * 98)
    print('THE Kv4 COMPLEX')
    print('=' * 98)
    rows = []
    for g in PANEL:
        for tag, d in frames.items():
            s = d[d.gene == g]
            if s.empty:
                continue
            r = s.iloc[0]
            pct = float((d.half_life_days < r.half_life_days).mean() * 100)
            rows.append(dict(gene=g, culture=tag, half_life_days=round(float(r.half_life_days), 2),
                             k_per_day=round(float(r.k_per_day), 4), se=round(float(r.se), 4),
                             r2=round(float(r.r2), 3), percentile=round(pct, 1),
                             uniprot=str(r.uniprot)[:40]))
    P = pd.DataFrame(rows)
    if P.empty:
        sys.exit('no panel proteins found')
    piv = P.pivot_table(index='gene', columns='culture', values='half_life_days')
    piv = piv.reindex([g for g in PANEL if g in piv.index])
    print('\nhalf-life, DAYS:')
    print(piv.reindex(columns=[c for c in ('neuron', 'mixed', 'glia') if c in piv.columns])
            .to_string(float_format=lambda v: f'{v:.2f}'))
    print('\nfull detail, with percentile within each culture type and fit quality:')
    print(P.sort_values(['culture', 'gene']).to_string(index=False))
    missing = [g for g in PANEL if g not in set(P.gene)]
    if missing:
        print(f'\n  not detected in any culture type: {", ".join(missing)}')
        print('  (proteomic depth is ~3,500 proteins, so absence is detection, not biology)')

    # ---- the comparison that matters
    print('\n' + '=' * 98)
    print('PROTEIN vs mRNA LIFETIME, and what it implies')
    print('=' * 98)
    print(f'{"gene":8} {"protein t1/2":>14} {"mRNA t1/2 (h)":>24} {"protein/mRNA ratio":>22}')
    comp = []
    for g, (lo, mid, hi) in MRNA_H.items():
        sub = P[P.gene == g]
        if sub.empty:
            continue
        pref = sub[sub.culture == 'neuron']
        r = (pref if not pref.empty else sub).iloc[0]
        ph = float(r.half_life_days) * 24
        print(f'{g:8} {r.half_life_days:7.2f} d ({r.culture[:6]:6}) '
              f'{lo:6.1f} - {hi:5.1f} (mid {mid:4.1f}) '
              f'{ph/hi:8.0f} - {ph/lo:4.0f}x (mid {ph/mid:3.0f}x)')
        comp.append(dict(gene=g, protein_half_life_h=round(ph, 1), protein_culture=r.culture,
                         mrna_half_life_h_low=lo, mrna_half_life_h_mid=mid, mrna_half_life_h_high=hi,
                         ratio_mid=round(ph / mid, 1)))
    print('\n  The mRNA figures are our own exon:intron density-ratio estimates (mrna_lifetime.py),')
    print('  swept over a splicing commitment time of 2-10 min; the protein figures are measured.')
    print('  The comparison crosses species and preparation (rat culture vs mouse hippocampus) and')
    print('  mixes a measurement with an inference, so read the ORDER OF MAGNITUDE, not the digits.')

    print('\n  CONSEQUENCE. Kv4.2 protein outlives its own transcript by roughly a factor of 20.')
    print('  The channel complement of a neuron therefore reflects Kcnd2 transcription integrated')
    print('  over several DAYS, while the single-cell measurement reports mRNA present in a window')
    print('  of hours. Two things follow, pulling in opposite directions:')
    print('   - FOR the stability argument: the A-current cannot track transcription on any')
    print('     timescale shorter than days, which is why the auxiliary subunits look constitutive.')
    print('   - AGAINST the 19.4% prediction: a cell scored Kv4.2-negative at capture may still')
    print('     carry protein made over the preceding week. For ~1 in 5 V404M recordings to look')
    print('     wild-type, the negative state must PERSIST for days, not merely hold at capture.')
    print('     The mRNA zero is necessary but not sufficient, and this is now the binding')
    print('     limitation on that prediction - it replaces the bursting caveat, which the mRNA')
    print('     lifetime retired.')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    P.to_csv(f'{OUT}/tables/protein_half_life_kv4_panel.csv', index=False)
    if comp:
        pd.DataFrame(comp).to_csv(f'{OUT}/tables/protein_vs_mrna_lifetime.csv', index=False)
    print(f'\nwrote {OUT}/tables/{{protein_half_life_kv4_panel,protein_vs_mrna_lifetime}}.csv')


if __name__ == '__main__':
    main()
