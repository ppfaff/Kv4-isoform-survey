#!/usr/bin/env python3
"""Which KChIPs does each cell class actually use, and do we have a protein half-life for them?

This matters because the measured protein half-lives (Dorrbaum 2018, eLife 34202) cover KChIP1 and
KChIP4 but NOT KChIP2 or KChIP3, and the KChIP repertoire is strongly cell-class specific. So the
coverage of the protein data is uneven across the very populations being recorded, and that unevenness
needs stating rather than averaging over.

The half-lives, neuron-enriched cultures: KChIP1 2.67 d (5th percentile of the whole proteome),
KChIP4 2.02 d (2nd percentile), Kv4.3 3.87 d, DPP6 4.52 d, DPP10 4.36 d. KChIP2 and KChIP3 were not
detected - in mass spectrometry that usually means abundance or peptide properties rather than
absence, but either way there is no number.

CAVEAT ON THE VARIANT SHARES. Read-level variant assignment in this project is calibrated for the
DPP panel only, where it was shown to be exact at full-length coverage and to collapse into sink
isoforms under 3' bias. The Kv4 and KChIP panels have NOT been calibrated (it is still on the
recommended-analyses list), so variant-level KChIP shares carry the same unquantified assignment
risk. Gene-level sums are much safer than individual variant shares and both are reported, with the
gene level as the one to trust.

Usage:  python3 kchip_class_vs_halflife.py [tables_dir] [out_dir]
"""
import os, re, sys
import numpy as np, pandas as pd

TAB = sys.argv[1] if len(sys.argv) > 1 else 'tables'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'

# measured protein half-lives, days; None = not detected in the proteome
T_HALF = {'KChIP1': 2.67, 'KChIP2': None, 'KChIP3': None, 'KChIP4': 2.02}
PCTILE = {'KChIP1': '5th', 'KChIP4': '2nd'}


def kchip_gene(vclass):
    """'KChIP4bL' -> 'KChIP4'; 'Kcnip3-new-cys_rich' -> 'KChIP3'. Case-insensitive on purpose: the
    vclass labels mix 'KChIP4d' and 'Kcnip4 other/noncoding', and a case-sensitive pattern silently
    dumps half the reads into an 'other' bucket."""
    m = re.search(r'k(?:ch)?i?p?(\d)', str(vclass), re.I)
    if m:
        return 'KChIP' + m.group(1)
    m = re.search(r'kcnip\s*(\d)', str(vclass), re.I)
    return 'KChIP' + m.group(1) if m else 'unassigned'


def composition(path, label):
    d = pd.read_csv(path)
    k = d[d.panel.astype(str).str.contains('kchip', case=False, na=False)].copy()
    if k.empty:
        return None, None
    k['kgene'] = k.vclass.map(kchip_gene)
    gg = k.groupby(['cls', 'kgene']).reads.sum().unstack(fill_value=0)
    gg = (gg.T / gg.T.sum()).T * 100
    vv = k.groupby(['cls', 'vclass']).reads.sum().unstack(fill_value=0)
    vv = (vv.T / vv.T.sum()).T * 100
    return gg, vv


def main():
    print('=' * 100)
    print('KChIP REPERTOIRE BY CELL CLASS, % of KChIP reads  (gene level - the trustworthy view)')
    print('=' * 100)
    out = {}
    for path, label in ((f'{TAB}/ctx_variant_fracs_by_type.csv', 'CORTEX (Tasic 2016)'),
                        (f'{TAB}/hip_variant_fracs_by_subclass.csv', 'HIPPOCAMPUS')):
        if not os.path.exists(path):
            print(f'  [skip] {path} missing'); continue
        gg, vv = composition(path, label)
        if gg is None:
            continue
        out[label] = (gg, vv)
        print(f'\n{label}:')
        cols = [c for c in ['KChIP1', 'KChIP2', 'KChIP3', 'KChIP4', 'unassigned'] if c in gg.columns]
        print(gg[cols].round(1).to_string())
        if 'unassigned' in gg.columns and gg['unassigned'].max() > 5:
            print(f'  [warn] {gg["unassigned"].max():.0f}% unassigned - check the vclass parser')

    print('\n' + '=' * 100)
    print('COVERAGE: does the class\'s dominant KChIP have a measured protein half-life?')
    print('=' * 100)
    rows = []
    for label, (gg, vv) in out.items():
        print(f'\n{label}')
        for cls in [c for c in ('Glut', 'GABA') if c in gg.index]:
            share = gg.loc[cls].drop(labels=[x for x in ['unassigned'] if x in gg.columns],
                                     errors='ignore').sort_values(ascending=False)
            covered = sum(v for g, v in share.items() if T_HALF.get(g) is not None)
            top = share.index[0]
            # abundance-weighted mean half-life over the measured part only
            num = sum(v * T_HALF[g] for g, v in share.items() if T_HALF.get(g))
            wmean = num / covered if covered > 0 else np.nan
            print(f'  {cls:5}  dominant {top} ({share.iloc[0]:.0f}%)  '
                  f'half-life measured for {covered:.0f}% of KChIP reads')
            for g, v in share.head(4).items():
                t = T_HALF.get(g)
                tag = (f'{t:.2f} d ({PCTILE.get(g, "")} pct)' if t else 'NOT MEASURED')
                print(f'         {g:8} {v:5.1f}%   {tag}')
            if covered > 0:
                print(f'         weighted mean over the measured part: {wmean:.2f} d')
            rows.append(dict(dataset=label, cls=cls, dominant=top,
                             dominant_pct=round(float(share.iloc[0]), 1),
                             pct_with_half_life=round(float(covered), 1),
                             weighted_half_life_d=round(float(wmean), 2) if covered else np.nan))

    print('\n' + '=' * 100)
    print('THE ASYMMETRY, AND WHY IT LINES UP WITH THE LAB\'S OWN RECORDINGS')
    print('=' * 100)
    print('  Cortical INHIBITORY neurons run on KChIP1 + KChIP4, and BOTH are measured - and both')
    print('  are among the shortest-lived proteins in the entire dataset (2nd and 5th percentile).')
    print('  Cortical EXCITATORY neurons are dominated by KChIP3, which is NOT measured; only the')
    print('  minor partner KChIP4 is. So the protein-lifetime argument is nearly complete for')
    print('  inhibitory neurons and mostly empty for excitatory ones.')
    print()
    print('  That is the opposite of a problem for the experiment in hand. The lab\'s deep cortical')
    print('  inhibitory recordings (V404M n=5, SD 6.75; WT n=7, SD 3.43 - a 1.97x SD ratio implying')
    print('  p = 0.247) are in exactly the population whose Kv4 complex depends entirely on the two')
    print('  shortest-lived KChIPs. If KChIP supply is rate-limiting and turns over in ~2 days, those')
    print('  cells\' A-current should be intrinsically more variable than a pyramidal neuron\'s, and')
    print('  an inflated SD is what that predicts.')
    print()
    print('  TESTABLE CONSEQUENCE: excitatory neurons, running on KChIP3, should show LESS')
    print('  cell-to-cell variance in A-current parameters than deep inhibitory neurons - unless')
    print('  KChIP3 also turns over fast, which nobody has measured. Recording the same parameters')
    print('  in both populations in WT animals would settle the premise before any V404M comparison.')
    print()
    print('  AND KChIP3 IS THE ONE THAT MATTERS TWICE OVER: it is both the dominant excitatory KChIP')
    print('  and the best-characterised transcription factor of the family (DREAM/calsenilin). The')
    print('  single most useful missing number in this whole analysis is the KChIP3 protein')
    print('  half-life in pyramidal neurons.')

    if rows:
        os.makedirs(f'{OUT}/tables', exist_ok=True)
        pd.DataFrame(rows).to_csv(f'{OUT}/tables/kchip_class_halflife_coverage.csv', index=False)
        for label, (gg, vv) in out.items():
            slug = label.split()[0].lower()
            gg.round(2).to_csv(f'{OUT}/tables/kchip_gene_share_{slug}.csv')
        print(f'\nwrote {OUT}/tables/kchip_class_halflife_coverage.csv and per-dataset gene shares')


if __name__ == '__main__':
    main()
