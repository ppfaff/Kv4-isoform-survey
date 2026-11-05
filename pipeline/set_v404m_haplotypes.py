#!/usr/bin/env python3
"""Align the knock-in allele sequence to the reference and record both haplotypes for allele counting.

Knock-in alleles usually carry silent changes alongside the intended one (here a silent change creating an AgeI
site, used for genotyping). Those extra changes matter for read counting: any anchor overlapping one of them
stops matching the mutant allele. This script takes the mutant sequence, aligns it ungapped to the reference
context recorded by find_v404m_codon.py, and reports every difference with its codon, amino acid and whether it
is silent — then writes a pair of haplotype windows spanning all the differences, which the extractor counts
directly. Counting whole haplotypes is immune to the anchor problem, because each allele is matched by its own
sequence.

Usage:  python3 set_v404m_haplotypes.py <mutant sequence> [flank]
        python3 set_v404m_haplotypes.py CATTGCGCTACCGGTGCCAATGATCGTGTCG
Appends wt_haplotype / mut_haplotype (and the difference list) to tables/v404m_codon.tsv.
"""
import sys, os

COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
CODON_TABLE = {}
for _i, _b1 in enumerate('TCAG'):
    for _j, _b2 in enumerate('TCAG'):
        for _k, _b3 in enumerate('TCAG'):
            CODON_TABLE[_b1 + _b2 + _b3] = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'[_i * 16 + _j * 4 + _k]


def rc(s):
    return s.translate(COMP)[::-1]


def best_offset(ref, qry, min_overlap=18):
    """Ungapped alignment: offset of qry[0] relative to ref[0] maximising (matches - mismatches)."""
    best = None
    for off in range(-len(qry) + min_overlap, len(ref) - min_overlap + 1):
        a0, a1 = max(0, off), min(len(ref), off + len(qry))
        if a1 - a0 < min_overlap:
            continue
        r = ref[a0:a1]; q = qry[a0 - off:a1 - off]
        m = sum(x == y for x, y in zip(r, q))
        score = 2 * m - (len(r) - m)
        if best is None or score > best[0]:
            best = (score, off, a0, a1, m, len(r) - m)
    return best


def main(mut_seq, flank=4, tsv='tables/v404m_codon.tsv'):
    flank = int(flank)
    cfg = dict(l.rstrip('\n').split('\t', 1) for l in open(tsv) if '\t' in l)
    ref = cfg['mrna_window'].upper()                 # left_anchor + codon + right_anchor
    codon_start = len(cfg['left_anchor'])            # index of codon within ref
    mut = ''.join(c for c in mut_seq.upper() if c in 'ACGTN')

    cand = [('as given', mut), ('reverse complement', rc(mut))]
    scored = [(best_offset(ref, q), lab, q) for lab, q in cand]
    scored = [s for s in scored if s[0] is not None]
    if not scored:
        sys.exit('could not align the mutant sequence to the reference context')
    (score, off, a0, a1, nm, nmm), lab, mut = max(scored, key=lambda s: s[0][0])
    print(f'mutant sequence aligned {lab}: offset {off}, overlap {a1 - a0} nt, {nm} matches, {nmm} mismatches')

    diffs = []
    for i in range(a0, a1):
        r, q = ref[i], mut[i - off]
        if r != q:
            cs = codon_start + 3 * ((i - codon_start) // 3) if i >= codon_start else codon_start - 3 * (((codon_start - i) + 2) // 3)
            wt_cod = ref[cs:cs + 3]
            mt_cod = ''.join(mut[j - off] if a0 <= j < a1 else ref[j] for j in range(cs, cs + 3))
            res = 404 + (cs - codon_start) // 3
            diffs.append((i, r, q, res, wt_cod, mt_cod, CODON_TABLE.get(wt_cod, '?'), CODON_TABLE.get(mt_cod, '?')))
    if not diffs:
        sys.exit('no differences found — is this the wild-type sequence?')

    print(f'\n{len(diffs)} difference(s) from the reference:')
    print(f'{"ref_pos":>8} {"WT":>3} {"MUT":>4} {"residue":>8} {"WT codon":>9} {"MUT codon":>10} {"AA":>7}   effect')
    for i, r, q, res, wc, mc, wa, ma in diffs:
        eff = 'silent' if wa == ma else f'MISSENSE {wa}{res}{ma}'
        print(f'{i:>8} {r:>3} {q:>4} {res:>8} {wc:>9} {mc:>10} {wa:>3}->{ma:<3}   {eff}')

    lo = max(0, min(d[0] for d in diffs) - flank)
    hi = min(len(ref), max(d[0] for d in diffs) + 1 + flank)
    wt_hap = ref[lo:hi]
    mut_hap = ''.join(mut[j - off] if a0 <= j < a1 else ref[j] for j in range(lo, hi))
    print(f'\nhaplotype windows (reference positions {lo}-{hi - 1}, {hi - lo} nt):')
    print(f'  wild type : {wt_hap}')
    print(f'  knock-in  : {mut_hap}')
    print('  ' + ' ' * 12 + ''.join('^' if a != b else ' ' for a, b in zip(wt_hap, mut_hap)))

    la, ra = cfg['left_anchor'], cfg['right_anchor']
    la_hit = any(d[0] < codon_start for d in diffs)
    ra_hit = any(d[0] >= codon_start + 3 for d in diffs)
    print(f'\nanchor status: left {"BROKEN by a knock-in change — blind to the mutant allele" if la_hit else "intact"}; '
          f'right {"BROKEN by a knock-in change — blind to the mutant allele" if ra_hit else "intact"}')
    if la_hit and ra_hit:
        print('  both anchors are affected, so haplotype counting is the only valid readout here')

    with open(tsv, 'a') as f:
        f.write(f'wt_haplotype\t{wt_hap}\n')
        f.write(f'mut_haplotype\t{mut_hap}\n')
        f.write(f'wt_haplotype_rc\t{rc(wt_hap)}\n')
        f.write(f'mut_haplotype_rc\t{rc(mut_hap)}\n')
        f.write(f'haplotype_codon_offset\t{codon_start - lo}\n')
        f.write('haplotype_diffs\t' + ';'.join(f'{d[3]}:{d[4]}>{d[5]}({d[6]}>{d[7]})' for d in diffs) + '\n')
    print(f'\nappended wt_haplotype / mut_haplotype to {tsv}')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(*sys.argv[1:])
