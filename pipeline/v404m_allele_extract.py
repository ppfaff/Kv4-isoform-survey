#!/usr/bin/env python3
"""Read out the sequence at a codon directly from raw reads, given flanking anchors (no alignment).

Reads anchor definitions from tables/v404m_codon.tsv and FASTQ sequence lines from stdin (normally pre-filtered
with `grep -F`, so only candidate reads arrive here). Each anchor reports the codon *and the flank on its far
side*, so the two anchors between them see both sides of the codon:

    left/fwd   read contains LEFT_ANCHOR       -> codon + FLANK bases downstream
    left/rev   read contains rc(LEFT_ANCHOR)   -> same, from a reverse-complement read
    right/fwd  read contains RIGHT_ANCHOR      -> FLANK bases upstream + codon
    right/rev  read contains rc(RIGHT_ANCHOR)  -> same, from a reverse-complement read

Knock-in alleles usually carry silent changes besides the intended one (here, one introducing a diagnostic
restriction site). Such a change may fall inside an anchor, in which case that anchor stops matching the mutant
allele and would under-report it. Two independent anchors make that visible: the affected one loses roughly half
its reads in heterozygotes while the other still sees both alleles, and the captured flank from the intact anchor
spells out the real mutant haplotype, which can then be used to redefine the anchors exactly.

Usage:  zcat R1.fastq.gz R2.fastq.gz | grep -F -e ... | python3 v404m_allele_extract.py SAMPLE [codon.tsv] [FLANK]
Prints TSV rows: sample, anchor, orientation, triplet, count, window
  (window is the captured sequence; the codon is its first 3 bases for the left anchor, its last 3 for the right)
"""
import sys

COMP = str.maketrans('ACGTN', 'TGCAN')
FLANK = 12


def rc(s):
    return s.translate(COMP)[::-1]


def main(sample, tsv='tables/v404m_codon.tsv', flank=FLANK):
    flank = int(flank)
    width = 3 + flank
    cfg = dict(l.rstrip('\n').split('\t', 1) for l in open(tsv) if '\t' in l)
    L, R = cfg['left_anchor'], cfg['right_anchor']
    Lrc, Rrc = rc(L), rc(R)
    counts = {}

    def add(anchor, orient, win):
        if len(win) == width and set(win) <= set('ACGT'):
            key = (anchor, orient, win)
            counts[key] = counts.get(key, 0) + 1

    # whole-allele haplotypes, when set_v404m_haplotypes.py has recorded them: each allele is matched by its own
    # sequence, so a knock-in change inside an anchor cannot hide it. This is the preferred readout.
    haps = []           # (sequence to search for, orientation, forward-orientation haplotype used as the key)
    off = int(cfg.get('haplotype_codon_offset', 0))
    for key in ('wt_haplotype', 'mut_haplotype'):
        h = cfg.get(key, '')
        if h:
            haps.append((h, 'fwd', h))
            haps.append((rc(h), 'rev', h))

    for line in sys.stdin:
        s = line.strip()
        if not s or s[0] in '@+' or not set(s[:4]) <= set('ACGTN'):
            continue
        for seq, orient, hap in haps:
            if seq in s:
                k = ('haplotype', orient, hap)
                counts[k] = counts.get(k, 0) + 1
        i = s.find(L)
        if i >= 0:
            add('left', 'fwd', s[i + len(L):i + len(L) + width])
        i = s.find(Lrc)
        if i >= width:
            add('left', 'rev', rc(s[i - width:i]))
        i = s.find(R)
        if i >= width:
            add('right', 'fwd', s[i - width:i])
        i = s.find(Rrc)
        if i >= 0:
            add('right', 'rev', rc(s[i + len(Rrc):i + len(Rrc) + width]))

    for (anchor, orient, win), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        trip = win[off:off + 3] if anchor == 'haplotype' else (win[:3] if anchor == 'left' else win[-3:])
        print(f'{sample}\t{anchor}\t{orient}\t{trip}\t{n}\t{win}')


if __name__ == '__main__':
    main(*sys.argv[1:])
