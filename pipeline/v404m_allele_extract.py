#!/usr/bin/env python3
"""Tabulate the trinucleotide at a codon from raw reads, given flanking anchors (no alignment).

Reads anchor definitions from tables/v404m_codon.tsv and FASTQ sequence lines from stdin (normally pre-filtered
with `grep -F`, so only candidate reads arrive here). For every read it reports the three bases sitting at the
codon position, resolved four ways:

    left/fwd   read contains LEFT_ANCHOR            -> the 3 bases after it
    left/rev   read contains rc(LEFT_ANCHOR)        -> the 3 bases before it, reverse-complemented
    right/fwd  read contains RIGHT_ANCHOR           -> the 3 bases before it
    right/rev  read contains rc(RIGHT_ANCHOR)       -> the 3 bases after it, reverse-complemented

The left and right anchors are independent: agreement between them means neither anchor overlaps an engineered
change in the knock-in allele. Every observed triplet is reported, not just the two expected ones, so an
unexpected mutant sequence (extra silent edits) shows up rather than being silently dropped.

Usage:  zcat R1.fastq.gz R2.fastq.gz | grep -F -e ... | python3 v404m_allele_extract.py SAMPLE [codon.tsv]
Prints TSV rows: sample, anchor, orientation, triplet, count
"""
import sys

COMP = str.maketrans('ACGTN', 'TGCAN')


def rc(s):
    return s.translate(COMP)[::-1]


def main(sample, tsv='tables/v404m_codon.tsv'):
    cfg = dict(l.rstrip('\n').split('\t', 1) for l in open(tsv) if '\t' in l)
    L, R = cfg['left_anchor'], cfg['right_anchor']
    Lrc, Rrc = rc(L), rc(R)
    counts = {}

    def add(anchor, orient, trip):
        if len(trip) == 3 and set(trip) <= set('ACGT'):
            counts[(anchor, orient, trip)] = counts.get((anchor, orient, trip), 0) + 1

    for line in sys.stdin:
        s = line.strip()
        if not s or s[0] in '@+' or not set(s[:4]) <= set('ACGTN'):
            continue
        i = s.find(L)
        if i >= 0:
            add('left', 'fwd', s[i + len(L):i + len(L) + 3])
        i = s.find(Lrc)
        if i >= 3:
            add('left', 'rev', rc(s[i - 3:i]))
        i = s.find(R)
        if i >= 3:
            add('right', 'fwd', s[i - 3:i])
        i = s.find(Rrc)
        if i >= 0:
            add('right', 'rev', rc(s[i + len(Rrc):i + len(Rrc) + 3]))

    for (anchor, orient, trip), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f'{sample}\t{anchor}\t{orient}\t{trip}\t{n}')


if __name__ == '__main__':
    main(*sys.argv[1:])
