#!/usr/bin/env python3
"""Build exon / intron annotation for the 9 Kv4-complex genes from the GENCODE vM36 GTF.

Writes (into ./ref):
  kv4_exons.saf     merged exon blocks per gene (union over all transcripts)         -> featureCounts -F SAF
  kv4_introns.saf   gaps between merged exon blocks, one row per intron, named        -> featureCounts -F SAF
                    <gene>_intron_<k> with k counted 5'->3' in the gene's orientation
  kv4_introns.bed   same intervals as BED (for samtools bedcov per-intron depth)
  kv4_genes.bed     gene spans +/- 2 kb (for extracting the small gene-region BAMs)

Usage: python3 make_intron_saf.py ref/gencode.vM36.primary_assembly.annotation.gtf.gz
"""
import sys, gzip, re, os
from collections import defaultdict

GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10']
PAD = 2000


def merge(iv):
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def main(gtf):
    exons = defaultdict(list); span = {}; strand = {}; chrom = {}
    op = gzip.open if gtf.endswith('.gz') else open
    with op(gtf, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            if p[2] not in ('gene', 'exon'):
                continue
            m = re.search(r'gene_name "([^"]+)"', p[8])
            if not m or m.group(1) not in GENES:
                continue
            g = m.group(1); s, e = int(p[3]), int(p[4])
            if p[2] == 'gene':
                span[g] = (s, e); strand[g] = p[6]; chrom[g] = p[0]
            else:
                exons[g].append((s, e))
    os.makedirs('ref', exist_ok=True)
    with open('ref/kv4_exons.saf', 'w') as fe, open('ref/kv4_introns.saf', 'w') as fi, \
         open('ref/kv4_introns.bed', 'w') as fb, open('ref/kv4_genes.bed', 'w') as fg:
        fe.write('GeneID\tChr\tStart\tEnd\tStrand\n'); fi.write('GeneID\tChr\tStart\tEnd\tStrand\n')
        for g in GENES:
            if g not in span:
                print(f'WARNING: {g} not found in GTF', file=sys.stderr); continue
            c, st = chrom[g], strand[g]
            blocks = merge(exons[g])
            for s, e in blocks:
                fe.write(f'{g}\t{c}\t{s}\t{e}\t{st}\n')
            gaps = [(blocks[i][1] + 1, blocks[i + 1][0] - 1) for i in range(len(blocks) - 1) if blocks[i + 1][0] - blocks[i][1] > 1]
            order = gaps if st == '+' else gaps[::-1]
            for k, (s, e) in enumerate(order, 1):
                fi.write(f'{g}_intron_{k}\t{c}\t{s}\t{e}\t{st}\n')
                fb.write(f'{c}\t{s - 1}\t{e}\t{g}_intron_{k}\t0\t{st}\n')
            gs, ge = span[g]
            fg.write(f'{c}\t{max(0, gs - PAD - 1)}\t{ge + PAD}\t{g}\t0\t{st}\n')
            print(f'{g}: {len(blocks)} merged exon blocks, {len(gaps)} introns, span {ge - gs + 1:,} bp ({st})', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1])
