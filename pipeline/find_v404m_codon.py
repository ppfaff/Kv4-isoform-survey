#!/usr/bin/env python3
"""Locate codon 404 of Kcnd2 (the V404M position) in GRCm39 and emit everything the allele counting needs.

Builds the CDS of the canonical Kcnd2 transcript from the GENCODE vM36 GTF plus the genome FASTA (via samtools
faidx), translates it, checks that residue 404 really is valine, and reports:
  * the genomic coordinates of the three bases of codon 404 (and whether the codon is split across an exon junction)
  * the mRNA-sense sequence context, as two 16-nt anchors immediately 5' and 3' of the codon

The two anchors are used independently by count_v404m_allele.sh so that an engineered silent change inside one of
them cannot silently bias the allele ratio: if both anchors give the same answer they are both intact, and if one
recovers far fewer reads than the other, that anchor overlaps an edit in the knock-in allele.

Usage (in the work dir, with ref/ from the dendrite run):
    python3 find_v404m_codon.py [ref/gencode.vM36.primary_assembly.annotation.gtf.gz] \
                                [ref/GRCm39.primary_assembly.genome.fa] [GENE] [RESIDUE]
Writes tables/v404m_codon.tsv (key/value) and prints a summary.
"""
import sys, os, gzip, re, subprocess

GENE = 'Kcnd2'
RESIDUE = 404
EXPECT_AA = 'V'
ANCHOR = 16

CODON_TABLE = {}
_bases = 'TCAG'
_aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, b1 in enumerate(_bases):
    for j, b2 in enumerate(_bases):
        for k, b3 in enumerate(_bases):
            CODON_TABLE[b1 + b2 + b3] = _aas[i * 16 + j * 4 + k]
COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def rc(s):
    return s.translate(COMP)[::-1]


def faidx(fa, chrom, start, end):
    """1-based inclusive region as uppercase sequence."""
    out = subprocess.run(['samtools', 'faidx', fa, f'{chrom}:{start}-{end}'],
                         capture_output=True, text=True, check=True).stdout
    return ''.join(l.strip() for l in out.splitlines()[1:]).upper()


def load_cds(gtf, gene):
    """Return {transcript_name: {'chrom','strand','tags','cds':[(start,end,frame)...]}} for one gene."""
    tx = {}
    op = gzip.open if gtf.endswith('.gz') else open
    with op(gtf, 'rt') as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.rstrip('\n').split('\t')
            if p[2] != 'CDS' or f'gene_name "{gene}"' not in p[8]:
                continue
            m = re.search(r'transcript_name "([^"]+)"', p[8])
            if not m:
                continue
            t = tx.setdefault(m.group(1), {'chrom': p[0], 'strand': p[6], 'tags': set(), 'cds': []})
            t['cds'].append((int(p[3]), int(p[4]), int(p[7]) if p[7] != '.' else 0))
            t['tags'].update(re.findall(r'tag "([^"]+)"', p[8]))
    return tx


def pick(tx):
    """Canonical transcript: Ensembl_canonical / MANE if tagged, else the longest CDS."""
    def key(item):
        name, t = item
        ln = sum(e - s + 1 for s, e, _ in t['cds'])
        return (('MANE_Select' in t['tags']) * 2 + ('Ensembl_canonical' in t['tags']), ln)
    return max(tx.items(), key=key)


def main(gtf='ref/gencode.vM36.primary_assembly.annotation.gtf.gz',
         fa='ref/GRCm39.primary_assembly.genome.fa', gene=GENE, residue=RESIDUE):
    residue = int(residue)
    if not os.path.exists(fa + '.fai'):
        subprocess.run(['samtools', 'faidx', fa], check=True)
    tx = load_cds(gtf, gene)
    if not tx:
        sys.exit(f'{gene} CDS not found in {gtf}')
    name, t = pick(tx)
    strand = t['strand']; chrom = t['chrom']
    exons = sorted(t['cds'], key=lambda e: e[0], reverse=(strand == '-'))
    # CDS sequence in transcript order, and a CDS-index -> genomic-position map
    seq, pos = [], []
    for s, e, _ in exons:
        block = faidx(fa, chrom, s, e)
        if strand == '-':
            block = rc(block)
            coords = list(range(e, s - 1, -1))
        else:
            coords = list(range(s, e + 1))
        seq.append(block); pos.extend(coords)
    cds = ''.join(seq)
    offset = exons[0][2]        # phase of the first CDS block; 0 for a complete transcript
    if offset:
        cds = cds[offset:]; pos = pos[offset:]
    prot = ''.join(CODON_TABLE.get(cds[i:i + 3], 'X') for i in range(0, len(cds) - 2, 3))
    i0 = (residue - 1) * 3
    codon = cds[i0:i0 + 3]
    aa = CODON_TABLE.get(codon, 'X')
    gpos = pos[i0:i0 + 3]
    split = (max(gpos) - min(gpos)) != 2
    left = cds[max(0, i0 - ANCHOR):i0]
    right = cds[i0 + 3:i0 + 3 + ANCHOR]

    print(f'{gene}: transcript {name} ({chrom}, strand {strand}), CDS {len(cds)} nt, protein {len(prot)} aa', file=sys.stderr)
    print(f'  protein starts {prot[:12]}... ; residue {residue} = {aa} (codon {codon})', file=sys.stderr)
    print(f'  genomic positions of codon {residue}: {chrom}:{gpos[0]},{gpos[1]},{gpos[2]}'
          + ('   *** SPLIT ACROSS AN EXON JUNCTION ***' if split else ''), file=sys.stderr)
    print(f'  mRNA context: {left} [{codon}] {right}', file=sys.stderr)
    if aa != EXPECT_AA:
        print(f'  WARNING: residue {residue} is {aa}, not the expected {EXPECT_AA} — check the residue numbering '
              f'(this script numbers from the initiator Met of {name})', file=sys.stderr)

    os.makedirs('tables', exist_ok=True)
    with open('tables/v404m_codon.tsv', 'w') as f:
        for k, v in [('gene', gene), ('transcript', name), ('chrom', chrom), ('strand', strand),
                     ('residue', residue), ('aa', aa), ('wt_codon', codon),
                     ('gpos1', gpos[0]), ('gpos2', gpos[1]), ('gpos3', gpos[2]),
                     ('codon_split_across_exons', int(split)),
                     ('left_anchor', left), ('right_anchor', right),
                     ('left_anchor_rc', rc(left)), ('right_anchor_rc', rc(right)),
                     ('mrna_window', left + codon + right)]:
            f.write(f'{k}\t{v}\n')
    with open('tables/v404m_codon.pos', 'w') as f:
        for g in sorted(gpos):
            f.write(f'{chrom}\t{g}\n')
    print('\nwrote tables/v404m_codon.tsv and tables/v404m_codon.pos', file=sys.stderr)


if __name__ == '__main__':
    main(*sys.argv[1:])
