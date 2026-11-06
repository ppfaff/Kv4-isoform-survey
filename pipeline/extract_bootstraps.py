#!/usr/bin/env python3
"""Pull the Kv4-complex rows out of one cell's salmon bootstrap replicates.

salmon writes inferential replicates as:
    <quant>/aux_info/bootstrap/names.tsv.gz    one tab-separated line of transcript names
    <quant>/aux_info/bootstrap/bootstraps.gz   raw float64, numBootstraps x numTranscripts, bootstrap-major
Only the Kv4-complex transcripts are kept, so each cell reduces from ~50 MB of raw replicates to a few kB and
the originals can be deleted immediately.

Output (one row per bootstrap replicate per panel variant):
    cell  panel  variant  replicate  reads  frac
where frac is the variant's share of its panel WITHIN that replicate -- the quantity whose across-replicate
spread is the assignment uncertainty we are trying to measure, and which is directly comparable to the
across-cell spread reported by variant_distribution_shape.py.

Usage:  python3 extract_bootstraps.py <cell name> <quant dir> [map.tsv] > out.tsv
The map is a two-column TSV (transcript, vclass). Generate it once on the analysis machine with
    python3 -c "import sys; sys.path.insert(0,'scripts'); from hip_common import final_map; \
                [print(f'{k}\\t{v}') for k,v in final_map.items()]" > tables/transcript_vclass_map.tsv
and copy it alongside the salmon index.
"""
import sys, os, gzip
import numpy as np

PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}
VAR_PANEL = {v: p for p, vs in PANEL_V.items() for v in vs}


def main():
    cell = sys.argv[1]
    qdir = sys.argv[2]
    mapfile = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(qdir), '..', '..',
                                                                 'tables', 'transcript_vclass_map.tsv')
    if not os.path.exists(mapfile):
        mapfile = 'tables/transcript_vclass_map.tsv'
    tx2v = {}
    with open(mapfile) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1] in VAR_PANEL:
                tx2v[parts[0]] = parts[1]
    if not tx2v:
        sys.exit(f'no usable transcript->vclass entries in {mapfile}')

    bdir = os.path.join(qdir, 'aux_info', 'bootstrap')
    with gzip.open(os.path.join(bdir, 'names.tsv.gz'), 'rt') as fh:
        names = fh.read().strip().split('\t')
    raw = gzip.open(os.path.join(bdir, 'bootstraps.gz'), 'rb').read()
    vals = np.frombuffer(raw, dtype=np.float64)
    ntx = len(names)
    if ntx == 0 or vals.size % ntx != 0:
        sys.exit(f'bootstrap size {vals.size} not divisible by {ntx} transcripts')
    nboot = vals.size // ntx
    B = vals.reshape(nboot, ntx)

    # columns of interest, grouped by variant (a variant can have several transcripts)
    cols = {}
    for i, n in enumerate(names):
        # salmon keeps the full GENCODE name; the map may be keyed on the bare transcript id
        v = tx2v.get(n)
        if v is None and '|' in n:
            for part in n.split('|'):
                v = tx2v.get(part)
                if v:
                    break
        if v is None:
            base = n.split('|')[0].split('.')[0]
            v = tx2v.get(base)
        if v:
            cols.setdefault(v, []).append(i)
    if not cols:
        sys.exit('no Kv4-complex transcripts matched the bootstrap names')

    print('cell\tpanel\tvariant\treplicate\treads\tfrac')
    for panel, variants in PANEL_V.items():
        present = [v for v in variants if v in cols]
        if not present:
            continue
        M = np.column_stack([B[:, cols[v]].sum(axis=1) for v in present])   # (nboot, n_variants)
        tot = M.sum(axis=1)
        for bi in range(nboot):
            if tot[bi] <= 0:
                continue
            for vi, v in enumerate(present):
                print(f'{cell}\t{panel}\t{v}\t{bi}\t{M[bi, vi]:.3f}\t{M[bi, vi] / tot[bi]:.6f}')


if __name__ == '__main__':
    main()
