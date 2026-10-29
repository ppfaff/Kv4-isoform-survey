#!/usr/bin/env python3
"""Pull a handful of marker genes per cell out of the Allen SSv4 (mouse ctx+hpf, Yao 2021) expression_matrix.hdf5
(the released RSEM gene-level exon counts), for joining to our isoform table by sample_name.

Setup (Lambda, kv4 env):
    micromamba install -y -n kv4 -c conda-forge h5py          # if 'import h5py' fails
    mkdir -p ~/allen_meta/allen_matrix && cd ~/allen_meta/allen_matrix
    curl -fL -o expression_matrix.hdf5 https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_smart-seq/expression_matrix.hdf5
    python3 ../extract_markers.py expression_matrix.hdf5 ../tables/allen_markers.csv

Output: one row per cell (all 73k cells in the release), columns = raw exon counts for each marker gene,
plus total_exon_counts for CPM normalisation. Our 8,131 hippocampal cells are a subset; the join happens later.
"""
import sys
import numpy as np, pandas as pd, h5py

GENES = [
    # immature / adult-born granule-cell markers
    'Dcx', 'Calb2', 'Sox11', 'Sox4', 'Neurod1', 'Igfbpl1', 'Tubb3',
    # mature granule-cell / DG identity markers
    'Calb1', 'Prox1', 'C1ql2', 'Pdzd2', 'Dsp',
    # cross-check against our salmon quantification
    'Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10',
    # housekeeping sanity
    'Actb', 'Gapdh',
]

def main(h5path, out_csv):
    with h5py.File(h5path, 'r') as f:
        genes = np.array([g.decode() if isinstance(g, bytes) else str(g) for g in f['gene_names'][:]])
        samples = np.array([s.decode() if isinstance(s, bytes) else str(s) for s in f['sample_names'][:]])
        n_genes, n_samples = len(genes), len(samples)
        # pick whichever sparse layout is indexed by gene (p has n_genes+1 entries), so one gene = one slice
        grp = None
        for name in ['data/exon', 'data/t_exon']:
            if name in f and f[name]['p'].shape[0] == n_genes + 1:
                grp = f[name]; break
        if grp is None:
            raise SystemExit('could not find a gene-indexed sparse group; h5ls the file and tell me the layout')
        p = grp['p'][:]
        gidx = {g: i for i, g in enumerate(genes)}
        cols = {}
        for g in GENES:
            if g not in gidx:
                print('not in matrix:', g); continue
            k = gidx[g]
            a, b = int(p[k]), int(p[k + 1])
            rows = grp['i'][a:b].astype(np.int64); vals = grp['x'][a:b]   # R writes the index arrays as doubles
            col = np.zeros(n_samples, dtype=np.float64); col[rows] = vals
            cols[g] = col
        tot = f['data/total_exon_counts'][:] if 'data/total_exon_counts' in f else None
    M = pd.DataFrame(cols, index=pd.Index(samples, name='sample_name'))
    if tot is not None:
        M['total_exon_counts'] = tot
    M.to_csv(out_csv)
    print(f'wrote {out_csv}: {M.shape[0]} cells x {M.shape[1]} columns')
    print(M.describe().T[['mean', '50%', 'max']].round(1).to_string())

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
