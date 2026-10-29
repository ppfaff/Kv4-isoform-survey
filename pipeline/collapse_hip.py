#!/usr/bin/env python3
"""Collapse salmon quant.sf files from ALL machines into the hippocampal isoform tables.

The pipeline's built-in collapse step only sees the quant/ directory on the machine it ran on.
After the 3-way split, Alienware and PFAFF10 each hold ~900-1000 cells' worth of quant/ output
that Lambda never saw. Rather than copying ~15 GB of quant directories around, each machine
runs this script on its own quant/ to produce one small per-cell CSV (one row per transcript
per cell, 9 genes only -> a few MB), and those CSVs are concatenated on Lambda for the final
subclass-level table.

Usage, on EACH machine, in ~/allen_meta:
    python3 collapse_hip.py hip_filtered_runlist.tsv tables/hip_percell_$(hostname).csv

Then on Lambda, after copying the other two machines' CSVs into tables/:
    python3 collapse_hip.py --merge hip_filtered_runlist.tsv tables/hip_percell_*.csv
which writes tables/hip_isoform_per_cell_annotated.csv, tables/hip_isoform_by_subclass.csv,
tables/run_qc.csv and prints a completeness summary against the master runlist.
"""
import sys, glob, gzip, json, os
import pandas as pd

GENES = {'Kcnd1', 'Kcnd2', 'Kcnd3', 'Dpp6', 'Dpp10', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4'}


def load_tx_map(fa='ref/gencode.vM36.transcripts.fa.gz'):
    tx = {}
    with gzip.open(fa, 'rt') as f:
        for l in f:
            if l[0] == '>':
                p = l[1:].strip().split('|')
                if p[5] in GENES:
                    tx[p[0]] = (p[5], p[4], p[7] if len(p) > 7 else '')
    return tx


def collapse_local(runlist_path, out_csv):
    tx = load_tx_map()
    rows, qc = [], []
    quants = sorted(glob.glob('quant/*/quant.sf'))
    for q in quants:
        sample = q.split('/')[1]
        d = pd.read_csv(q, sep='\t', index_col=0)
        try:
            m = json.load(open(f'quant/{sample}/aux_info/meta_info.json'))
            qc.append((sample, m.get('num_processed'), m.get('percent_mapped')))
        except Exception:
            qc.append((sample, None, None))
        hit = d.index.intersection(list(tx))
        for t in hit:
            g, name, btype = tx[t]
            rows.append((sample, g, name, t, btype, d.at[t, 'TPM'], d.at[t, 'NumReads'], d.at[t, 'EffectiveLength']))
    T = pd.DataFrame(rows, columns=['sample_name', 'gene', 'transcript', 'ensmust', 'biotype', 'TPM', 'reads', 'eff_len'])
    Q = pd.DataFrame(qc, columns=['sample_name', 'num_processed', 'percent_mapped'])
    T = T.merge(Q, on='sample_name', how='left')
    os.makedirs(os.path.dirname(out_csv) or '.', exist_ok=True)
    T.to_csv(out_csv, index=False)
    print(f'{len(quants)} cells collapsed -> {out_csv} ({len(T)} rows)')


def merge(runlist_path, csvs):
    parts = [pd.read_csv(c) for c in csvs]
    T = pd.concat(parts, ignore_index=True)
    dup = T.duplicated(subset=['sample_name', 'ensmust']).sum()
    if dup:
        print(f'WARNING: {dup} duplicate (cell, transcript) rows across machines — keeping first')
        T = T.drop_duplicates(subset=['sample_name', 'ensmust'])
    runlist = pd.read_csv(runlist_path, sep='\t')
    J = T.merge(runlist[['sample_name', 'subclass_label', 'cluster_label', 'region_label']], on='sample_name', how='left')
    J['gene_tpm'] = J.groupby(['sample_name', 'gene']).TPM.transform('sum')
    J['gene_reads'] = J.groupby(['sample_name', 'gene']).reads.transform('sum')
    J['frac'] = J.TPM / J.gene_tpm.where(J.gene_tpm > 0)
    os.makedirs('tables', exist_ok=True)
    J.to_csv('tables/hip_isoform_per_cell_annotated.csv', index=False)
    J[['sample_name', 'num_processed', 'percent_mapped']].drop_duplicates().to_csv('tables/run_qc.csv', index=False)
    out = J.groupby(['gene', 'transcript', 'ensmust', 'subclass_label']).agg(
        mean_TPM=('TPM', 'mean'), mean_frac=('frac', 'mean'),
        pct_cells_detected=('reads', lambda x: (x >= 3).mean() * 100),
        n=('sample_name', 'nunique')).reset_index()
    out.to_csv('tables/hip_isoform_by_subclass.csv', index=False)

    done = set(J.sample_name)
    want = set(runlist.sample_name)
    missing = sorted(want - done)
    print(f'cells in runlist: {len(want)}   quantified: {len(done & want)}   missing: {len(missing)}')
    if missing:
        pd.Series(missing, name='sample_name').to_csv('tables/missing_cells.txt', index=False, header=False)
        print('  -> tables/missing_cells.txt')
    print('per-subclass cell counts:')
    print(J.drop_duplicates('sample_name').subclass_label.value_counts().to_string())
    print('wrote tables/hip_isoform_per_cell_annotated.csv, hip_isoform_by_subclass.csv, run_qc.csv')


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--merge':
        merge(args[1], args[2:])
    elif len(args) == 2:
        collapse_local(*args)
    else:
        print(__doc__); sys.exit(1)
