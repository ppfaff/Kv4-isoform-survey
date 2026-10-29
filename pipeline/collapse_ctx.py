#!/usr/bin/env python3
"""Collapse GSE71585 (Tasic 2016 cortex) salmon quant.sf files across machines.

Per machine (in the kv4_isoforms / kv4_cortex directory, needs ref/gencode.vM36.transcripts.fa.gz):
    python3 collapse_ctx.py tables/ctx_percell_$(hostname).csv

Merge on one machine (needs tables/ena_runs.tsv and anno.csv there):
    python3 collapse_ctx.py --merge tables/ctx_percell_*.csv
writes tables/kv4_isoform_tpm_per_cell.csv, tables/kv4_isoform_per_cell_annotated.csv,
tables/kv4_isoform_by_type.csv, tables/run_qc.csv — identical layout to the v4 script's output.
"""
import sys, glob, gzip, json, os, re
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


def collapse_local(out_csv):
    tx = load_tx_map()
    rows, qc = [], []
    quants = sorted(glob.glob('quant/*/quant.sf'))
    for q in quants:
        run = q.split('/')[1]
        d = pd.read_csv(q, sep='\t', index_col=0)
        try:
            m = json.load(open(f'quant/{run}/aux_info/meta_info.json'))
            qc.append((run, m.get('num_processed'), m.get('percent_mapped')))
        except Exception:
            qc.append((run, None, None))
        for t in d.index.intersection(list(tx)):
            g, name, btype = tx[t]
            rows.append((run, g, name, t, btype, d.at[t, 'TPM'], d.at[t, 'NumReads'], d.at[t, 'EffectiveLength']))
    T = pd.DataFrame(rows, columns=['run', 'gene', 'transcript', 'ensmust', 'biotype', 'TPM', 'reads', 'eff_len'])
    T = T.merge(pd.DataFrame(qc, columns=['run', 'num_processed', 'percent_mapped']), on='run', how='left')
    os.makedirs(os.path.dirname(out_csv) or '.', exist_ok=True)
    T.to_csv(out_csv, index=False)
    print(f'{len(quants)} runs collapsed -> {out_csv} ({len(T)} rows)')


def merge(csvs):
    T = pd.concat([pd.read_csv(c) for c in csvs], ignore_index=True)
    dup = T.duplicated(subset=['run', 'ensmust']).sum()
    if dup:
        print(f'WARNING: {dup} duplicate (run, transcript) rows across machines — keeping first')
        T = T.drop_duplicates(subset=['run', 'ensmust'])
    T[['run', 'num_processed', 'percent_mapped']].drop_duplicates().to_csv('tables/run_qc.csv', index=False)

    ena = pd.read_csv('tables/ena_runs.tsv', sep='\t').rename(columns={'run_accession': 'run'})
    T = T.merge(ena[['run', 'sample_alias', 'sample_title', 'sample_accession']], on='run', how='left')
    pat = re.compile(r'^[A-Za-z0-9\-]+_tdT(positive|negative)_cell_\d+$')

    def cellname(r):
        for c in ['sample_title', 'sample_alias']:
            v = str(r[c])
            if pat.match(v):
                return v
        return str(r['sample_title'])
    T['sample_name'] = T.apply(cellname, axis=1)
    T.to_csv('tables/kv4_isoform_tpm_per_cell.csv', index=False)
    print(f'runs in ENA project: {ena.run.nunique()}   quantified: {T.run.nunique()}   '
          f'missing: {ena.run.nunique() - T.run.nunique()}')

    if os.path.exists('anno.csv'):
        a = pd.read_csv('anno.csv')
        a = a[(a.core_intermediate != 'unclassified') & (a.cre_driver_1 != 'CAV2-Cre')]
        J = T.merge(a[['sample_name', 'primary_type_id', 'primary_type_label', 'broad_type', 'dissection']], on='sample_name')
        print('matched annotated cells:', J.sample_name.nunique(), 'of', a.sample_name.nunique())
        J['gene_tpm'] = J.groupby(['sample_name', 'gene']).TPM.transform('sum')
        J['gene_reads'] = J.groupby(['sample_name', 'gene']).reads.transform('sum')
        J['frac'] = J.TPM / J.gene_tpm.where(J.gene_tpm > 0)
        J.to_csv('tables/kv4_isoform_per_cell_annotated.csv', index=False)
        out = J.groupby(['gene', 'transcript', 'ensmust', 'primary_type_id', 'primary_type_label']).agg(
            mean_TPM=('TPM', 'mean'), mean_frac=('frac', 'mean'),
            pct_cells_detected=('reads', lambda x: (x >= 3).mean() * 100),
            n=('sample_name', 'nunique')).reset_index()
        out.to_csv('tables/kv4_isoform_by_type.csv', index=False)
        missing = sorted(set(a.sample_name) - set(J.sample_name))
        if missing:
            pd.Series(missing).to_csv('tables/missing_cells.txt', index=False, header=False)
            print(f'  {len(missing)} annotated cells not quantified -> tables/missing_cells.txt')
    print('wrote tables/')


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0] == '--merge':
        merge(args[1:])
    elif len(args) == 1:
        collapse_local(args[0])
    else:
        print(__doc__); sys.exit(1)
