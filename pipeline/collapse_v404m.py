#!/usr/bin/env python3
"""Collapse the V404M salmon runs into the matrices the analysis needs.

Reads quant/<sample>/quant.sf for every sample in tables/v404m_samples.tsv and writes:
  tables/v404m_gene_counts.csv.gz      gene x sample estimated read counts (all genes) — differential expression
  tables/v404m_gene_tpm.csv.gz         gene x sample TPM (all genes)
  tables/v404m_isoform_per_sample.csv  transcript-level rows for the 9 Kv4-complex genes, with variant class
  tables/v404m_gene_id_map.csv         ENSMUSG -> symbol, so the STAR count matrix can be used interchangeably
  tables/v404m_salmon_qc.csv           reads processed, percent mapped, inferred library type per sample

Run in the work dir after run_v404m_salmon.sh:   python3 collapse_v404m.py
"""
import sys, os, glob, gzip, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
for cand in [os.path.join(HERE, '..', 'scripts'), HERE, os.getcwd(),
             os.path.expanduser('~/Kv4-isoform-survey/scripts'), os.path.expanduser('~/kv4-isoform-survey/scripts')]:
    if os.path.exists(os.path.join(cand, 'hip_common.py')):
        sys.path.insert(0, cand); break
from hip_common import final_map, PANEL_OF  # noqa: E402

GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Dpp6', 'Dpp10', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']


def load_tx_map(fa='ref/gencode.vM36.transcripts.fa.gz'):
    """ENSMUST -> (gene_symbol, transcript_name, biotype, ENSMUSG) for every transcript."""
    tx = {}
    with gzip.open(fa, 'rt') as f:
        for l in f:
            if l[0] == '>':
                p = l[1:].strip().split('|')
                tx[p[0]] = (p[5], p[4], p[7] if len(p) > 7 else '', p[1])
    return tx


def main(samples='tables/v404m_samples.tsv'):
    S = pd.read_csv(samples, sep='\t')
    tx = load_tx_map()
    gene_of = {t: v[0] for t, v in tx.items()}
    counts, tpms, rows, qc, missing = {}, {}, [], [], []
    for _, r in S.iterrows():
        q = f'quant/{r["sample"]}/quant.sf'
        if not os.path.exists(q):
            missing.append(r['sample']); continue
        d = pd.read_csv(q, sep='\t', index_col=0)
        g = d.index.map(gene_of)
        counts[r['sample']] = d.NumReads.groupby(g).sum()
        tpms[r['sample']] = d.TPM.groupby(g).sum()
        try:
            m = json.load(open(f'quant/{r["sample"]}/aux_info/meta_info.json'))
            qc.append((r['sample'], r['tissue'], r['genotype'], r['sex'], m.get('num_processed'),
                       m.get('percent_mapped'), ','.join(m.get('library_types', []))))
        except Exception:
            qc.append((r['sample'], r['tissue'], r['genotype'], r['sex'], None, None, ''))
        for t in d.index:
            info = tx.get(t)
            if info and info[0] in GENES:
                rows.append((r['sample'], r['tissue'], r['genotype'], r['sex'], info[0], info[1], t, info[2],
                             d.at[t, 'TPM'], d.at[t, 'NumReads']))
    if missing:
        print(f'WARNING: {len(missing)} samples have no quant.sf yet: {" ".join(missing)}', file=sys.stderr)
    os.makedirs('tables', exist_ok=True)

    C = pd.DataFrame(counts).fillna(0); C.index.name = 'gene'
    T = pd.DataFrame(tpms).fillna(0); T.index.name = 'gene'
    C.to_csv('tables/v404m_gene_counts.csv.gz', compression='gzip')
    T.to_csv('tables/v404m_gene_tpm.csv.gz', compression='gzip')
    pd.DataFrame(qc, columns=['sample', 'tissue', 'genotype', 'sex', 'num_processed', 'percent_mapped',
                              'library_types']).to_csv('tables/v404m_salmon_qc.csv', index=False)
    gm = pd.DataFrame(sorted({(v[3].split('.')[0], v[0]) for v in tx.values()}), columns=['gene_id', 'symbol'])
    gm.to_csv('tables/v404m_gene_id_map.csv', index=False)

    I = pd.DataFrame(rows, columns=['sample', 'tissue', 'genotype', 'sex', 'gene', 'transcript', 'ensmust',
                                    'biotype', 'TPM', 'reads'])
    I['vclass'] = I.transcript.map(final_map)
    I['panel'] = I.gene.map(PANEL_OF)
    I['gene_reads'] = I.groupby(['sample', 'gene']).reads.transform('sum')
    I['gene_tpm'] = I.groupby(['sample', 'gene']).TPM.transform('sum')
    I['frac'] = I.reads / I.gene_reads.where(I.gene_reads > 0)
    I.to_csv('tables/v404m_isoform_per_sample.csv', index=False)

    print(f'{C.shape[1]} samples, {C.shape[0]} genes -> tables/v404m_gene_counts.csv.gz')
    print(pd.DataFrame(qc, columns=['sample', 'tissue', 'genotype', 'sex', 'num_processed', 'percent_mapped',
                                    'library_types']).groupby(['tissue', 'genotype']).size().to_string())
    q = pd.read_csv('tables/v404m_salmon_qc.csv')
    print('\nlibrary types:', q.library_types.value_counts().to_dict())
    print('percent mapped: median %.1f, range %.1f-%.1f' % (q.percent_mapped.median(), q.percent_mapped.min(), q.percent_mapped.max()))
    print('reads processed: median %.1fM' % (q.num_processed.median() / 1e6))


if __name__ == '__main__':
    main(*sys.argv[1:])
