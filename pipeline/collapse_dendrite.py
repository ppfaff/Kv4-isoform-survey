#!/usr/bin/env python3
"""Collapse the dendrite/soma salmon runs into per-sample tables.

Reads quant/<run>/quant.sf for every run in tables/dendrite_samples.tsv, keeps the 9 Kv4-complex genes at
transcript level, sums runs that belong to the same library (same dataset/unit/region/compartment —
Middleton samples were sometimes split over lanes), and writes:

  tables/dendrite_isoform_per_sample.csv   one row per transcript per sample: dataset, sample_id, region,
                                           compartment, unit, gene, transcript, ensmust, biotype, vclass,
                                           TPM, reads, gene_tpm, gene_reads, frac (of gene reads)
  tables/dendrite_gene_tpm.csv.gz          gene-level TPM for ALL genes x samples (for marker/QC checks and
                                           any global normalisation the analysis wants)
  tables/dendrite_run_qc.csv               num_processed, percent_mapped, inferred library type per run

Run in the work directory on Lambda after run_dendrite_isoforms.sh:   python3 collapse_dendrite.py
"""
import sys, os, glob, gzip, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scripts'))
from hip_common import final_map, PANEL_OF  # noqa: E402

GENES = ['Kcnd1', 'Kcnd2', 'Kcnd3', 'Dpp6', 'Dpp10', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']


def load_tx_map(fa='ref/gencode.vM36.transcripts.fa.gz'):
    """ENSMUST -> (gene_name, transcript_name, biotype) for all transcripts; gene map for all genes."""
    tx, gene_of = {}, {}
    with gzip.open(fa, 'rt') as f:
        for l in f:
            if l[0] == '>':
                p = l[1:].strip().split('|')
                gene_of[p[0]] = p[5]
                if p[5] in GENES:
                    tx[p[0]] = (p[5], p[4], p[7] if len(p) > 7 else '')
    return tx, gene_of


def main(samples='tables/dendrite_samples.tsv'):
    S = pd.read_csv(samples, sep='\t', dtype=str).fillna('')
    bad = S[(S.region == '?') | (S.compartment == '?') | (S.unit == '?')]
    if len(bad):
        sys.exit(f'{len(bad)} runs still have "?" in region/compartment/unit in {samples}; fix those first:\n'
                 + bad[['dataset', 'run_accession', 'title']].to_string(index=False))
    S['sample_id'] = S.dataset + '_' + S.region + '_' + S.compartment + '_' + S.unit
    tx, gene_of = load_tx_map()
    rows, qc, gtpm = [], [], {}
    missing = []
    for _, r in S.iterrows():
        q = f'quant/{r.run_accession}/quant.sf'
        if not os.path.exists(q):
            missing.append(r.run_accession); continue
        d = pd.read_csv(q, sep='\t', index_col=0)
        try:
            m = json.load(open(f'quant/{r.run_accession}/aux_info/meta_info.json'))
            qc.append((r.run_accession, r.sample_id, m.get('num_processed'), m.get('percent_mapped'), ','.join(m.get('library_types', []))))
        except Exception:
            qc.append((r.run_accession, r.sample_id, None, None, ''))
        # gene-level TPM for everything (summed over runs of the same sample later)
        g = d.TPM.groupby(d.index.map(gene_of)).sum()
        gtpm.setdefault(r.sample_id, []).append(g)
        for t in d.index.intersection(list(tx)):
            gname, tname, btype = tx[t]
            rows.append((r.dataset, r.sample_id, r.region, r.compartment, r.unit, r.run_accession,
                         gname, tname, t, btype, d.at[t, 'TPM'], d.at[t, 'NumReads']))
    if missing:
        print(f'WARNING: {len(missing)} runs have no quant.sf yet: {" ".join(missing)}', file=sys.stderr)
    T = pd.DataFrame(rows, columns=['dataset', 'sample_id', 'region', 'compartment', 'unit', 'run', 'gene', 'transcript',
                                    'ensmust', 'biotype', 'TPM', 'reads'])
    Q = pd.DataFrame(qc, columns=['run', 'sample_id', 'num_processed', 'percent_mapped', 'library_types'])
    # merge runs of the same sample: reads add; TPM is averaged weighted by processed reads
    T = T.merge(Q[['run', 'num_processed']], on='run', how='left')
    T['num_processed'] = T.num_processed.fillna(1).astype(float)
    T['w'] = T.TPM * T.num_processed
    A = T.groupby(['dataset', 'sample_id', 'region', 'compartment', 'unit', 'gene', 'transcript', 'ensmust', 'biotype'], as_index=False).agg(
        reads=('reads', 'sum'), w=('w', 'sum'), n=('num_processed', 'sum'))
    A['TPM'] = A.w / A.n
    A = A.drop(columns=['w', 'n'])
    A['vclass'] = A.transcript.map(final_map)
    A['panel'] = A.gene.map(PANEL_OF)
    A['gene_tpm'] = A.groupby(['sample_id', 'gene']).TPM.transform('sum')
    A['gene_reads'] = A.groupby(['sample_id', 'gene']).reads.transform('sum')
    A['frac'] = A.reads / A.gene_reads.where(A.gene_reads > 0)
    os.makedirs('tables', exist_ok=True)
    A.to_csv('tables/dendrite_isoform_per_sample.csv', index=False)
    Q.to_csv('tables/dendrite_run_qc.csv', index=False)
    G = pd.DataFrame({sid: pd.concat(v, axis=1).mean(axis=1) for sid, v in gtpm.items()})
    G.index.name = 'gene'
    G.to_csv('tables/dendrite_gene_tpm.csv.gz', compression='gzip')
    print(f'{A.sample_id.nunique()} samples, {len(A)} transcript rows -> tables/dendrite_isoform_per_sample.csv')
    print(A.groupby(['dataset', 'region', 'compartment']).sample_id.nunique().to_string())
    print('\nlibrary types inferred by salmon:'); print(Q.library_types.value_counts().to_string())


if __name__ == '__main__':
    main(*sys.argv[1:])
