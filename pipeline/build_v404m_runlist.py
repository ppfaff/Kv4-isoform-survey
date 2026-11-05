#!/usr/bin/env python3
"""Build the sample sheet for the Kv4.2-V404M vs wild-type bulk RNA-seq (3-month littermates, 4 tissues).

Input: the core's sample-info CSV (ID, Tissue, Genotype, Sex) plus the directory holding the FASTQs, named
<ID>_R1_001.fastq.gz / <ID>_R2_001.fastq.gz (the I1/I2 index reads are not used). Checks every file exists and
reports its size, and orders samples by tissue priority so the most informative tissue finishes first.

Usage (WSL on Lambda, after mounting T:):
    python3 build_v404m_runlist.py "/mnt/t/Laurence David/RNA_Seq_VM_3mo/RNA-seq_KV404M3Mo_Sample_Info.csv" \
                                   "/mnt/t/Laurence David/RNA_Seq_VM_3mo"
Writes tables/v404m_samples.tsv:
    sample, tissue, genotype, sex, r1, r2, r1_gb, r2_gb
"""
import sys, os
import pandas as pd

TISSUE_ORDER = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']


def main(csv_path, fastq_dir):
    S = pd.read_csv(csv_path).dropna(how='all')
    S.columns = [c.strip().lower() for c in S.columns]
    S = S.rename(columns={'id': 'sample'})
    for c in ['sample', 'tissue', 'genotype', 'sex']:
        S[c] = S[c].astype(str).str.strip()
    S['r1'] = S['sample'].map(lambda s: os.path.join(fastq_dir, f'{s}_R1_001.fastq.gz'))
    S['r2'] = S['sample'].map(lambda s: os.path.join(fastq_dir, f'{s}_R2_001.fastq.gz'))
    for side in ['r1', 'r2']:
        S[side + '_gb'] = S[side].map(lambda p: round(os.path.getsize(p) / 2**30, 2) if os.path.exists(p) else 0.0)
    missing = S[(S.r1_gb == 0) | (S.r2_gb == 0)]
    S['tissue_rank'] = S.tissue.map({t: i for i, t in enumerate(TISSUE_ORDER)}).fillna(99)
    S = S.sort_values(['tissue_rank', 'genotype', 'sample']).drop(columns='tissue_rank')
    os.makedirs('tables', exist_ok=True)
    S.to_csv('tables/v404m_samples.tsv', sep='\t', index=False)

    print(S.groupby(['tissue', 'genotype', 'sex']).size().unstack(fill_value=0).to_string(), file=sys.stderr)
    print(f'\n{len(S)} samples, {S.r1_gb.sum() + S.r2_gb.sum():.0f} GB of R1+R2 -> tables/v404m_samples.tsv', file=sys.stderr)
    if len(missing):
        print(f'\nWARNING: {len(missing)} samples missing a FASTQ:', file=sys.stderr)
        print(missing[['sample', 'tissue', 'genotype', 'r1_gb', 'r2_gb']].to_string(index=False), file=sys.stderr)
    # per-tissue genotype balance, which the analysis has to respect
    b = S.groupby(['tissue', 'genotype']).size().unstack(fill_value=0)
    thin = b[(b.min(axis=1) < 4)]
    if len(thin):
        print('\nNOTE: these tissues have fewer than 4 samples in one genotype; treat their statistics as weak:', file=sys.stderr)
        print(thin.to_string(), file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
