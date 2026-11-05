#!/usr/bin/env python3
"""Build the sample sheet for the Kv4.2-V404M vs wild-type bulk RNA-seq (3-month littermates, 4 tissues).

Input: the core's sample-info CSV (ID, Tissue, Genotype, Sex) plus the directory holding the FASTQs (the I1/I2
index reads are not used). Checks every file exists and reports its size, and orders samples by tissue priority so
the most informative tissue finishes first.

Files are matched on the animal number (PPF_<n>), not on the full ID, because a library that was under-sequenced
in the first run may have been topped up on later lanes under a different S-number: PPF_75704 exists as
S69_L006 (0.4 GB), S8_L002 (1.0 GB) and S52_L008 (0.7 GB). All parts of a sample are listed, ';'-separated, and
the pipeline concatenates them when staging (concatenated gzip members are a valid gzip stream). Samples made of
more than one part, and samples whose total size is an outlier, are both reported so they can be checked.

Usage (WSL on Lambda, after mounting T:):
    python3 build_v404m_runlist.py "/mnt/t/Laurence David/RNA_Seq_VM_3mo/RNA-seq_KV404M3Mo_Sample_Info.csv" \
                                   "/mnt/t/Laurence David/RNA_Seq_VM_3mo"
Writes tables/v404m_samples.tsv:
    sample, tissue, genotype, sex, r1, r2, r1_gb, r2_gb
"""
import sys, os, re, glob
import pandas as pd

TISSUE_ORDER = ['hippocampus', 'cortex', 'thalamus', 'cerebellum']


def parts_for(fastq_dir, sample_id, mate):
    """Every FASTQ belonging to this animal for one mate, sorted, as a ';'-joined string."""
    m = re.match(r'(PPF_\d+)', sample_id)
    prefix = m.group(1) if m else sample_id
    hits = sorted(glob.glob(os.path.join(fastq_dir, f'{prefix}_S*_R{mate}_001.fastq.gz')))
    if not hits:   # fall back to the exact ID, in case the naming differs
        p = os.path.join(fastq_dir, f'{sample_id}_R{mate}_001.fastq.gz')
        hits = [p] if os.path.exists(p) else []
    return ';'.join(hits)


def gb(paths):
    return round(sum(os.path.getsize(p) for p in paths.split(';') if p and os.path.exists(p)) / 2**30, 2)


def main(csv_path, fastq_dir):
    S = pd.read_csv(csv_path).dropna(how='all')
    S.columns = [c.strip().lower() for c in S.columns]
    S = S.rename(columns={'id': 'sample'})
    for c in ['sample', 'tissue', 'genotype', 'sex']:
        S[c] = S[c].astype(str).str.strip()
    S['r1'] = S['sample'].map(lambda s: parts_for(fastq_dir, s, 1))
    S['r2'] = S['sample'].map(lambda s: parts_for(fastq_dir, s, 2))
    S['n_parts'] = S.r1.map(lambda p: len([x for x in p.split(';') if x]))
    S['n_parts_r2'] = S.r2.map(lambda p: len([x for x in p.split(';') if x]))
    for side in ['r1', 'r2']:
        S[side + '_gb'] = S[side].map(gb)
    missing = S[(S.n_parts == 0) | (S.n_parts != S.n_parts_r2)]
    S['tissue_rank'] = S.tissue.map({t: i for i, t in enumerate(TISSUE_ORDER)}).fillna(99)
    S = S.sort_values(['tissue_rank', 'genotype', 'sample']).drop(columns='tissue_rank')
    os.makedirs('tables', exist_ok=True)
    S.to_csv('tables/v404m_samples.tsv', sep='\t', index=False)

    print(S.groupby(['tissue', 'genotype', 'sex']).size().unstack(fill_value=0).to_string(), file=sys.stderr)
    print(f'\n{len(S)} samples, {S.r1_gb.sum() + S.r2_gb.sum():.0f} GB of R1+R2 -> tables/v404m_samples.tsv', file=sys.stderr)
    if len(missing):
        print(f'\nWARNING: {len(missing)} samples with no FASTQ, or an R1/R2 part-count mismatch:', file=sys.stderr)
        print(missing[['sample', 'tissue', 'genotype', 'n_parts', 'n_parts_r2', 'r1_gb', 'r2_gb']].to_string(index=False), file=sys.stderr)

    multi = S[S.n_parts > 1]
    if len(multi):
        print(f'\n{len(multi)} sample(s) made of more than one lane/run — the parts will be concatenated '
              f'(check these are top-ups of the same library, not different samples):', file=sys.stderr)
        for _, r in multi.iterrows():
            print(f'  {r["sample"]}  ({r.tissue}, {r.genotype})  total R1 {r.r1_gb} GB', file=sys.stderr)
            for p in r.r1.split(';'):
                print(f'      {os.path.basename(p)}  {os.path.getsize(p) / 2**30:.2f} GB', file=sys.stderr)

    med = S.r1_gb.median()
    low = S[S.r1_gb < 0.6 * med]
    if len(low):
        print(f'\n{len(low)} sample(s) well below the median depth ({med:.1f} GB of R1) — they will still be used, '
              f'but watch them in the QC/PCA and consider dropping any that sit apart:', file=sys.stderr)
        print(low[['sample', 'tissue', 'genotype', 'sex', 'n_parts', 'r1_gb']].to_string(index=False), file=sys.stderr)
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
