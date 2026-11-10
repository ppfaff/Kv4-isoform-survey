#!/usr/bin/env python3
"""Build the ThalamoSeq (Phillips et al. 2019, SRP213359 / GSE133911) runlist and sample sheet from
tables/thalamoseq_runs.tsv (pysradb metadata SRP213359 --detailed).

Writes  tables/thalamoseq_samples.tsv   run_accession, gsm, title, nucleus, label_method (DiI/AAV), target, rep, Mreads, readlen
        tables/thalamoseq_runlist.tsv   run_accession <TAB> fastq_ftp (ENA fallback paths; blank if unknown)
Optional: --ena  queries ENA for the FASTQ FTP paths (needs network to ebi.ac.uk; run on Lambda).
Usage: python3 build_thalamoseq_runlist.py [--ena]
"""
import sys, io, os
import pandas as pd
ENA = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=run_accession,fastq_ftp&format=tsv'

t = pd.read_csv('tables/thalamoseq_runs.tsv', sep='\t', dtype=str)
t['title'] = t.experiment_title.str.extract(r'GSM\d+: (.+?); Mus')[0]
t['gsm'] = t.experiment_title.str.extract(r'(GSM\d+)')[0]
p = t.title.str.extract(r'^(?P<nucleus>[A-Za-z0-9]+)_(?P<label_method>DiI|AAV)\.(?P<target>[a-z0-9]+)_(?P<rep>\d+)$')
S = pd.concat([t[['run_accession', 'gsm', 'title']], p], axis=1)
S['Mreads'] = (t.run_total_spots.astype(float) / 1e6).round(1)
S['readlen'] = (t.run_total_bases.astype(float) / t.run_total_spots.astype(float)).round(0).astype(int)
S['fastq_ftp'] = ''
if '--ena' in sys.argv:
    import urllib.request
    rows = []
    for acc in S.run_accession:
        try:
            with urllib.request.urlopen(ENA.format(acc=acc), timeout=60) as r:
                d = pd.read_csv(io.StringIO(r.read().decode()), sep='\t', dtype=str).fillna('')
                rows.append((acc, d.fastq_ftp.iloc[0] if len(d) else ''))
        except Exception as e:
            rows.append((acc, '')); print(acc, 'ENA lookup failed:', e, file=sys.stderr)
    S = S.drop(columns='fastq_ftp').merge(pd.DataFrame(rows, columns=['run_accession', 'fastq_ftp']), on='run_accession')
os.makedirs('tables', exist_ok=True)
S.to_csv('tables/thalamoseq_samples.tsv', sep='\t', index=False)
S[['run_accession', 'fastq_ftp']].to_csv('tables/thalamoseq_runlist.tsv', sep='\t', index=False, header=False)
print(f'{len(S)} runs, {S.nucleus.nunique()} nuclei -> tables/thalamoseq_samples.tsv, tables/thalamoseq_runlist.tsv', file=sys.stderr)
print(S.groupby(['nucleus', 'label_method', 'target']).size().to_string(), file=sys.stderr)
