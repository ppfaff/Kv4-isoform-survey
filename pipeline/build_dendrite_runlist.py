#!/usr/bin/env python3
"""Build the run list and sample sheet for the dendrite-vs-soma datasets.

Datasets
  farris    Farris et al. 2019 (GEO GSE116342 / BioProject PRJNA478217): adult mouse hippocampus,
            laser-microdissected cell-body layer (CB) vs dendritic layer (DE) for CA1, CA2, CA3, DG;
            3 animals -> 24 stranded paired-end libraries.
  middleton Middleton et al. 2019 (GEO GSE115480 / BioProject PRJNA475141 / SRP150011): 16 cultured
            hippocampal neurons (E18 + 15 DIV), each micro-dissected into a soma and a dendrite sample
            -> 32 stranded paired-end libraries.

Metadata comes from two places, because ENA's sample_title is sometimes blank for GEO submissions:
  * ENA filereport (run -> experiment, biosample, layout, read counts, fastq paths)
  * GEO series matrix (GSM -> title, characteristics, and the SRX it was submitted under)
They are joined on the SRX experiment accession. The sample sheet carries the parsed design so
the downstream scripts never have to interpret titles themselves.

Outputs (in ./tables):
  dendrite_runlist.tsv   run_accession <TAB> fastq_ftp      (what run_dendrite_isoforms.sh consumes)
  dendrite_samples.tsv   one row per run with: dataset, gsm, srx, biosample, title, region, compartment,
                         unit (animal id or neuron id), layout, read_count, plus the raw characteristics
Rows whose title could not be parsed are written with region/compartment = '?' and listed on stderr —
fix those by hand in dendrite_samples.tsv before running the analysis (the collapse step refuses '?').

Usage (on Lambda, in the work dir, with network):
    python3 build_dendrite_runlist.py            # both datasets
    python3 build_dendrite_runlist.py farris     # one dataset
"""
import sys, re, io, gzip, os, urllib.request
import pandas as pd

DATASETS = {
    'farris':    {'project': 'PRJNA478217', 'gse': 'GSE116342'},
    'middleton': {'project': 'PRJNA475141', 'gse': 'GSE115480'},
}
ENA = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession={p}&result=read_run'
       '&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,'
       'sample_alias,sample_title,library_layout,library_strategy,instrument_model,read_count,base_count,fastq_ftp'
       '&format=tsv&download=true')
GEO = 'https://ftp.ncbi.nlm.nih.gov/geo/series/{stem}nnn/{gse}/matrix/{gse}_series_matrix.txt.gz'


def fetch(url, binary=False):
    req = urllib.request.Request(url, headers={'User-Agent': 'kv4-isoform-survey/1.0'})
    with urllib.request.urlopen(req, timeout=120) as r:
        data = r.read()
    return data if binary else data.decode()


def ena_runs(project):
    t = pd.read_csv(io.StringIO(fetch(ENA.format(p=project))), sep='\t', dtype=str).fillna('')
    t = t[t.library_strategy.str.upper().isin(['RNA-SEQ', 'OTHER', ''])]
    return t


def geo_samples(gse):
    """Parse a GEO series matrix into one row per GSM: title, characteristics (joined), SRX."""
    stem = gse[:-3]
    raw = gzip.decompress(fetch(GEO.format(stem=stem, gse=gse), binary=True)).decode(errors='replace')
    rows = {}
    for line in raw.splitlines():
        if not line.startswith('!Sample_'):
            continue
        key, _, rest = line.partition('\t')
        vals = [v.strip().strip('"') for v in rest.split('\t')]
        rows.setdefault(key, []).append(vals)
    gsm = rows['!Sample_geo_accession'][0]
    title = rows['!Sample_title'][0]
    chars = ['; '.join(filter(None, [c[i] for c in rows.get('!Sample_characteristics_ch1', [])]))
             for i in range(len(gsm))]
    srx = []
    for i in range(len(gsm)):
        s = ''
        for rel in rows.get('!Sample_relation', []):
            m = re.search(r'SRA:\s*\S*?(SRX\d+)', rel[i])
            if m:
                s = m.group(1)
        srx.append(s)
    source = [c[0] for c in rows.get('!Sample_source_name_ch1', [[''] * len(gsm)])] if rows.get('!Sample_source_name_ch1') else [''] * len(gsm)
    return pd.DataFrame({'gsm': gsm, 'title': title, 'characteristics': chars, 'source_name': source, 'srx': srx})


# ---------------------------------------------------------------- design parsers
REGION_RE = re.compile(r'\b(CA1|CA2|CA3|DG)\b', re.I)


def parse_farris(text):
    """Return (region, compartment, animal) from a Farris title/characteristics string.
    Compartment: cell-body layer -> 'soma'; dendritic layer (s. radiatum / molecular layer / neuropil) -> 'dendrite'."""
    t = text.replace('_', ' ')
    rm = REGION_RE.search(t)
    region = rm.group(1).upper() if rm else '?'
    tl = t.lower()
    if re.search(r'\bcb\b|cell.?body|pyramidal|granule.?cell|somat', tl):
        comp = 'soma'
    elif re.search(r'\bde\b|dendrit|radiatum|molecular|neuropil|\bsr\b|\bml\b', tl):
        comp = 'dendrite'
    else:
        comp = '?'
    am = re.search(r'(?:animal|mouse|rep(?:licate)?|biological)[^0-9]*(\d+)', tl) or re.search(r'(\d{3,4})\b', t)
    animal = am.group(1) if am else '?'
    return region, comp, animal


def parse_middleton(text):
    """Return ('culture', compartment, neuron id) from a Middleton title/characteristics string."""
    t = text.replace('_', ' ')
    tl = t.lower()
    if re.search(r'dendrit|neurite|process', tl):
        comp = 'dendrite'
    elif re.search(r'soma|cell.?body', tl):
        comp = 'soma'
    else:
        comp = '?'
    nm = (re.search(r'(?:neuron|cell|pair|n)[^0-9a-z]*(\d+)', tl) or re.search(r'(\d+)', tl))
    neuron = nm.group(1) if nm else '?'
    return 'culture', comp, neuron


PARSERS = {'farris': parse_farris, 'middleton': parse_middleton}


def build(name):
    cfg = DATASETS[name]
    runs = ena_runs(cfg['project'])
    print(f'{name}: {len(runs)} ENA runs for {cfg["project"]}', file=sys.stderr)
    try:
        geo = geo_samples(cfg['gse'])
        print(f'{name}: {len(geo)} GEO samples for {cfg["gse"]} ({(geo.srx != "").sum()} with SRX links)', file=sys.stderr)
    except Exception as e:  # GEO ftp occasionally refuses; ENA titles alone may still be enough
        print(f'{name}: GEO series matrix unavailable ({e}); using ENA titles only', file=sys.stderr)
        geo = pd.DataFrame(columns=['gsm', 'title', 'characteristics', 'source_name', 'srx'])
    m = runs.merge(geo, left_on='experiment_accession', right_on='srx', how='left').fillna('')
    # text to parse: GEO title + characteristics + source, falling back to ENA title/alias
    m['text'] = (m.title + ' | ' + m.characteristics + ' | ' + m.source_name + ' | ' + m.sample_title + ' | ' + m.sample_alias)
    parsed = m.text.map(PARSERS[name])
    m['dataset'] = name
    m['region'] = [p[0] for p in parsed]
    m['compartment'] = [p[1] for p in parsed]
    m['unit'] = [p[2] for p in parsed]
    return m


def main(names):
    os.makedirs('tables', exist_ok=True)
    parts = [build(n) for n in names]
    S = pd.concat(parts, ignore_index=True)
    cols = ['dataset', 'run_accession', 'experiment_accession', 'sample_accession', 'gsm', 'title', 'region', 'compartment', 'unit',
            'library_layout', 'instrument_model', 'read_count', 'base_count', 'characteristics', 'sample_title', 'sample_alias', 'fastq_ftp']
    S = S[[c for c in cols if c in S.columns]]
    S.to_csv('tables/dendrite_samples.tsv', sep='\t', index=False)
    S[['run_accession', 'fastq_ftp']].to_csv('tables/dendrite_runlist.tsv', sep='\t', index=False, header=False)
    bad = S[(S.region == '?') | (S.compartment == '?') | (S.unit == '?')]
    print(f'\n{len(S)} runs -> tables/dendrite_runlist.tsv; sample sheet -> tables/dendrite_samples.tsv', file=sys.stderr)
    print(S.groupby(['dataset', 'region', 'compartment']).size().to_string(), file=sys.stderr)
    if len(bad):
        print(f'\n{len(bad)} runs could not be fully parsed — edit region/compartment/unit in tables/dendrite_samples.tsv:', file=sys.stderr)
        print(bad[['dataset', 'run_accession', 'gsm', 'title', 'characteristics', 'sample_title']].to_string(index=False), file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1:] or list(DATASETS))
