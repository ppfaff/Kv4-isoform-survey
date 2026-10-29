#!/usr/bin/env python3
"""Build a stratified VISp runlist from the NeMO directory listing + Allen metadata.csv, in the same
format the hippocampal pipeline consumes (url, size, exp_component_name, sample_name, subclass_label,
cluster_label, region_label). Exact file sizes come from one HEAD request per selected file.

Usage (Lambda, in ~/allen_meta, with metadata.csv present):
    curl -s https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/VISp/ > visp_listing.html
    python3 build_visp_runlist.py visp_listing.html metadata.csv visp_runlist.tsv [cells_per_subclass=100]
"""
import sys, re, csv, random, time, urllib.request
from concurrent.futures import ThreadPoolExecutor

BASE = 'https://data.nemoarchive.org/biccn/grant/u19_zeng/zeng/transcriptome/scell/SSv4/mouse/raw/VISp/'
DRY = '--dry' in sys.argv; args = [a for a in sys.argv[1:] if a != '--dry']
listing, meta_csv, out_tsv = args[0], args[1], args[2]
PER_SUBCLASS = int(args[3]) if len(args) > 3 else 100
random.seed(1)

names = re.findall(r'href="([^"]+\.fastq\.tar)"', open(listing).read())
stems = {n[:-len('.fastq.tar')]: n for n in names}
print('files in VISp listing:', len(stems))

meta = {}
with open(meta_csv, newline='') as f:
    for row in csv.DictReader(f):
        meta[row['exp_component_name']] = row
matched = [s for s in stems if s in meta]
print('matched to metadata:', len(matched), '  unmatched:', len(stems) - len(matched))
visp = [s for s in matched if meta[s]['region_label'] == 'VISp']
print('of which region_label == VISp:', len(visp))

# stratify by subclass
by_sub = {}
for s in visp:
    by_sub.setdefault(meta[s]['subclass_label'], []).append(s)
chosen = []
print('\nsubclass           available  chosen')
for sub in sorted(by_sub, key=lambda k: -len(by_sub[k])):
    pool = by_sub[sub]
    pick = pool if len(pool) <= PER_SUBCLASS else random.sample(pool, PER_SUBCLASS)
    chosen += pick
    print(f'{sub:18s} {len(pool):9d} {len(pick):7d}')
print('total chosen:', len(chosen))

def head_size(stem, tries=4):
    url = BASE + stems[stem]
    for t in range(tries):
        try:
            req = urllib.request.Request(url, method='HEAD')
            with urllib.request.urlopen(req, timeout=60) as r:
                return stem, int(r.headers['Content-Length'])
        except Exception as e:
            time.sleep(5 * (t + 1))
    return stem, ''

print('\nfetching exact sizes (HEAD, 2 at a time)...' if not DRY else '\n--dry: skipping HEAD requests, sizes left blank')
sizes = {s: '' for s in chosen}
if not DRY:
  with ThreadPoolExecutor(max_workers=2) as ex:
    for i, (stem, sz) in enumerate(ex.map(head_size, chosen), 1):
        sizes[stem] = sz
        if i % 100 == 0: print(f'  {i}/{len(chosen)}')
missing = [s for s in chosen if sizes[s] == '']
if missing and not DRY: print('WARNING: no size for', len(missing), 'files (left blank; size check skipped for those)')

with open(out_tsv, 'w') as out:
    out.write('url\tsize\texp_component_name\tsample_name\tsubclass_label\tcluster_label\tregion_label\n')
    for s in chosen:
        m = meta[s]
        out.write('\t'.join([BASE + stems[s], str(sizes[s]), s, m['sample_name'], m['subclass_label'], m['cluster_label'], m['region_label']]) + '\n')
tot = sum(v for v in sizes.values() if v != '')
print(f'wrote {out_tsv}: {len(chosen)} cells, {tot/1e9:.0f} GB')
