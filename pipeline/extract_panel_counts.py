#!/usr/bin/env python3
"""Extract per-cell RAW READ COUNTS for the Kv4-complex panel from salmon quant.sf files.

Run this where the quant directory lives - for the Allen/Tasic cortex set that is inside WSL at
~/allen_meta - then copy the single small CSV it writes to a Windows folder Claude can reach.

    cd ~/allen_meta
    curl -sSLO https://raw.githubusercontent.com/ppfaff/Kv4-isoform-survey/main/pipeline/extract_panel_counts.py
    python3 extract_panel_counts.py quant panel_counts_per_cell.csv
    cp panel_counts_per_cell.csv /mnt/c/Users/paulp/

Pure standard library - no pandas, no numpy, nothing to install.

WHY RAW COUNTS. Every per-cell table currently in the repo carries CPM or TPM. A ratio of two
normalised values cannot be given a sensible confidence interval, which is why the Dirichlet guard on
the KChIP:Kv4 ratio is inert: CPM was used as pseudo-counts and overstates the evidence in shallow
libraries. salmon's NumReads column is the quantity that supports a posterior, and quant.sf has it.

It also writes num_processed and percent_mapped per run from aux_info/meta_info.json, so a low ratio
can be judged against how deeply that cell was sequenced.

IF A TRANSCRIPT IS MISSING. The panel map embedded below was built from transcripts detected in the
dendrite dataset against GENCODE vM36. If the reference FASTA is present (ref/gencode.vM36.transcripts.fa.gz
or passed as a third argument) the script rebuilds the map from it instead, which is authoritative and
catches any panel transcript the embedded list lacks. It reports which source it used.
"""
import csv, glob, gzip, json, os, sys

PANEL_GENES = {'Kcnd1', 'Kcnd2', 'Kcnd3', 'Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4', 'Dpp6', 'Dpp10'}

# ENSMUST (version-stripped) -> (gene, transcript, vclass, panel)
EMBEDDED = {
    'ENSMUST00000112603': ('Dpp10', 'Dpp10-201', 'DPP10d', 'DPP'),
    'ENSMUST00000112606': ('Dpp10', 'Dpp10-202', 'DPP10c', 'DPP'),
    'ENSMUST00000140361': ('Dpp10', 'Dpp10-203', 'DPP10d', 'DPP'),
    'ENSMUST00000187202': ('Dpp10', 'Dpp10-204', 'Dpp10 other/fragment', 'DPP'),
    'ENSMUST00000187286': ('Dpp10', 'Dpp10-205', 'DPP10a', 'DPP'),
    'ENSMUST00000239072': ('Dpp10', 'Dpp10-206', 'DPP10a', 'DPP'),
    'ENSMUST00000071500': ('Dpp6', 'Dpp6-201', 'DPP6a', 'DPP'),
    'ENSMUST00000101471': ('Dpp6', 'Dpp6-202', 'DPP6S', 'DPP'),
    'ENSMUST00000120555': ('Dpp6', 'Dpp6-203', 'DPP6K', 'DPP'),
    'ENSMUST00000122171': ('Dpp6', 'Dpp6-204', 'DPP6L', 'DPP'),
    'ENSMUST00000134175': ('Dpp6', 'Dpp6-205', 'Dpp6 other/fragment', 'DPP'),
    'ENSMUST00000136574': ('Dpp6', 'Dpp6-206', 'DPP6a', 'DPP'),
    'ENSMUST00000148039': ('Dpp6', 'Dpp6-207', 'DPP6S', 'DPP'),
    'ENSMUST00000149678': ('Dpp6', 'Dpp6-208', 'DPP6a', 'DPP'),
    'ENSMUST00000200274': ('Dpp6', 'Dpp6-209', 'DPP6L', 'DPP'),
    'ENSMUST00000009875': ('Kcnd1', 'Kcnd1-201', 'Kv4.1', 'Kv4'),
    'ENSMUST00000125440': ('Kcnd1', 'Kcnd1-202', '', 'Kv4'),
    'ENSMUST00000081542': ('Kcnd2', 'Kcnd2-201', 'Kv4.2', 'Kv4'),
    'ENSMUST00000079169': ('Kcnd3', 'Kcnd3-201', 'Kv4.3L', 'Kv4'),
    'ENSMUST00000098761': ('Kcnd3', 'Kcnd3-202', 'Kv4.3S', 'Kv4'),
    'ENSMUST00000118360': ('Kcnd3', 'Kcnd3-203', 'Kv4.3L', 'Kv4'),
    'ENSMUST00000141694': ('Kcnd3', 'Kcnd3-204', '', 'Kv4'),
    'ENSMUST00000147209': ('Kcnd3', 'Kcnd3-205', '', 'Kv4'),
    'ENSMUST00000065970': ('Kcnip1', 'Kcnip1-201', 'KChIP1b', 'KChIP'),
    'ENSMUST00000101368': ('Kcnip1', 'Kcnip1-202', 'Kcnip1-new_short', 'KChIP'),
    'ENSMUST00000109340': ('Kcnip1', 'Kcnip1-203', 'KChIP1a', 'KChIP'),
    'ENSMUST00000135034': ('Kcnip1', 'Kcnip1-204', 'Kcnip1 other/fragment', 'KChIP'),
    'ENSMUST00000154760': ('Kcnip1', 'Kcnip1-205', 'Kcnip1 other/fragment', 'KChIP'),
    'ENSMUST00000238436': ('Kcnip1', 'Kcnip1-206', 'KChIP1e', 'KChIP'),
    'ENSMUST00000238903': ('Kcnip1', 'Kcnip1-207', 'KChIP1c', 'KChIP'),
    'ENSMUST00000026247': ('Kcnip2', 'Kcnip2-201', 'KChIP2b', 'KChIP'),
    'ENSMUST00000079431': ('Kcnip2', 'Kcnip2-202', 'KChIP2b', 'KChIP'),
    'ENSMUST00000086993': ('Kcnip2', 'Kcnip2-203', 'Kcnip2 other/fragment', 'KChIP'),
    'ENSMUST00000111906': ('Kcnip2', 'Kcnip2-204', 'KChIP2b', 'KChIP'),
    'ENSMUST00000159210': ('Kcnip2', 'Kcnip2-205', 'Kcnip2 other/fragment', 'KChIP'),
    'ENSMUST00000159245': ('Kcnip2', 'Kcnip2-206', 'Kcnip2 other/fragment', 'KChIP'),
    'ENSMUST00000159446': ('Kcnip2', 'Kcnip2-207', 'Kcnip2 other/fragment', 'KChIP'),
    'ENSMUST00000161886': ('Kcnip2', 'Kcnip2-208', 'KChIP2c', 'KChIP'),
    'ENSMUST00000162528': ('Kcnip2', 'Kcnip2-209', 'KChIP2a', 'KChIP'),
    'ENSMUST00000162661': ('Kcnip2', 'Kcnip2-210', 'KChIP2c', 'KChIP'),
    'ENSMUST00000239049': ('Kcnip2', 'Kcnip2-211', 'KChIP2x', 'KChIP'),
    'ENSMUST00000028850': ('Kcnip3', 'Kcnip3-201', 'Kcnip3-new-cys_rich', 'KChIP'),
    'ENSMUST00000088538': ('Kcnip3', 'Kcnip3-202', 'KChIP3x', 'KChIP'),
    'ENSMUST00000103215': ('Kcnip3', 'Kcnip3-203', 'KChIP3a', 'KChIP'),
    'ENSMUST00000137625': ('Kcnip3', 'Kcnip3-204', 'Kcnip3 other/fragment', 'KChIP'),
    'ENSMUST00000249876': ('Kcnip3', 'Kcnip3-205', 'Kcnip3 other/fragment', 'KChIP'),
    'ENSMUST00000087395': ('Kcnip4', 'Kcnip4-201', 'KChIP4bL', 'KChIP'),
    'ENSMUST00000101214': ('Kcnip4', 'Kcnip4-202', 'Kcnip4 other/noncoding', 'KChIP'),
    'ENSMUST00000101215': ('Kcnip4', 'Kcnip4-203', 'Kcnip4 other/noncoding', 'KChIP'),
    'ENSMUST00000166924': ('Kcnip4', 'Kcnip4-204', 'KChIP4e', 'KChIP'),
    'ENSMUST00000172363': ('Kcnip4', 'Kcnip4-205', 'Kcnip4 other/noncoding', 'KChIP'),
    'ENSMUST00000175660': ('Kcnip4', 'Kcnip4-206', 'KChIP4d', 'KChIP'),
    'ENSMUST00000176191': ('Kcnip4', 'Kcnip4-207', 'KChIP4b', 'KChIP'),
    'ENSMUST00000176978': ('Kcnip4', 'Kcnip4-208', 'KChIP4a', 'KChIP'),
}


def map_from_fasta(fa):
    """Authoritative map straight from the GENCODE transcript FASTA headers."""
    out = {}
    op = gzip.open if fa.endswith('.gz') else open
    with op(fa, 'rt') as f:
        for line in f:
            if line[:1] != '>':
                continue
            p = line[1:].strip().split('|')
            if len(p) > 5 and p[5] in PANEL_GENES:
                out[p[0].split('.')[0]] = (p[5], p[4] if len(p) > 4 else '', '', '')
    return out


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else 'quant'
    out_csv = sys.argv[2] if len(sys.argv) > 2 else 'panel_counts_per_cell.csv'
    fa = sys.argv[3] if len(sys.argv) > 3 else 'ref/gencode.vM36.transcripts.fa.gz'

    tx, src = dict(EMBEDDED), 'embedded map (%d transcripts)' % len(EMBEDDED)
    if os.path.exists(fa):
        full = map_from_fasta(fa)
        if full:
            for k, v in full.items():
                if k not in tx:                       # keep embedded vclass labels where we have them
                    tx[k] = v
            src = 'reference FASTA %s (%d panel transcripts; %d added beyond the embedded map)' % (
                fa, len(full), len(set(full) - set(EMBEDDED)))
    print('transcript map: ' + src)

    quants = sorted(glob.glob(os.path.join(root, '*', 'quant.sf')))
    if not quants:
        quants = sorted(glob.glob(os.path.join(root, '**', 'quant.sf'), recursive=True))
    if not quants:
        sys.exit('no quant.sf found under %r - pass the quant directory as the first argument' % root)
    print('found %d quant.sf files under %s' % (len(quants), root))

    n_rows = 0
    missing_meta = 0
    with open(out_csv, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['run', 'ensmust', 'gene', 'transcript', 'vclass', 'panel',
                    'NumReads', 'TPM', 'EffectiveLength', 'num_processed', 'percent_mapped'])
        for i, q in enumerate(quants, 1):
            run = os.path.basename(os.path.dirname(q))
            npz = pm = ''
            mj = os.path.join(os.path.dirname(q), 'aux_info', 'meta_info.json')
            try:
                m = json.load(open(mj))
                npz, pm = m.get('num_processed', ''), m.get('percent_mapped', '')
            except Exception:
                missing_meta += 1
            with open(q) as f:
                r = csv.DictReader(f, delimiter='\t')
                for row in r:
                    key = row['Name'].split('|')[0].split('.')[0]
                    hit = tx.get(key)
                    if hit is None:
                        continue
                    g, t, vc, pan = hit
                    w.writerow([run, row['Name'], g, t, vc, pan,
                                row['NumReads'], row['TPM'], row['EffectiveLength'], npz, pm])
                    n_rows += 1
            if i % 200 == 0:
                print('  %d/%d runs...' % (i, len(quants)), flush=True)

    print('wrote %s: %d rows from %d runs' % (out_csv, n_rows, len(quants)))
    if missing_meta:
        print('  note: %d runs had no aux_info/meta_info.json, so depth is blank for those' % missing_meta)
    exp = len(quants) * len(tx)
    if n_rows < 0.5 * exp:
        print('  WARNING: %d rows is well under the %d expected if every panel transcript appeared in'
              ' every run - check that the quant used the same reference.' % (n_rows, exp))


if __name__ == '__main__':
    main()
