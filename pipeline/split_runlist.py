import csv, glob, sys

RUNLIST = 'hip_filtered_runlist.tsv'
# core-weighted shares: Alienware 8, PFAFF-LAMBDA 32, PFAFF10 8
WEIGHTS = {'alienware': 8, 'lambda': 32, 'pfaff10': 8}
MACHINES = ['alienware', 'lambda', 'pfaff10']  # alienware's slice is taken first (it's furthest along)

# cells already fully quantified here (Alienware) -> keep them credited to alienware,
# no need to ever touch them on another machine
done = {d.split('/')[1] for d in glob.glob('quant/*/quant.sf')}

with open(RUNLIST, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    header = next(r)
    rows = list(r)

sample_col = header.index('sample_name')
already_done_rows = [row for row in rows if row[sample_col] in done]
remaining = [row for row in rows if row[sample_col] not in done]

total_weight = sum(WEIGHTS.values())
n = len(remaining)
counts, acc = {}, 0
for i, m in enumerate(MACHINES):
    if i == len(MACHINES) - 1:
        counts[m] = n - acc
    else:
        c = round(n * WEIGHTS[m] / total_weight)
        counts[m] = c
        acc += c

out, idx = {}, 0
for m in MACHINES:
    out[m] = remaining[idx:idx + counts[m]]
    idx += counts[m]
# alienware also keeps its already-completed rows (harmless - quant.sf check skips them on restart)
out['alienware'] = already_done_rows + out['alienware']

for m, rws in out.items():
    fname = f'hip_runlist_{m}.tsv'
    with open(fname, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for row in rws:
            f.write('\t'.join(row) + '\n')
    print(f'{m}: {len(rws)} cells -> {fname}')

print(f'\ntotal rows in: {len(rows)}  (already done: {len(already_done_rows)}, remaining split: {n})')
