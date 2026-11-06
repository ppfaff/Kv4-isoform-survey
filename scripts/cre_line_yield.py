#!/usr/bin/env python3
"""Which Cre lines actually yield a given transcriptomic type, and at what purity?

There is no Tacr3-Cre, and cluster-level driver lines generally do not exist below the subclass, so the
practical question for a recording experiment is which EXISTING line gives the best yield of the target
type. Tasic 2016 answers it empirically: every cell in that dataset was FACS-sorted from a named Cre line
crossed to a tdTomato reporter, and the sample titles carry the line (e.g. Calb2_tdTpositive_cell_14).

Reports both directions, which answer different questions:
  COMPOSITION  of each type: which lines contributed its cells. Answers "can line X reach this type at all".
  YIELD        within each line: what fraction of that line's tdT+ cells landed in the target type. This is
               the number that decides how many cells must be patched to collect n of the target, and it is
               the one that matters at the rig.

IMPORTANT CAVEAT ON YIELD. Tasic 2016 sampled deliberately to cover transcriptomic diversity, with
stratification by Cre line and by dissected layer, so the cell counts are not a natural frequency. Within a
single line the cluster composition approximates what that line labels, but it is distorted by their layer
stratification, and a line sampled only from one layer will look more specific than it is. Treat yield as a
rank ordering between lines, not as an expected hit rate. tdT-negative cells are excluded.

Usage:
    CTX_DATA=/path/kv4_isoform_per_cell_annotated.csv python3 cre_line_yield.py [out_dir]
    TARGETS='Pvalb Tacr3,Sst Chodl,L6a Sla' python3 cre_line_yield.py
    ENA=/path/filereport.tsv python3 cre_line_yield.py       # use a saved report instead of fetching
"""
import os, re, sys, subprocess
import pandas as pd

OUT      = sys.argv[1] if len(sys.argv) > 1 else '.'
PROJECT  = os.environ.get('PROJECT', 'PRJNA291624')
TARGETS  = [t.strip() for t in os.environ.get(
    'TARGETS', 'Pvalb Tacr3,Sst Chodl,L6a Sla,L5a Batf3,L2/3 Ptgs2,Ndnf Car4').split(',')]
MIN_LINE = int(os.environ.get('MIN_LINE', 10))      # minimum cells from a line to report it


def get_ena():
    path = os.environ.get('ENA')
    if path and os.path.exists(path):
        return pd.read_csv(path, sep='\t')
    url = (f'https://www.ebi.ac.uk/ena/portal/api/filereport?accession={PROJECT}'
           '&result=read_run&fields=run_accession,sample_title,sample_alias&format=tsv')
    r = subprocess.run(['curl', '-fsSL', url], capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout.strip():
        sys.exit(f'ENA fetch failed; save the report and pass ENA=<path>\n{r.stderr[:400]}')
    from io import StringIO
    return pd.read_csv(StringIO(r.stdout), sep='\t')


def parse_line(title):
    """'Calb2_tdTpositive_cell_14' -> ('Calb2', True). Returns (None, None) if it does not parse."""
    m = re.match(r'^(.*?)_tdT(positive|negative)_cell', str(title), re.I)
    if not m:
        return None, None
    return m.group(1), m.group(2).lower() == 'positive'


def main():
    ctx = os.environ.get('CTX_DATA')
    if not ctx or not os.path.exists(ctx):
        sys.exit('set CTX_DATA to the Tasic 2016 annotated per-cell CSV')
    d = pd.read_csv(ctx, usecols=lambda c: c in (
        'sample_name', 'primary_type_label', 'cluster_label', 'subclass', 'cluster',
        'dissection', 'layer', 'broad_type'))
    tcol = next((c for c in ('primary_type_label', 'cluster_label', 'cluster', 'subclass')
                 if c in d.columns), None)
    if tcol is None:
        sys.exit(f'no type column found; columns were {list(d.columns)}')
    lcol = next((c for c in ('dissection', 'layer') if c in d.columns), None)
    bcol = 'broad_type' if 'broad_type' in d.columns else None
    keep = ['sample_name', tcol] + [c for c in (lcol, bcol) if c]
    cells = d[keep].drop_duplicates().rename(columns={tcol: 'cell_type'})
    if lcol:
        cells = cells.rename(columns={lcol: 'layer'})
    if bcol:
        cells = cells.rename(columns={bcol: 'broad_type'})

    E = get_ena()
    E['cre_line'], E['tdT_pos'] = zip(*E['sample_title'].map(parse_line))
    E = E[E.tdT_pos == True].copy()
    # the per-cell table may key on the sample title or the run accession
    for key in ('sample_title', 'run_accession', 'sample_alias'):
        m = cells.merge(E[[key, 'cre_line']], left_on='sample_name', right_on=key, how='inner')
        if len(m) > 0.5 * len(cells):
            print(f'joined on {key}: {len(m):,} of {len(cells):,} cells matched to a Cre line')
            break
    else:
        sys.exit(f'could not join; sample_name looks like "{cells.sample_name.iloc[0]}" '
                 f'and ENA titles look like "{E.sample_title.iloc[0]}"')

    os.makedirs(f'{OUT}/tables', exist_ok=True)
    m.to_csv(f'{OUT}/tables/cre_line_per_cell.csv', index=False)

    print('\n' + '=' * 96)
    print('COMPOSITION: which Cre lines contributed the cells of each target type')
    print('=' * 96)
    for t in TARGETS:
        sub = m[m.cell_type == t]
        if sub.empty:
            print(f'\n{t}: not present under that name'); continue
        vc = sub.cre_line.value_counts()
        print(f'\n{t}  (n = {len(sub)})')
        for line, n in vc.items():
            print(f'    {line:22} {n:4d} cells  ({n/len(sub):5.1%} of this type)')

    print('\n' + '=' * 96)
    print('YIELD: of each line\'s tdT+ cells, what fraction landed in the target type')
    print('  rank ordering only - Tasic 2016 stratified its sampling, so these are not hit rates')
    print('=' * 96)
    tot = m.cre_line.value_counts()
    rows = []
    for t in TARGETS:
        sub = m[m.cell_type == t]
        if sub.empty:
            continue
        for line, n in sub.cre_line.value_counts().items():
            if tot[line] >= MIN_LINE:
                rows.append(dict(cell_type=t, cre_line=line, n_in_type=n,
                                 n_line_total=int(tot[line]), yield_pct=n / tot[line] * 100))
    Y = pd.DataFrame(rows).sort_values(['cell_type', 'yield_pct'], ascending=[True, False])
    Y.to_csv(f'{OUT}/tables/cre_line_yield.csv', index=False)
    for t, sub in Y.groupby('cell_type'):
        print(f'\n{t}')
        for _, r in sub.head(6).iterrows():
            print(f'    {r.cre_line:22} {r.yield_pct:5.1f}%  ({int(r.n_in_type)} of '
                  f'{int(r.n_line_total)} cells from this line)')

    # ---- line x layer: does restricting the dissection improve yield?
    if 'layer' in m.columns:
        print('\n' + '=' * 96)
        print('LINE x LAYER: does restricting the dissected layer improve yield?')
        print('  A pan-GABAergic line (Gad2, Slc32a1/Viaat) is only worth using if some layer')
        print('  concentrates the target enough to beat the subclass-specific line.')
        print('=' * 96)
        for t in TARGETS:
            sub = m[m.cell_type == t]
            if sub.empty:
                continue
            print(f'\n{t}')
            g = (m.groupby(['cre_line', 'layer'])
                   .agg(n_total=('cell_type', 'size'),
                        n_target=('cell_type', lambda x, tt=t: int((x == tt).sum())))
                   .reset_index())
            g = g[(g.n_total >= MIN_LINE) & (g.n_target > 0)]
            g['yield_pct'] = g.n_target / g.n_total * 100
            for _, r in g.sort_values('yield_pct', ascending=False).head(8).iterrows():
                print(f'    {r.cre_line:20} {str(r.layer)[:18]:18} {r.yield_pct:5.1f}%  '
                      f'({int(r.n_target)} of {int(r.n_total)})')
        # what a pan-GABAergic line would give, by layer
        if 'broad_type' in m.columns:
            gaba = m[m.broad_type.astype(str).str.contains('GABA', case=False, na=False)]
            if not gaba.empty:
                print('\n' + '-' * 96)
                print('If a pan-GABAergic driver (Viaat/Slc32a1, Gad2) is used, what each layer yields:')
                print('-' * 96)
                for t in [x for x in TARGETS if x in set(gaba.cell_type)]:
                    sub = (gaba.assign(hit=gaba.cell_type == t)
                               .groupby('layer').agg(n=('hit', 'size'), k=('hit', 'sum')))
                    sub = sub[sub.n >= MIN_LINE]
                    if sub.empty:
                        continue
                    best = (sub.k / sub.n * 100).sort_values(ascending=False)
                    line = '  '.join(f'{str(i)[:14]} {v:.0f}% ({int(sub.k[i])}/{int(sub.n[i])})'
                                     for i, v in best.head(4).items())
                    print(f'  {t:16} {line}')

    print('\n' + '=' * 96)
    print('All lines in the dataset, by total cells (for reference)')
    print('=' * 96)
    print(tot[tot >= MIN_LINE].to_string())
    print(f'\nwrote {OUT}/tables/cre_line_{{per_cell,yield}}.csv')


if __name__ == '__main__':
    main()
