#!/usr/bin/env python3
"""Known-ratio calibration of salmon's variant assignment, in silico.

This is the computational version of the standardizing constructs in Jerng & Pfaffinger 2012 (PLoS ONE
7:e38205): there, tandem plasmids carrying GAPDH exons 5-8 followed by exons 1-3 of two DPLP variants gave a
guaranteed 1:1:1 molar ratio, so each primer set's amplification efficiency could be measured against a known
answer and cross-variant ratios reported as true molar ratios. salmon has the identical problem and no such
control: its per-isoform estimates come from an effective-length model, a mappability model and an EM that
apportions ambiguous reads, and nothing in this project has ever checked them against a known ratio.

Why this rather than bootstraps: bootstraps resample reads UNDER THE SAME MODEL, so they measure precision
and are structurally blind to systematic model error. The 4-fold Kv4.3L underestimate traced to 3'UTR model
asymmetry would not appear in a bootstrap at all. A known input ratio measures accuracy.

THE TEST. Simulate equal MOLAR amounts of each variant class (so reads scale with transcript length, as a
plasmid standard would), quantify against the project's own salmon index, and compare. Equal molarity in
means equal TPM out, since TPM is molarity-like -- so any departure from a flat TPM profile is bias,
measured in the same units as every result in this project.

It sweeps the 3' coverage bias, because that is the suspected mechanism. The DPP6 and DPP10 N-terminal
variants are distinguished almost entirely by 5' sequence, and the two datasets here disagree in exactly the
way that predicts: full-length SMART-seq (Tasic 2016) puts DPP10a at 62.5% of the DPP panel in excitatory
cells, matching the known biology that DPP10a is the dominant N-type variant in cortical pyramidal neurons,
while 3'-biased polyA bulk (V404M) reports 12.4% and piles the reads onto DPP10d and DPP6S. If DPP10a
recovery collapses as tau shrinks and holds at tau=inf, the mechanism is confirmed and the sweep gives the
correction factor.

Usage:
    python3 simulate_variant_calibration.py --fasta ref/gencode.vM36.transcripts.fa.gz \
        --index ref/salmon_idx_tx --map tables/transcript_vclass_map.tsv --out calib
Options:
    --mol N          molecules per variant class (default 20000)
    --readlen N      read length (default 100)
    --taus a,b,c     3'-bias decay constants in nt; 'inf' is uniform full-length coverage
                     (default: inf,2000,1000,500,250)
    --layout         se | pe | both   (default both: SMART-seq is SE 100, the V404M bulk is PE)
"""
import os, sys, gzip, argparse, subprocess, random, json
import numpy as np

PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}
VAR_PANEL = {v: p for p, vs in PANEL_V.items() for v in vs}
RC = str.maketrans('ACGTN', 'TGCAN')


def revcomp(s):
    return s.translate(RC)[::-1]


def load_map(path):
    tx2v = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[1] in VAR_PANEL:
                tx2v[p[0]] = p[1]
    if not tx2v:
        sys.exit(f'no usable entries in {path}')
    return tx2v


def load_seqs(fasta, tx2v):
    """Pull the sequences of every transcript belonging to a variant class. GENCODE headers are
    pipe-delimited, so match on the full name, any pipe field, and the version-stripped id."""
    want, seqs, name, buf = {}, {}, None, []
    op = gzip.open if fasta.endswith('.gz') else open
    with op(fasta, 'rt') as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None and name in want:
                    seqs[name] = ''.join(buf).upper()
                hdr = line[1:].strip()
                full = hdr.split()[0]
                cand = [full] + full.split('|') + [full.split('|')[0].split('.')[0]]
                v = next((tx2v[c] for c in cand if c in tx2v), None)
                name = full if v else None
                if v:
                    want[full] = v
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name is not None and name in want:
        seqs[name] = ''.join(buf).upper()
    return {t: (want[t], s) for t, s in seqs.items() if len(s) > 200}


def simulate(seqs, mol, readlen, tau, layout, rng, fraglen=250, fragsd=40, err=0.001):
    """Equal molarity per VARIANT CLASS, split evenly over that class's transcripts.

    tau = inf gives uniform (full-length) coverage. Finite tau anchors each fragment's 3' end at a distance
    from the transcript 3' end drawn from Exponential(tau), which is how a polyA-primed, partially degraded
    library concentrates coverage near the 3' end.
    """
    by_class = {}
    for t, (v, s) in seqs.items():
        by_class.setdefault(v, []).append((t, s))
    r1, r2, truth = [], [], {}
    for v, txs in by_class.items():
        per_tx = max(1, mol // len(txs))
        truth[v] = 0
        for t, s in txs:
            L = len(s)
            n = 0
            # A molecule of a long transcript yields proportionally more fragments than a short one, so
            # reads must scale with length for equal MOLARITY to be what goes in. Without this the
            # simulation puts in equal READS per class, under which a flat TPM profile would be wrong by
            # construction and the whole calibration would be measuring the bug instead of salmon.
            n_frag = max(1, int(round(per_tx * L / 1000.0)))
            for _ in range(n_frag):
                fl = max(readlen + 10, int(rng.normal(fraglen, fragsd)))
                if fl >= L:
                    start = 0; fl = L
                elif np.isinf(tau):
                    start = rng.integers(0, L - fl + 1)
                else:
                    d = rng.exponential(tau)          # distance of fragment 3' end from transcript 3' end
                    end = int(L - d)
                    start = end - fl
                    if start < 0 or end > L:
                        continue                      # fragment fell off the molecule
                frag = s[start:start + fl]
                if len(frag) < readlen:
                    continue
                a = frag[:readlen]
                b = revcomp(frag[-readlen:])
                if err > 0:
                    a = mutate(a, err, rng); b = mutate(b, err, rng)
                n += 1
                rid = f'{v}:{t}:{n}'
                r1.append((rid, a))
                if layout == 'pe':
                    r2.append((rid, b))
            truth[v] += n
    return r1, r2, truth


def mutate(seq, rate, rng):
    k = rng.binomial(len(seq), rate)
    if k == 0:
        return seq
    s = list(seq)
    for i in rng.choice(len(seq), size=k, replace=False):
        s[i] = 'ACGT'[rng.integers(0, 4)]
    return ''.join(s)


def write_fq(path, reads):
    with gzip.open(path, 'wt') as fh:
        for rid, s in reads:
            fh.write(f'@{rid}\n{s}\n+\n{"I" * len(s)}\n')


def run_salmon(index, r1p, r2p, outdir, threads):
    cmd = ['salmon', 'quant', '-i', index, '-l', 'A', '-p', str(threads), '-o', outdir, '-q']
    cmd += (['-1', r1p, '-2', r2p] if r2p else ['-r', r1p])
    p = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(os.path.join(outdir, 'quant.sf')):
        sys.stderr.write(p.stderr[-2000:] + '\n')
        return None
    return os.path.join(outdir, 'quant.sf')


def read_quant(sf, tx2v):
    tot = {}
    with open(sf) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        iN, iT = hdr.index('Name'), hdr.index('TPM')
        iR = hdr.index('NumReads')
        for line in fh:
            f = line.rstrip('\n').split('\t')
            n = f[iN]
            cand = [n] + n.split('|') + [n.split('|')[0].split('.')[0]]
            v = next((tx2v[c] for c in cand if c in tx2v), None)
            if v:
                a, b = tot.get(v, (0.0, 0.0))
                tot[v] = (a + float(f[iT]), b + float(f[iR]))
    return tot


def _emit_purity(conf, out):
    if not conf:
        return
    import pandas as pd
    C = pd.DataFrame(conf)
    C.to_csv(os.path.join(out, 'variant_cross_assignment.csv'), index=False)
    for (panel, layout, tau), sub in C.groupby(['panel', 'layout', 'tau']):
        M = sub.pivot_table(index='true_variant', columns='assigned_variant', values='pct')
        print('\n' + '=' * 88)
        print(f'{panel}  {layout}  tau={tau}: rows = TRUE origin, columns = where salmon put the reads (%)')
        print('=' * 88)
        print(M.round(1).to_string())
        off = {c: float(M[c].drop(index=c).sum()) if c in M.index else float(M[c].sum())
               for c in M.columns}
        sink = max(off, key=off.get)
        print(f'  diagonal (correctly assigned): '
              f'{", ".join(f"{i} {M.loc[i, i]:.0f}%" for i in M.index if i in M.columns)}')
        print(f'  biggest SINK: {sink} absorbs {off[sink]:.0f}% summed over other variants\' reads')
    print(f"\nwrote {out}/variant_cross_assignment.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--index', required=True)
    ap.add_argument('--map', default='tables/transcript_vclass_map.tsv')
    ap.add_argument('--out', default='calib')
    ap.add_argument('--mol', type=int, default=20000)
    ap.add_argument('--readlen', type=int, default=100)
    ap.add_argument('--taus', default='inf,2000,1000,500,250')
    ap.add_argument('--layout', default='both', choices=['se', 'pe', 'both'])
    ap.add_argument('--threads', type=int, default=6)
    ap.add_argument('--panels', default='DPP,Kv4,KChIP')
    ap.add_argument('--mode', default='both', choices=['mix', 'purity', 'both'],
                    help="mix: equimolar, is the panel profile flat. purity: one variant at "
                         "a time, where do its reads actually go (sink detection).")
    a = ap.parse_args()

    rng = np.random.default_rng(0)
    os.makedirs(a.out, exist_ok=True)
    tx2v = load_map(a.map)
    panels = a.panels.split(',')
    tx2v = {t: v for t, v in tx2v.items() if VAR_PANEL[v] in panels}
    seqs = load_seqs(a.fasta, tx2v)
    print(f'{len(seqs)} transcripts across {len(set(v for v, _ in seqs.values()))} variant classes')
    if not seqs:
        sys.exit('no sequences matched -- check --fasta and --map')

    taus = [float('inf') if t.strip() == 'inf' else float(t) for t in a.taus.split(',')]
    layouts = ['se', 'pe'] if a.layout == 'both' else [a.layout]
    rows = []

    # ---------------------------------------------------------------- purity: where do PURE reads go?
    # Tests directly whether one variant acts as the default bin for reads that are identifiably Dpp6 but
    # not variant-identifiable. Simulating one variant at a time makes the answer unambiguous: every read
    # has a known true origin, so any mass landing on another variant is leakage, and a column that is
    # heavy for EVERY row is a sink. The equimolar mix cannot separate "DPP6S is over-called" from
    # "the others are under-called"; this can.
    conf = []
    if a.mode in ('purity', 'both'):
        by_class = {}
        for t, (v, s) in seqs.items():
            by_class.setdefault(v, []).append(t)
        for layout in layouts:
            for tau in [t for t in taus if np.isinf(t) or t <= 1000]:
                for v in sorted(by_class):
                    if VAR_PANEL[v] not in panels:
                        continue
                    sub = {t: seqs[t] for t in by_class[v]}
                    tag = f'pure_{v}_{layout}_tau{"inf" if np.isinf(tau) else int(tau)}'
                    d = os.path.join(a.out, tag); os.makedirs(d, exist_ok=True)
                    r1, r2, truth = simulate(sub, a.mol, a.readlen, tau, layout, rng)
                    if not r1:
                        continue
                    r1p = os.path.join(d, 'R1.fq.gz'); write_fq(r1p, r1)
                    r2p = None
                    if layout == 'pe':
                        r2p = os.path.join(d, 'R2.fq.gz'); write_fq(r2p, r2)
                    sf = run_salmon(a.index, r1p, r2p, os.path.join(d, 'quant'), a.threads)
                    os.remove(r1p)
                    if r2p:
                        os.remove(r2p)
                    if not sf:
                        continue
                    got = read_quant(sf, tx2v)
                    panel = VAR_PANEL[v]
                    vs = [x for x in PANEL_V[panel] if x in set(w for w, _ in seqs.values())]
                    tot = sum(got.get(x, (0, 0))[1] for x in vs) or 1.0
                    for x in vs:
                        conf.append(dict(layout=layout, tau=('inf' if np.isinf(tau) else int(tau)),
                                         panel=panel, true_variant=v, assigned_variant=x,
                                         pct=got.get(x, (0.0, 0.0))[1] / tot * 100))
                    print(f'{tag}: {sum(truth.values()):,} reads -> '
                          f'{got.get(v,(0,0))[1]/tot*100:.1f}% stayed on {v}')
    if a.mode == 'purity':
        _emit_purity(conf, a.out)
        return
    for layout in layouts:
        for tau in taus:
            tag = f'{layout}_tau{"inf" if np.isinf(tau) else int(tau)}'
            d = os.path.join(a.out, tag); os.makedirs(d, exist_ok=True)
            r1, r2, truth = simulate(seqs, a.mol, a.readlen, tau, layout, rng)
            if not r1:
                print(f'{tag}: no reads simulated'); continue
            r1p = os.path.join(d, 'R1.fq.gz'); write_fq(r1p, r1)
            r2p = None
            if layout == 'pe':
                r2p = os.path.join(d, 'R2.fq.gz'); write_fq(r2p, r2)
            sf = run_salmon(a.index, r1p, r2p, os.path.join(d, 'quant'), a.threads)
            if not sf:
                print(f'{tag}: salmon failed'); continue
            got = read_quant(sf, tx2v)
            for panel in panels:
                vs = [v for v in PANEL_V[panel] if v in truth]
                if not vs:
                    continue
                tpm_tot = sum(got.get(v, (0, 0))[0] for v in vs) or 1.0
                exp = 100.0 / len(vs)                     # equal molarity -> equal TPM share
                for v in vs:
                    obs = got.get(v, (0.0, 0.0))[0] / tpm_tot * 100
                    rows.append(dict(layout=layout, tau=('inf' if np.isinf(tau) else int(tau)),
                                     panel=panel, variant=v, expected_pct=exp, observed_pct=obs,
                                     ratio=obs / exp if exp else np.nan,
                                     sim_reads=truth.get(v, 0)))
            os.remove(r1p)
            if r2p:
                os.remove(r2p)
            print(f'{tag}: done ({sum(truth.values()):,} reads)')

    import pandas as pd
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(a.out, 'variant_calibration.csv'), index=False)
    for panel in panels:
        P = R[R.panel == panel]
        if P.empty:
            continue
        print('\n' + '=' * 88)
        print(f'{panel}: observed % of panel TPM, equal molarity in (expected '
              f'{100/P.variant.nunique():.1f}% each)')
        print('=' * 88)
        piv = P.pivot_table(index='variant', columns=['layout', 'tau'], values='observed_pct')
        cols = sorted(piv.columns, key=lambda c: (c[0], 1e9 if c[1] == 'inf' else -int(c[1])))
        print(piv[cols].round(1).to_string())
    _emit_purity(conf, a.out)
    print(f'\nwrote {a.out}/variant_calibration.csv')
    print('\nA flat column means unbiased assignment at that coverage profile. Watch DPP10a across the tau')
    print('sweep: if it is recovered at tau=inf and collapses as tau shrinks, the 3\'-bias mechanism is')
    print('confirmed and the ratio column is the correction factor for the polyA bulk data.')


if __name__ == '__main__':
    main()
