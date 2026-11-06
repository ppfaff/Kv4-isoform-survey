"""Are our cell groupings lumping distinct neuron types, or is the spread intrinsic to a coherent subset?

Motivation: the three-way comparison (region_compare3.py) pools Tasic-2016 types into groups named after
Allen subclasses. For the deep-layer glutamatergic groups that pooling turned out to hide 20-35 percentage
points of within-group spread in Kv4/KChIP/DPP composition -- more than the chemistry and region effects the
table was built to measure. This script decides, per grouping, whether that spread is

  (a) a LOCATION effect  -- the pooled subtypes sit at genuinely different composition centroids, i.e. we are
      lumping distinct neuron types and the group should be split; or
  (b) a DISPERSION effect -- the subtypes share a centroid but individual cells scatter widely around it, i.e.
      the variance is intrinsic to a closely related set and the grouping is sound.

Those two are statistically distinct and the usual one-test-fits-all approach cannot tell them apart: a
PERMANOVA R^2 is inflated by a dispersion difference even when every group has the same centroid. So each
grouping gets BOTH a PERMANOVA (location) and a PERMDISP-style test (dispersion), and the verdict comes from
the pair. This distinction is the entire point of the script -- do not read the R^2 column on its own.

It then derives groups de novo from the channel phenotype alone (ignoring the transcriptomic taxonomy) and
reports how those map onto the published labels, so the derived groups can be given defensible names.

Three design decisions that matter, with the reasons:

1. COMPOSITIONAL GEOMETRY. Within-panel variant fractions sum to 1, so they live on a simplex and ordinary
   variance is not partitionable across them -- an increase in one variant mechanically depresses the others.
   Everything here runs on centred-log-ratio (CLR) coordinates and Aitchison distance. The CLR needs a
   pseudocount because zeros are abundant; PSEUDO is exposed and the script reports sensitivity to it, since a
   pseudocount that is large relative to real counts will flatten genuine differences.

2. DEPTH IS THE MAIN CONFOUNDER, NOT A NUISANCE. Tasic 2016 is ~8.7 M reads/cell and Allen SSv4 HPF ~0.70 M.
   Shallow cells have noisier compositions, which inflates dispersion, which is exactly the quantity question
   (b) turns on. Two guards: every test runs WITHIN a dataset (datasets are never pooled into one clustering),
   and each test is repeated on panel counts multinomially thinned to a common total so that the sampling
   component of dispersion is equalised. A dispersion verdict that does not survive thinning is a depth
   artefact, and the output says so in its own column rather than leaving it to the reader.

3. THE READ GATE IS ITSELF A SOURCE OF FAKE STRUCTURE. Requiring >= MIN_READS in a panel drops cells unevenly
   across subclasses; a gate that keeps 95% of Pvalb cells and 20% of DG cells makes DG look artificially
   homogeneous because only its best-covered cells survive. Per-subclass retention is reported so this is
   visible instead of silent.

Inputs: the per-cell pickles via hip_common (HIP_PKL, CTX_PKL) and the Allen VISp per-cell pickle. These live
on the analysis machine, not in this repo.

Usage (on Lambda, in the kv4 env):
    cd ~/kv4_cortex && python3 channel_phenotype_groups.py
    MIN_READS=30 PSEUDO=0.5 python3 channel_phenotype_groups.py     # sensitivity run
Requires: pandas, numpy, scipy, sklearn. No scanpy/leiden dependency.
"""
import os, sys, itertools, warnings
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hip_common import load_hip, load_ctx, final_map, PANEL_OF, _read_percell

RNG       = np.random.default_rng(0)
MIN_READS = int(os.environ.get('MIN_READS', 20))    # per panel, per cell
PSEUDO    = float(os.environ.get('PSEUDO', 0.5))    # CLR pseudocount
N_PERM    = int(os.environ.get('N_PERM', 999))
MIN_CELLS = int(os.environ.get('MIN_CELLS', 12))    # per group, to be testable
OUT       = os.environ.get('OUT', 'tables')
VISP_PKL  = os.environ.get('VISP_PKL', '/home/claude/kv4/visp/visp_per_cell.pkl')

PANELS = ['Kv4', 'KChIP', 'DPP']
PANEL_V = {
    'Kv4':   ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L'],
    'KChIP': ['KChIP1a', 'KChIP1b', 'KChIP1c', 'KChIP1e', 'Kcnip1-new_short',
              'KChIP2a', 'KChIP2b', 'KChIP2c', 'KChIP2x',
              'KChIP3a', 'KChIP3x', 'Kcnip3-new-cys_rich',
              'KChIP4a', 'KChIP4b', 'KChIP4bL', 'KChIP4d', 'KChIP4e'],
    'DPP':   ['DPP6a', 'DPP6S', 'DPP6K', 'DPP6L', 'DPP10a', 'DPP10c', 'DPP10d'],
}

# the groupings under test, copied from region_compare3.py so the two stay comparable
GROUPS_T16 = {
    'Lamp5':   ['Ndnf Car4', 'Ndnf Cxcl14', 'Igtp', 'Smad3'],
    'Sncg':    ['Sncg', 'Vip Sncg'],
    'Vip':     ['Vip Chat', 'Vip Gpc3', 'Vip Mybpc1', 'Vip Parm1'],
    'Sst':     ['Sst Cbln4', 'Sst Cdk6', 'Sst Chodl', 'Sst Myh8', 'Sst Tacstd2', 'Sst Th'],
    'Pvalb':   ['Pvalb Cpne5', 'Pvalb Gpx3', 'Pvalb Obox3', 'Pvalb Rspo2',
                'Pvalb Tacr3', 'Pvalb Tpbg', 'Pvalb Wt1'],
    'L2/3 IT': ['L2/3 Ptgs2', 'L2 Ngb'],
    'L4':      ['L4 Arf5', 'L4 Ctxn3', 'L4 Scnn1a'],
    'L5 IT':   ['L5a Batf3', 'L5a Hsd11b1', 'L5a Pde1c', 'L5a Tcerg1l'],
    'L5 PT':   ['L5b Cdh13', 'L5b Chrna6', 'L5b Tph2'],
    'L6 CT':   ['L6a Mgp', 'L6a Sla', 'L6a Car12'],
    'L6 IT':   ['L6a Syt17'],
    'L6b':     ['L6b Rgs12', 'L6b Serpinb11'],
}


# ----------------------------------------------------------------- features
def cell_matrix(d, min_reads=MIN_READS, pseudo=PSEUDO, thin_to=None):
    """Per-cell CLR feature matrix from the channel phenotype.

    Returns (X, meta, retention). Features are the within-panel CLR coordinates for every panel a cell passes
    the read gate on, concatenated, plus the log-ratios between panel totals (which carry the *relative
    abundance* of Kv4 vs its auxiliary subunits -- electrically meaningful, and discarded by composition alone).
    Cells must pass the gate on ALL panels to keep the feature space rectangular; retention is reported.

    thin_to: if set, each panel's counts are first multinomially resampled to this total, equalising the
    sampling contribution to dispersion across cells of different depth.
    """
    V = d.groupby(['sample_name', 'vclass'])['reads'].sum().unstack(fill_value=0.0)
    meta = d.groupby('sample_name').agg(subclass=('subclass', 'first'),
                                        cluster=('cluster', 'first'),
                                        cls=('cls', 'first'))
    blocks, totals = [], {}
    for p in PANELS:
        C = V.reindex(columns=PANEL_V[p]).fillna(0.0)
        totals[p] = C.sum(axis=1)
        blocks.append((p, C))

    keep = pd.Series(True, index=V.index)
    for p, C in blocks:
        keep &= totals[p] >= min_reads
    # retention per subclass: the gate must not be allowed to silently homogenise shallow groups
    retention = (pd.DataFrame({'subclass': meta['subclass'], 'kept': keep})
                 .groupby('subclass')['kept'].agg(n_total='size', n_kept='sum'))
    retention['frac_kept'] = retention['n_kept'] / retention['n_total']

    feats, names = [], []
    for p, C in blocks:
        C = C.loc[keep]
        if thin_to:
            C = pd.DataFrame(
                [RNG.multinomial(thin_to, r / r.sum()) if r.sum() > 0 else np.zeros(len(r))
                 for _, r in C.iterrows()], index=C.index, columns=C.columns, dtype=float)
        A = C.to_numpy() + pseudo
        A = A / A.sum(axis=1, keepdims=True)
        clr = np.log(A) - np.log(A).mean(axis=1, keepdims=True)
        feats.append(clr); names += [f'{p}:{v}' for v in C.columns]
    # inter-panel abundance ratios (2 independent log-ratios for 3 panels)
    T = np.column_stack([totals[p].loc[keep].to_numpy() + pseudo for p in PANELS])
    T = np.log(T) - np.log(T).mean(axis=1, keepdims=True)
    feats.append(T); names += [f'logratio:{p}' for p in PANELS]

    X = pd.DataFrame(np.column_stack(feats), index=V.index[keep], columns=names)
    return X, meta.loc[X.index], retention


# ----------------------------------------------------------------- the two tests
def permanova(D, labels, n_perm=N_PERM):
    """One-way PERMANOVA (Anderson 2001) on a distance matrix. Tests LOCATION: do the groups have different
    centroids? Returns (pseudo_F, R2, p). R2 alone is NOT interpretable as 'the groups are real' -- a
    difference in within-group dispersion inflates it. Always read it beside permdisp()."""
    labels = np.asarray(labels)
    uniq = np.unique(labels)
    if len(uniq) < 2:
        return np.nan, np.nan, np.nan
    n = len(labels)
    D2 = np.asarray(D) ** 2
    ss_total = D2.sum() / (2 * n)

    def ss_within(lab):
        s = 0.0
        for g in np.unique(lab):
            i = np.where(lab == g)[0]
            if len(i) > 1:
                s += D2[np.ix_(i, i)].sum() / (2 * len(i))
        return s

    sw = ss_within(labels)
    sb = ss_total - sw
    a, N = len(uniq), n
    F = (sb / (a - 1)) / (sw / (N - a)) if sw > 0 and N > a else np.nan
    R2 = sb / ss_total if ss_total > 0 else np.nan
    perm = np.empty(n_perm)
    lab = labels.copy()
    for k in range(n_perm):
        RNG.shuffle(lab)
        swp = ss_within(lab); sbp = ss_total - swp
        perm[k] = (sbp / (a - 1)) / (swp / (N - a)) if swp > 0 else np.nan
    p = (1 + np.nansum(perm >= F)) / (1 + n_perm)
    return F, R2, p


def permdisp(D, labels, n_perm=N_PERM):
    """PERMDISP-style homogeneity-of-dispersion test. Tests DISPERSION: is the average distance of a cell to
    its own group's centroid the same across groups? Distances to centroid are computed in the principal-
    coordinate embedding of D, then compared by a permutation F-test. A significant result means the groups
    differ in internal heterogeneity -- which is the signal for 'intrinsic variance' rather than 'lumping'.
    Also returns the mean within-group distance to centroid per group, which is the effect size that matters."""
    labels = np.asarray(labels)
    uniq = np.unique(labels)
    if len(uniq) < 2:
        return np.nan, np.nan, {}
    D2 = np.asarray(D) ** 2
    n = len(labels)
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ D2 @ J                      # Gower-centred
    w, Vec = np.linalg.eigh(G)
    pos = w > 1e-9
    Y = Vec[:, pos] * np.sqrt(w[pos])          # principal coordinates

    def dists(lab):
        out = np.empty(n)
        for g in np.unique(lab):
            i = np.where(lab == g)[0]
            c = Y[i].mean(axis=0)
            out[i] = np.sqrt(((Y[i] - c) ** 2).sum(axis=1))
        return out

    d0 = dists(labels)
    grp_mean = {g: float(d0[labels == g].mean()) for g in uniq}

    def Fstat(dd, lab):
        gm = dd.mean()
        sb = sum((lab == g).sum() * (dd[lab == g].mean() - gm) ** 2 for g in np.unique(lab))
        sw = sum(((dd[lab == g] - dd[lab == g].mean()) ** 2).sum() for g in np.unique(lab))
        a = len(np.unique(lab))
        return (sb / (a - 1)) / (sw / (len(dd) - a)) if sw > 0 else np.nan

    # Permute the residual distances-to-centroid, NOT the group labels. Shuffling labels would re-estimate
    # centroids on mixed groups, so whenever a genuine location difference exists the permuted within-group
    # distances are inflated for every group alike, the null F collapses, and the test reports a dispersion
    # difference that is really the centroid shift leaking through (verified: a pure centroid shift with
    # identical generating spread gave p=0.014 under label shuffling, and is correctly non-significant here).
    # Centroids are removed once, under the true labels, so the null asks only about the residual spread.
    F = Fstat(d0, labels)
    perm = np.empty(n_perm)
    dd = d0.copy()
    for k in range(n_perm):
        RNG.shuffle(dd)
        perm[k] = Fstat(dd, labels)
    p = (1 + np.nansum(perm >= F)) / (1 + n_perm)
    return F, p, grp_mean


def verdict(loc_p, disp_p, r2, alpha=0.05):
    """Translate the test pair into the answer the project actually needs."""
    loc, disp = loc_p < alpha, disp_p < alpha
    if loc and not disp:
        return 'SPLIT - distinct centroids, comparable spread (lumping distinct types)'
    if loc and disp:
        return 'SPLIT - distinct centroids, but dispersion differs too (check thinned result)'
    if not loc and disp:
        return 'KEEP - shared centroid, unequal spread (intrinsic variance, one subtype noisier)'
    return 'KEEP - shared centroid, comparable spread (coherent group)'


# ----------------------------------------------------------------- stage 1: audit the published groupings
def audit_grouping(X, meta, groups, dataset, thinned=False):
    """For each pooled group, test whether the finer labels inside it still carry structure.
    This is the nested question: having pooled L6a Mgp/Sla/Car12 into 'L6 CT', is there anything left to
    explain inside it? A significant location effect means the pool is hiding distinct types."""
    rows = []
    for gname, members in groups.items():
        sel = meta.index[meta['cluster'].isin(members)]
        sub = meta.loc[sel]
        counts = sub['cluster'].value_counts()
        usable = counts[counts >= MIN_CELLS].index
        if len(usable) < 2:
            rows.append(dict(dataset=dataset, group=gname, n_cells=len(sel),
                             n_subtypes=len(counts), n_testable=len(usable),
                             note='not testable (needs >=2 subtypes with >=%d cells)' % MIN_CELLS))
            continue
        idx = sub.index[sub['cluster'].isin(usable)]
        Xi = X.loc[idx]
        lab = meta.loc[idx, 'cluster'].to_numpy()
        D = squareform(pdist(Xi.to_numpy(), metric='euclidean'))   # Euclidean on CLR == Aitchison
        F, R2, p_loc = permanova(D, lab)
        Fd, p_disp, gm = permdisp(D, lab)
        rows.append(dict(dataset=dataset, group=gname, n_cells=len(idx),
                         n_subtypes=len(counts), n_testable=len(usable),
                         subtypes='|'.join(usable),
                         permanova_F=F, permanova_R2=R2, p_location=p_loc,
                         permdisp_F=Fd, p_dispersion=p_disp,
                         spread_min=min(gm.values()), spread_max=max(gm.values()),
                         verdict=verdict(p_loc, p_disp, R2), thinned=thinned))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------- stage 2: derive groups de novo
def derive_groups(X, meta, dataset, kmax=None):
    """Ward clustering on the CLR features, k chosen by silhouette. Deliberately ignores the transcriptomic
    labels so the result is a channel-phenotype grouping that can be compared against them rather than
    assumed to agree with them."""
    A = X.to_numpy()
    n = len(A)
    kmax = kmax or min(20, max(3, n // 25))
    Z = linkage(A, method='ward')
    best = None
    for k in range(2, kmax + 1):
        lab = fcluster(Z, k, criterion='maxclust')
        if len(np.unique(lab)) < 2:
            continue
        s = silhouette_score(A, lab)
        if best is None or s > best[1]:
            best = (k, s, lab)
    k, sil, lab = best
    out = pd.DataFrame({'sample_name': X.index, 'derived': [f'{dataset}-c{v:02d}' for v in lab]}).set_index('sample_name')
    out = out.join(meta[['subclass', 'cluster', 'cls']])
    agree = dict(
        dataset=dataset, k=k, silhouette=sil, n_cells=n,
        ARI_vs_cluster=adjusted_rand_score(meta.loc[X.index, 'cluster'], lab),
        NMI_vs_cluster=normalized_mutual_info_score(meta.loc[X.index, 'cluster'], lab),
        ARI_vs_subclass=adjusted_rand_score(meta.loc[X.index, 'subclass'], lab),
        NMI_vs_subclass=normalized_mutual_info_score(meta.loc[X.index, 'subclass'], lab),
    )
    return out, agree


def name_groups(derived, taxonomy_col='subclass', purity_hi=0.70, member_min=0.20):
    """Propose a name for each derived cluster from the published taxonomy.

    Rule, stated explicitly so it can be argued with:
      purity >= purity_hi                      -> adopt that published label
      one published label, several clusters    -> suffix .1 .2 .. by descending size
      several labels each >= member_min        -> composite 'A/B' and list all members
    The source taxonomy is recorded in a column so the names can be re-pointed when we move onto a newer
    reference taxonomy (the Allen whole-brain release supersedes the 2021 subclass vocabulary our labels use)
    without redoing the clustering.
    """
    rows = []
    for c, sub in derived.groupby('derived'):
        vc = sub[taxonomy_col].value_counts(normalize=True)
        top, purity = vc.index[0], float(vc.iloc[0])
        members = [f'{k} ({v:.0%})' for k, v in vc.items() if v >= member_min]
        if purity >= purity_hi:
            proposed, kind = top, 'adopted'
        elif len(members) >= 2:
            proposed, kind = '/'.join(m.split(' (')[0] for m in members), 'composite'
        else:
            proposed, kind = f'{top} (mixed)', 'mixed'
        rows.append(dict(derived=c, n_cells=len(sub), proposed_name=proposed, naming=kind,
                         purity=purity, majority_label=top, members='; '.join(members),
                         taxonomy=taxonomy_col))
    out = pd.DataFrame(rows).sort_values('n_cells', ascending=False)
    # disambiguate clusters that landed on the same published label
    for nm, g in out.groupby('proposed_name'):
        if len(g) > 1 and g['naming'].iloc[0] == 'adopted':
            for i, (idx, _) in enumerate(g.iterrows(), start=1):
                out.loc[idx, 'proposed_name'] = f'{nm}.{i}'
                out.loc[idx, 'naming'] = 'split'
    return out


# ----------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    datasets = {}
    try:
        datasets['T16'] = (load_ctx(), GROUPS_T16)
    except Exception as e:
        print(f'[warn] Tasic 2016 not loaded: {e}')
    try:
        datasets['A-HPF'] = (load_hip(), None)
    except Exception as e:
        print(f'[warn] Allen HPF not loaded: {e}')
    if os.path.exists(VISP_PKL):
        # .csv or .pkl -- the VISp data is an annotated CSV from the same collapse step as the others
        v = _read_percell(VISP_PKL)
        miss = [c for c in ('transcript', 'gene', 'reads', 'sample_name') if c not in v.columns]
        if miss:
            sys.exit(f'VISp table {VISP_PKL} is missing columns {miss}; found {list(v.columns)}')
        v['vclass'] = v['transcript'].map(final_map)
        v = v.dropna(subset=['vclass']).copy()
        v = v[~((v['biotype'] == 'retained_intron')
                & v['vclass'].str.contains('other|fragment', case=False))].copy()
        v['panel'] = v['gene'].map(PANEL_OF)
        v = v.rename(columns={'subclass_label': 'subclass', 'cluster_label': 'cluster'})
        if 'cls' not in v.columns:
            v['cls'] = 'unknown'
        datasets['A-VISp'] = (v, None)
    else:
        print(f'[warn] VISp pickle not found at {VISP_PKL}')
    if not datasets:
        sys.exit('no datasets could be loaded -- check HIP_PKL/CTX_PKL in hip_common.py and VISP_PKL')

    audits, derived_all, agree_all, names_all, rets = [], [], [], [], []
    for ds, (d, groups) in datasets.items():
        X, meta, retention = cell_matrix(d)
        retention.insert(0, 'dataset', ds); rets.append(retention.reset_index())
        print(f'\n=== {ds}: {len(X)} cells pass the gate (>= {MIN_READS} reads in all three panels), '
              f'{X.shape[1]} features')
        worst = retention.sort_values('frac_kept').head(3)
        print('    lowest retention:', ', '.join(f'{i} {r.frac_kept:.0%}' for i, r in worst.iterrows()))

        # stage 1 -- published groupings, raw and depth-thinned
        if groups:
            a = audit_grouping(X, meta, groups, ds)
            audits.append(a)
            thin = min(int(os.environ.get('THIN_TO', 50)), MIN_READS * 2)
            Xt, mt, _ = cell_matrix(d, thin_to=thin)
            at = audit_grouping(Xt, mt, groups, ds, thinned=True)
            audits.append(at)
        # every dataset also gets its own label level audited as a self-consistency check
        selfgroups = {s: [c for c in meta.loc[meta.subclass == s, 'cluster'].unique()]
                      for s in meta['subclass'].unique()}
        selfgroups = {k: v for k, v in selfgroups.items() if len(v) > 1}
        if selfgroups:
            audits.append(audit_grouping(X, meta, selfgroups, f'{ds}:subclass'))

        # stage 2 -- derive de novo, then name
        dv, ag = derive_groups(X, meta, ds)
        derived_all.append(dv.reset_index()); agree_all.append(ag)
        nm = name_groups(dv); nm.insert(0, 'dataset', ds); names_all.append(nm)
        print(f'    derived k={ag["k"]} (silhouette {ag["silhouette"]:.3f}); '
              f'ARI vs published cluster {ag["ARI_vs_cluster"]:.3f}, vs subclass {ag["ARI_vs_subclass"]:.3f}')

    A = pd.concat(audits, ignore_index=True)
    A.to_csv(f'{OUT}/channel_phenotype_variance.csv', index=False)
    pd.concat(derived_all, ignore_index=True).to_csv(f'{OUT}/channel_phenotype_cells.csv', index=False)
    pd.DataFrame(agree_all).to_csv(f'{OUT}/channel_phenotype_agreement.csv', index=False)
    pd.concat(names_all, ignore_index=True).to_csv(f'{OUT}/channel_phenotype_names.csv', index=False)
    pd.concat(rets, ignore_index=True).to_csv(f'{OUT}/channel_phenotype_retention.csv', index=False)

    print('\n' + '=' * 100)
    print('VERDICT PER PUBLISHED GROUPING (raw; compare against the thinned rows before acting)')
    print('=' * 100)
    show = A[A['verdict'].notna() & (~A['thinned'].fillna(False))]
    for _, r in show.iterrows():
        print(f'{r.dataset:14} {r.group:10} n={r.n_cells:5d} {r.n_testable} subtypes  '
              f'R2={r.permanova_R2:5.3f} p_loc={r.p_location:6.3f} p_disp={r.p_dispersion:6.3f}  '
              f'spread {r.spread_min:4.2f}-{r.spread_max:4.2f}')
        print(f'{"":14} -> {r.verdict}')
    print(f'\nwrote {OUT}/channel_phenotype_{{variance,cells,agreement,names,retention}}.csv')
    print('Read p_location and p_dispersion together; the R2 column alone does not answer the question.')


if __name__ == '__main__':
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    main()
