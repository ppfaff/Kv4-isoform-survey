#!/usr/bin/env python3
"""Extract the Kv4-complex gene panel for thalamic neurons from the ABC Atlas 10x whole-mouse-brain data.

Cells are taken from the thalamus dissections (WMB-10Xv3-TH, WMB-10Xv2-TH) plus, optionally, any other
dissection that contributes cells to the thalamic taxonomy classes (TH Glut, MH-LH Glut, the RT-ZI and ZI
GABA subclasses of HY GABA) so that reticular-nucleus neurons dissected with hypothalamus are not lost.

Downloads go through AbcProjectCache (public S3 bucket, no credentials). The two TH log2 matrices are
~4 GB each; a full extra-matrix sweep (HY, MB) adds ~10 GB. Run on a machine with a fast link.

Usage:
  python3 extract_th_10x.py --cache /path/to/abc_atlas_cache --out tables [--extra-matrices] [--release 20250531]

Writes:
  OUT/th10x_cells.parquet        per-cell: metadata + panel genes (log2 CPM+1, 10x 'log2' matrices)
  OUT/th10x_cell_counts.csv      cells per (dataset, feature matrix, class, subclass)
  OUT/th10x_gene_panel.csv       gene -> Ensembl id used
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import anndata
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from th_common import PANEL, TH_ZI_CLASSES, write_cells  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument('--cache', required=True)
ap.add_argument('--out', default='tables')
ap.add_argument('--release', default='20250531', help='ABC Atlas manifest release')
ap.add_argument('--extra-matrices', action='store_true',
                help='also read non-TH matrices that hold cells of thalamic classes (RT-ZI in HY etc.)')
ap.add_argument('--matrices', nargs='*', default=None, help='explicit feature_matrix_labels to read')
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)

cache = AbcProjectCache.from_cache_dir(args.cache)
cache.load_manifest(f'releases/{args.release}/manifest.json')
print('manifest', cache.current_manifest)

cell = cache.get_metadata_dataframe(directory='WMB-10X', file_name='cell_metadata_with_cluster_annotation',
                                    dtype={'cell_label': str})
cell.set_index('cell_label', inplace=True)
print(f'{len(cell):,} cells in WMB-10X')

is_th_dissection = cell['feature_matrix_label'].str.endswith('-TH')
is_th_class = cell['class'].isin(TH_ZI_CLASSES) & ~cell['subclass'].str.contains(
    r'PVHd|DMH|AHN|ARH|TU-|SBPV|PVpo|LHA|MPN|Gsc|PVa|Tbx3|Onecut3', regex=True)
keep = is_th_dissection | (is_th_class if args.extra_matrices else (is_th_class & is_th_dissection))
cells = cell[keep].copy()
if args.matrices:
    cells = cells[cells['feature_matrix_label'].isin(args.matrices)]
print(f'{len(cells):,} cells selected; matrices: {sorted(cells.feature_matrix_label.unique())}')

gene = cache.get_metadata_dataframe(directory='WMB-10X', file_name='gene').set_index('gene_identifier')
gsel = gene[gene['gene_symbol'].isin(PANEL)]
missing = sorted(set(PANEL) - set(gsel['gene_symbol']))
if missing:
    print('WARNING genes not in ABC gene table:', missing)
gsel[['gene_symbol']].to_csv(f'{args.out}/th10x_gene_panel.csv')

frames = []
for fml, sub in cells.groupby('feature_matrix_label'):
    ds = sub['dataset_label'].iloc[0]
    t0 = time.time()
    path = cache.get_file_path(directory=ds, file_name=f'{fml}/log2')
    ad = anndata.read_h5ad(path, backed='r')
    idx = ad.obs_names.intersection(sub.index)
    expr = ad[idx, gsel.index].to_df()
    expr.columns = gsel.loc[expr.columns, 'gene_symbol']
    ad.file.close()
    frames.append(expr)
    print(f'  {fml}: {len(idx):,} cells, {time.time()-t0:.0f}s')

expr = pd.concat(frames)
meta_cols = ['dataset_label', 'feature_matrix_label', 'library_method', 'region_of_interest_acronym',
             'donor_label', 'donor_sex', 'donor_genotype', 'neurotransmitter', 'class', 'subclass',
             'supertype', 'cluster', 'cluster_alias']
meta_cols = [c for c in meta_cols if c in cells.columns]
out = cells.loc[expr.index, meta_cols].join(expr)
path = write_cells(out, f'{args.out}/th10x_cells')
out.groupby(['dataset_label', 'feature_matrix_label', 'class', 'subclass']).size().rename('n_cells') \
   .reset_index().to_csv(f'{args.out}/th10x_cell_counts.csv', index=False)
print(f'wrote {path}  ({len(out):,} cells x {len(expr.columns)} genes)')
