"""Shared constants for the thalamus (ABC Atlas / THALMANAC) extension of the Kv4 isoform survey.

Gene panel, thalamocortical-loop grouping of CCFv3 nuclei, and the cluster -> nucleus lookup that
THALMANAC ships (abc_merfish_analysis/resources/annotations_c2n_combined.csv).
"""
import os
import pandas as pd

# ---------------------------------------------------------------- gene panel
KV4 = ['Kcnd1', 'Kcnd2', 'Kcnd3']
KCHIP = ['Kcnip1', 'Kcnip2', 'Kcnip3', 'Kcnip4']
DPP = ['Dpp6', 'Dpp10']
COMPLEX = KV4 + KCHIP + DPP
# burst / rebound machinery of the thalamocortical loop, for context next to the A-current genes
BURST = ['Cacna1g', 'Cacna1h', 'Cacna1i', 'Hcn1', 'Hcn2', 'Hcn3', 'Hcn4',
         'Kcnc1', 'Kcnc2', 'Kcnc3', 'Kcna1', 'Kcna2', 'Kcnq2', 'Kcnq3', 'Scn1a', 'Scn8a',
         'Kcnn2', 'Kcnn3', 'Kcnk1', 'Kcnk2', 'Cacna1a', 'Cacna1b', 'Cacna1e', 'Gabra1', 'Gabra4', 'Gabrd', 'Gabrb3']
# identity markers: relay vs TRN vs interneuron; TRN subtypes; cortical layers (for the cortex-side join)
MARKERS = ['Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Pvalb', 'Sst', 'Calb1', 'Calb2', 'Tcf7l2', 'Prkcd', 'Grin2c',
           'Gnb3', 'Spp1', 'Ecel1', 'Ntsr1', 'Rorb', 'Foxp2', 'Fezf2', 'Syt6', 'Snap25', 'Rbfox3', 'Gfap', 'Mbp']
PANEL = COMPLEX + BURST + MARKERS

# ----------------------------------------------- thalamocortical loop grouping of CCFv3 nuclei
LOOP_GROUPS = {
    'TRN (RT)': ['RT'],
    'first-order sensory relay': ['VPM', 'VPL', 'VPMpc', 'VPLpc', 'LGd', 'MG'],
    'higher-order relay': ['PO', 'LP', 'POL', 'PoT', 'SGN', 'MD', 'Eth'],
    'motor relay': ['VAL', 'VM'],
    'intralaminar / midline': ['CM', 'CL', 'PCN', 'PF', 'IMD', 'IAD', 'IAM', 'RH', 'RE', 'Xi', 'SMT', 'PR'],
    'anterior / limbic': ['AD', 'AV', 'AM', 'LD', 'PVT', 'PT'],
    'epithalamus': ['MH', 'LH'],
    'ZI / LGv / IGL (non-relay)': ['ZI', 'LGv', 'IGL', 'SubG', 'IntG', 'SPA', 'SPFm', 'SPFp', 'PIL', 'PP'],
}
LOOP_ORDER = list(LOOP_GROUPS)
NUCLEUS_TO_GROUP = {n: g for g, ns in LOOP_GROUPS.items() for n in ns}

# taxonomy classes that THALMANAC treats as thalamic (TH+ZI) neurons, by CCN20230722 class name
TH_ZI_CLASSES = ['18 TH Glut', '17 MH-LH Glut', '12 HY GABA']
RT_SUBCLASS = '093 RT-ZI Gnb3 Gaba'


def cluster_to_nucleus():
    """cluster label -> nucleus (string; multi-nucleus entries like 'LP SGN' kept as given), from THALMANAC."""
    try:
        from importlib.resources import files
        p = files('abc_merfish_analysis') / 'resources' / 'annotations_c2n_combined.csv'
        df = pd.read_csv(p)
    except Exception:
        here = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(os.path.join(here, 'annotations_c2n_combined.csv'))
    df = df.dropna(subset=['nuclei'])
    # a cluster may be annotated to more than one nucleus; keep all, first listed = primary
    prim = df.groupby('cluster')['nuclei'].first()
    return prim.to_dict(), df


def loop_group_of_nucleus(nuc):
    """Map a nucleus string (possibly 'LP SGN' or 'PVT IMD') to a loop group using its first token."""
    if not isinstance(nuc, str) or not nuc:
        return 'unassigned'
    first = nuc.split()[0]
    return NUCLEUS_TO_GROUP.get(first, 'other')


def subclass_fallback_group(subclass):
    """Loop group guessed from the subclass name when a cluster has no nucleus annotation."""
    if not isinstance(subclass, str):
        return 'unassigned'
    s = subclass
    if 'RT-ZI' in s:
        return 'TRN (RT)'
    if 'ZI ' in s or 'LGv' in s or 'IGL' in s or 'SPA' in s or 'SPF' in s:
        return 'ZI / LGv / IGL (non-relay)'
    if 'MH ' in s or 'LH ' in s:
        return 'epithalamus'
    if 'AD ' in s or 'AV ' in s or 'PVT' in s:
        return 'anterior / limbic'
    if 'CM-' in s or 'PF ' in s or 'RE-' in s:
        return 'intralaminar / midline'
    if 'MG-' in s or 'POL' in s:
        return 'higher-order relay'
    if 'TH Prkcd' in s:
        return 'relay (Prkcd Grin2c, nucleus unresolved)'
    return 'other'


def write_cells(df, base):
    """Write a per-cell table as parquet if an engine is available, else csv.gz. Returns the path."""
    try:
        df.to_parquet(base + '.parquet')
        return base + '.parquet'
    except ImportError:
        df.to_csv(base + '.csv.gz')
        return base + '.csv.gz'


def read_cells(base):
    if os.path.exists(base + '.parquet'):
        return pd.read_parquet(base + '.parquet')
    return pd.read_csv(base + '.csv.gz', index_col=0)
