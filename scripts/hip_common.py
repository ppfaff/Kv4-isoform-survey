"""Shared definitions for the hippocampal (Allen SSv4, Yao 2021) Kv4-complex isoform analysis.
Variant nomenclature and colour/hatch scheme are identical to the cortex (Tasic 2016) work."""
import pandas as pd, numpy as np

HIP_PKL = '/home/claude/kv4/hip/hip_per_cell.pkl'

final_map = {
 'Dpp6-201':'DPP6a','Dpp6-206':'DPP6a','Dpp6-208':'DPP6a',
 'Dpp6-202':'DPP6S','Dpp6-207':'DPP6S',
 'Dpp6-203':'DPP6K',
 'Dpp6-204':'DPP6L','Dpp6-209':'DPP6L',
 'Dpp6-205':'Dpp6 other/fragment',
 'Dpp10-201':'DPP10d','Dpp10-203':'DPP10d',
 'Dpp10-202':'DPP10c',
 'Dpp10-206':'DPP10a','Dpp10-205':'DPP10a',
 'Dpp10-204':'Dpp10 other/fragment',
 'Kcnd3-202':'Kv4.3S','Kcnd3-201':'Kv4.3L','Kcnd3-203':'Kv4.3L',
 'Kcnd1-201':'Kv4.1','Kcnd2-201':'Kv4.2',
 'Kcnip1-201':'KChIP1b','Kcnip1-203':'KChIP1a','Kcnip1-207':'KChIP1c','Kcnip1-206':'KChIP1e',
 'Kcnip1-202':'Kcnip1-new_short','Kcnip1-204':'Kcnip1 other/fragment','Kcnip1-205':'Kcnip1 other/fragment',
 'Kcnip2-208':'KChIP2c','Kcnip2-210':'KChIP2c',
 'Kcnip2-201':'KChIP2b','Kcnip2-202':'KChIP2b','Kcnip2-204':'KChIP2b',
 'Kcnip2-209':'KChIP2a','Kcnip2-211':'KChIP2x',
 'Kcnip2-203':'Kcnip2 other/fragment','Kcnip2-205':'Kcnip2 other/fragment',
 'Kcnip2-206':'Kcnip2 other/fragment','Kcnip2-207':'Kcnip2 other/fragment',
 'Kcnip3-203':'KChIP3a','Kcnip3-202':'KChIP3x','Kcnip3-201':'Kcnip3-new-cys_rich',
 'Kcnip3-205':'Kcnip3 other/fragment','Kcnip3-204':'Kcnip3 other/fragment',
 'Kcnip4-201':'KChIP4bL','Kcnip4-207':'KChIP4b','Kcnip4-208':'KChIP4a','Kcnip4-204':'KChIP4e','Kcnip4-206':'KChIP4d',
 'Kcnip4-202':'Kcnip4 other/noncoding','Kcnip4-203':'Kcnip4 other/noncoding','Kcnip4-205':'Kcnip4 other/noncoding',
}
PANEL_OF = {'Kcnd1':'Kv4','Kcnd2':'Kv4','Kcnd3':'Kv4',
            'Kcnip1':'KChIP','Kcnip2':'KChIP','Kcnip3':'KChIP','Kcnip4':'KChIP',
            'Dpp6':'DPP','Dpp10':'DPP'}

GABA = ['Lamp5','Sncg','Vip','Sst','Sst Chodl','Pvalb']
NONN = ['Micro-PVM','Astro','Oligo']
# principal-cell / glutamatergic subclasses in anatomical-ish order: hippocampus proper -> subiculum -> adjoining cortex
GLUT = ['DG','CA3','CA2-IG-FC','CA1-ProS','SUB-ProS','NP SUB','CT SUB',
        'L2/3 IT RHP','L2/3 IT PPP','L2/3 IT ENTl','L4 RSP-ACA','L4/5 IT CTX','L5 PPP','NP PPP','L6 IT CTX','L6 CT CTX','L6b CTX']

def cls_of(sub):
    if sub in GABA: return 'GABA'
    if sub in NONN: return 'NonN'
    return 'Glut'

SUBCLASS_ORDER = GABA + GLUT + NONN

HEX = {
 'purple':'#7B4FA3','light_purple':'#B9A3D9','sky_blue':'#4FA8D8','bluish_green':'#4CAF50',
 'orange':'#E8822D','yellow':'#E8D44D',
 'terracotta':'#C1593C','mint_green':'#8FD9B6','charcoal':'#3A3A3A','maize':'#F0B429',
 'white':'#FFFFFF','blue_dot':'#2E6FBA','light_gray':'#CFCFCF','red':'#D6372E',
}
GENEBASE = {'Kcnip1':HEX['purple'],'Kcnip2':HEX['sky_blue'],'Kcnip3':HEX['bluish_green'],'Kcnip4':HEX['orange'],
            'Dpp6':HEX['purple'],'Dpp10':HEX['orange'],'Kcnd1':HEX['sky_blue'],'Kcnd3':HEX['purple']}
VSPEC = {
 'KChIP1a': (HEX['purple'], 'X', HEX['yellow']),
 'KChIP1b': (HEX['purple'], 'X', HEX['mint_green']),
 'KChIP1c': (HEX['purple'], 'X', HEX['white']),
 'KChIP1e': (HEX['purple'], 'bars', HEX['maize']),
 'Kcnip1-new_short': (HEX['purple'], 'solid', HEX['purple']),
 'Kcnip1 other/fragment': (GENEBASE['Kcnip1'], 'O', HEX['yellow']),
 'KChIP2a': (HEX['sky_blue'], 'X', HEX['terracotta']),
 'KChIP2b': (HEX['sky_blue'], 'X', HEX['charcoal']),
 'KChIP2c': (HEX['sky_blue'], 'solid', HEX['sky_blue']),
 'KChIP2x': (HEX['sky_blue'], 'bars', HEX['white']),
 'Kcnip2 other/fragment': (GENEBASE['Kcnip2'], 'O', HEX['yellow']),
 'KChIP3a': (HEX['bluish_green'], 'X', HEX['terracotta']),
 'KChIP3x': (HEX['bluish_green'], 'bars', HEX['white']),
 'Kcnip3-new-cys_rich': (HEX['bluish_green'], 'X', HEX['charcoal']),
 'Kcnip3 other/fragment': (GENEBASE['Kcnip3'], 'O', HEX['yellow']),
 'KChIP4a': (HEX['orange'], 'bars', HEX['white']),
 'KChIP4b': (HEX['orange'], 'solid', HEX['orange']),
 'KChIP4bL': (HEX['orange'], 'X', HEX['charcoal']),
 'KChIP4d': (HEX['orange'], 'X', HEX['terracotta']),
 'KChIP4e': (HEX['orange'], 'bars', HEX['maize']),
 'Kcnip4 other/noncoding': (GENEBASE['Kcnip4'], 'O', HEX['yellow']),
 'DPP6a': (HEX['purple'], 'dots', HEX['maize']),
 'DPP6S': (HEX['purple'], 'solid', HEX['purple']),
 'DPP6L': (HEX['purple'], 'X', HEX['white']),
 'DPP6K': (HEX['purple'], 'X', HEX['mint_green']),
 'Dpp6 other/fragment': (GENEBASE['Dpp6'], 'O', HEX['yellow']),
 'DPP10a': (HEX['orange'], 'dots', HEX['blue_dot']),
 'DPP10c': (HEX['orange'], 'solid', HEX['orange']),
 'DPP10d': (HEX['orange'], 'X', HEX['charcoal']),
 'Dpp10 other/fragment': (GENEBASE['Dpp10'], 'O', HEX['yellow']),
 'Kv4.1': (HEX['sky_blue'], 'solid', HEX['sky_blue']),
 'Kv4.2': (HEX['orange'], 'solid', HEX['orange']),
 'Kv4.3S': (HEX['light_purple'], 'solid', HEX['light_purple']),
 'Kv4.3L': (HEX['purple'], 'solid', HEX['purple']),
}
HATCH = {'bars':'---','dots':'..','X':'xx','O':'OO','solid':None}
CLS_COLOR = {'GABA':'#7A3FA0','Glut':'#2060B0','NonN':'#777777'}
BLANK = '__blank__'
LEGEND_GROUPS = {
 'Kv4': [['Kv4.1','Kv4.2','Kv4.3S','Kv4.3L']],
 'KChIP': [
   ['KChIP1a','KChIP1b','KChIP1c','KChIP1e','Kcnip1-new_short','Kcnip1 other/fragment', BLANK,
    'KChIP2a','KChIP2b','KChIP2c','KChIP2x','Kcnip2 other/fragment'],
   ['KChIP3a','KChIP3x','Kcnip3-new-cys_rich','Kcnip3 other/fragment', BLANK,
    'KChIP4a','KChIP4b','KChIP4bL','KChIP4d','KChIP4e','Kcnip4 other/noncoding'],
 ],
 'DPP': [
   ['DPP6a','DPP6K','DPP6L','DPP6S','Dpp6 other/fragment'],
   ['DPP10a','DPP10c','DPP10d','Dpp10 other/fragment'],
 ],
}


def load_hip():
    """Per-transcript-per-cell table with variant classes applied, same filtering as the cortex work:
    unmapped transcripts (no literature support) dropped; retained_intron transcripts dropped only when
    they fall in an unnamed other/fragment bucket."""
    d = pd.read_pickle(HIP_PKL)
    d['vclass'] = d['transcript'].map(final_map)
    d = d.dropna(subset=['vclass']).copy()
    is_other = d['vclass'].str.contains('other|fragment', case=False, regex=True)
    d = d[~((d['biotype'] == 'retained_intron') & is_other)].copy()
    d['panel'] = d['gene'].map(PANEL_OF)
    d['cls'] = d['subclass_label'].map(cls_of)
    d = d.rename(columns={'subclass_label': 'subclass', 'cluster_label': 'cluster', 'region_label': 'region'})
    return d


# ---------------------------------------------------------------- cortex (Tasic 2016, full 1,677 cells)
CTX_PKL = '/home/claude/kv4/ctx/ctx_per_cell.pkl'
CTX_CLS = {'GABA-ergic Neuron': 'GABA', 'Glutamatergic Neuron': 'Glut'}

def load_ctx():
    """Same shape as load_hip(): columns subclass (=Tasic primary type), cluster (=same), cls, region, vclass, panel."""
    d = pd.read_pickle(CTX_PKL)
    d['vclass'] = d['transcript'].map(final_map)
    d = d.dropna(subset=['vclass']).copy()
    is_other = d['vclass'].str.contains('other|fragment', case=False, regex=True)
    d = d[~((d['biotype'] == 'retained_intron') & is_other)].copy()
    d['panel'] = d['gene'].map(PANEL_OF)
    d['cls'] = d['broad_type'].map(CTX_CLS).fillna('NonN')
    d = d.rename(columns={'primary_type_label': 'subclass', 'dissection': 'region'})
    d['cluster'] = d['subclass']
    return d

def ctx_order(d):
    """Tasic types in taxonomy order (primary_type_id), GABA first, then Glut, then non-neuronal."""
    t = d[['subclass', 'primary_type_id', 'cls']].drop_duplicates().sort_values('primary_type_id')
    return [s for c in ['GABA', 'Glut', 'NonN'] for s in t[t.cls == c]['subclass']]

def load_dataset(name):
    if name == 'hip':
        d = load_hip(); return d, SUBCLASS_ORDER, 'Hippocampal formation', 'Allen SSv4 (Yao 2021)', lambda s: cls_of(s)
    d = load_ctx(); order = ctx_order(d); cmap = dict(d[['subclass', 'cls']].drop_duplicates().values)
    return d, order, 'Primary visual cortex', 'Tasic 2016 (GSE71585)', lambda s: cmap[s]
