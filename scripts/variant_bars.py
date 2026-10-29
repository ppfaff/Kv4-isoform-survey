import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

d = pd.read_pickle('/home/claude/kv4/iso/d.pkl')

final_map={
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
 # Kcnd1-202 (partial/no-CDS fragment) and Kcnd3-204/205 (no literature support) dropped entirely
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

d['vclass']=d['transcript'].map(final_map)
d=d.dropna(subset=['vclass'])

# drop retained_intron transcripts ONLY when they fall in an unnamed "other/fragment" bucket —
# never drop a transcript that carries a confirmed named variant identity (e.g. Dpp6-206/
# Dpp6-209/Kcnip2-204 are GENCODE-tagged retained_intron but were sequence-confirmed as
# real DPP6a/DPP6L/KChIP2b transcripts, so biotype is overridden by the confirmed mapping)
is_other = d['vclass'].str.contains('other|fragment', case=False, regex=True)
drop_mask = (d['biotype']=='retained_intron') & is_other
dropped = sorted(d.loc[drop_mask, 'transcript'].unique())
d = d[~drop_mask]
print('dropped retained_intron "other" transcripts:', dropped)

classmap={'Kcnd1':'Kv4','Kcnd2':'Kv4','Kcnd3':'Kv4',
          'Kcnip1':'KChIP','Kcnip2':'KChIP','Kcnip3':'KChIP','Kcnip4':'KChIP',
          'Dpp6':'DPP','Dpp10':'DPP'}
d['panel']=d['gene'].map(classmap)

agg = d.groupby(['type','primary_type_id','cls','panel','vclass'])['reads'].sum().reset_index()
agg['panel_total']=agg.groupby(['type','panel'])['reads'].transform('sum')
agg['frac']=agg['reads']/agg['panel_total'].replace(0,np.nan)

# ---------------- color scheme ----------------
HEX = {
 'purple':'#7B4FA3','light_purple':'#B9A3D9','sky_blue':'#4FA8D8','bluish_green':'#4CAF50',
 'orange':'#E8822D','yellow':'#E8D44D',
 'terracotta':'#C1593C','mint_green':'#8FD9B6','charcoal':'#3A3A3A','maize':'#F0B429',
 'white':'#FFFFFF','blue_dot':'#2E6FBA','light_gray':'#CFCFCF','red':'#D6372E',
}

GENEBASE = {'Kcnip1':HEX['purple'],'Kcnip2':HEX['sky_blue'],'Kcnip3':HEX['bluish_green'],'Kcnip4':HEX['orange'],
            'Dpp6':HEX['purple'],'Dpp10':HEX['orange'],'Kcnd1':HEX['sky_blue'],'Kcnd3':HEX['purple']}

# VSPEC[v] = (base_color, kind, pattern_color)
# kind: 'solid' | 'bars' (tiled horiz-bar hatch) | 'dots' (tiled dot hatch) | 'X' (single centered X) | 'O' (single centered yellow dot)
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

 'DPP10a': (HEX['orange'], 'dots', HEX['blue_dot']),   # changed: blue dot shade per feedback
 'DPP10c': (HEX['orange'], 'solid', HEX['orange']),
 'DPP10d': (HEX['orange'], 'X', HEX['charcoal']),
 'Dpp10 other/fragment': (GENEBASE['Dpp10'], 'O', HEX['yellow']),

 'Kv4.1': (HEX['sky_blue'], 'solid', HEX['sky_blue']),
 'Kv4.2': (HEX['orange'], 'solid', HEX['orange']),
 'Kv4.3S': (HEX['light_purple'], 'solid', HEX['light_purple']),
 'Kv4.3L': (HEX['purple'], 'solid', HEX['purple']),
}

PANELS = ['Kv4','KChIP','DPP']

# n cells per type that have made it through the pipeline so far (any of the 9 genes quantified)
ncells = d.groupby('type')['sample_name'].nunique()

order = agg[['type','primary_type_id','cls']].drop_duplicates().sort_values('primary_type_id')
types = order['type'].tolist()
cls_of = dict(zip(order['type'], order['cls']))
CLS_COLOR = {'GABA':'#7A3FA0','Glut':'#2060B0','NonN':'#777777'}

MIN_READS = 20
BAR_H = 0.75
GAP = 1.2  # extra vertical space inserted between GABA / Glut / NonN blocks

# assign each type a y-position from the bottom, inserting a gap whenever cls changes
cls_order = [cls_of[t] for t in types]
ypos = []
y = 0.0
for i, t in enumerate(types):
    if i > 0 and cls_order[i] != cls_order[i-1]:
        y += GAP
    ypos.append(y)
    y += 1.0
ypos = {t: yp for t, yp in zip(types, ypos)}
# plotted bottom-to-top, but types list is in top-to-bottom reading order -> flip
ymax = max(ypos.values())
ypos = {t: ymax - yp for t, yp in ypos.items()}

GLOBAL_ORDER = {}
for panel in PANELS:
    tot = agg[agg['panel']==panel].groupby('vclass')['reads'].sum().sort_values(ascending=False)
    GLOBAL_ORDER[panel] = [v for v in tot.index if v in VSPEC]

fig, axes = plt.subplots(1, 4, figsize=(19, 0.17*(ymax+2)+3.2), sharey=True,
                          gridspec_kw={'width_ratios':[3,3,3,1.1]})
bar_axes, ncell_ax = axes[:3], axes[3]

for ax, panel in zip(bar_axes, PANELS):
    sub = agg[agg['panel']==panel]
    for yi, t in enumerate(types):
        row = sub[sub['type']==t]
        total = row['reads'].sum()
        y = ypos[t]
        if total < MIN_READS or row.empty:
            ax.barh(y, 1.0, color='#EEEEEE', height=BAR_H, zorder=1)
            continue
        row_sorted = row[row['frac'] > 0].sort_values('frac', ascending=False)
        left = 0.0
        for _, r in row_sorted.iterrows():
            v, frac = r['vclass'], r['frac']
            base, kind, pcolor = VSPEC[v]
            hatch = {'bars':'---','dots':'..','X':'xx','O':'OO','solid':None}[kind]
            ax.barh(y, frac, left=left, color=base, edgecolor=pcolor, hatch=hatch, height=BAR_H, linewidth=0.3, zorder=2)
            left += frac
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.75, ymax + 0.75)
    ax.set_title(panel, fontsize=12, fontweight='bold')
    ax.set_yticks([ypos[t] for t in types])
    ax.set_xticks([])
    for spine in ['top','right','left','bottom']:
        ax.spines[spine].set_visible(False)

axes[0].set_yticklabels(types, fontsize=6)
for lbl, t in zip(axes[0].get_yticklabels(), types):
    lbl.set_color(CLS_COLOR[cls_of[t]])

# ---------------- 4th column: n cells quantified so far, per type ----------------
nmax = max(1, int(ncells.reindex(types).fillna(0).max()))
for t in types:
    n = int(ncells.get(t, 0))
    y = ypos[t]
    color = '#999999' if n < 3 else '#3A6EA5'
    ncell_ax.barh(y, n, color=color, height=BAR_H, zorder=2)
    ncell_ax.text(n + nmax*0.04, y, str(n), va='center', ha='left', fontsize=5.5)
ncell_ax.set_xlim(0, nmax*1.35)
ncell_ax.set_ylim(-0.75, ymax + 0.75)
ncell_ax.set_title('n cells\n(preview)', fontsize=10, fontweight='bold')
ncell_ax.set_xticks([])
ncell_ax.set_yticks([ypos[t] for t in types])
for spine in ['top','right','left','bottom']:
    ncell_ax.spines[spine].set_visible(False)

# legend order: grouped by gene, ordered alphabetically/numerically within gene, "other" last in each group.
# For KChIP and DPP (too many entries for one column) the legend is split into two side-by-side
# stacked groups: KChIP1+2 | KChIP3+4, and DPP6 | DPP10.
BLANK = '__blank__'  # sentinel: renders as an invisible spacer row in the legend

LEGEND_GROUPS = {
 'Kv4': [['Kv4.1','Kv4.2','Kv4.3S','Kv4.3L']],
 'KChIP': [
   ['KChIP1a','KChIP1b','KChIP1c','KChIP1e','Kcnip1-new_short','Kcnip1 other/fragment',
    BLANK,
    'KChIP2a','KChIP2b','KChIP2c','KChIP2x','Kcnip2 other/fragment'],
   ['KChIP3a','KChIP3x','Kcnip3-new-cys_rich','Kcnip3 other/fragment',
    BLANK,
    'KChIP4a','KChIP4b','KChIP4bL','KChIP4d','KChIP4e','Kcnip4 other/noncoding'],
 ],
 'DPP': [
   ['DPP6a','DPP6K','DPP6L','DPP6S','Dpp6 other/fragment'],
   ['DPP10a','DPP10c','DPP10d','Dpp10 other/fragment'],
 ],
}

def legend_patches(vlist):
    pats = []
    for v in vlist:
        if v == BLANK:
            pats.append(mpatches.Patch(facecolor='none', edgecolor='none', label=' '))
            continue
        if v not in VSPEC:
            continue
        base, kind, pcolor = VSPEC[v]
        hatch = {'bars':'---','dots':'..','X':'xx','O':'OO','solid':None}[kind]
        pats.append(mpatches.Patch(facecolor=base, edgecolor=pcolor, hatch=hatch, label=v, linewidth=0.3))
    return pats

for ax, panel in zip(bar_axes, PANELS):
    groups = LEGEND_GROUPS[panel]
    n = len(groups)
    xs = [0.5] if n == 1 else [0.28, 0.76]
    legends = []
    for xpos_, vlist in zip(xs, groups):
        pats = legend_patches(vlist)
        leg = ax.legend(handles=pats, loc='upper center', bbox_to_anchor=(xpos_, -0.01),
                         ncol=1, fontsize=7.5, frameon=False, handlelength=3.0, handleheight=2.4)
        legends.append(leg)
    for leg in legends[:-1]:
        ax.add_artist(leg)

fig.suptitle('Pilot v6: per-type variant composition — Kv4 / KChIP / DPP, + n-cells-quantified column\n'
             'repeating X hatch / repeating O hatch (yellow, "other/fragment"); retained-intron '
             'transcripts dropped from denominator; 460-cell preview, pooled-read fractions, most-common segment first',
             fontsize=10)
fig.subplots_adjust(left=0.20, right=0.98, top=0.90, bottom=0.40, wspace=0.08)
fig.savefig('/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/pilot_variant_scheme_v17.png', dpi=200)
print('done', len(types))
