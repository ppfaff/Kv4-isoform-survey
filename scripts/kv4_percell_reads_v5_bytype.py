import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

d = pd.read_pickle('/home/claude/kv4/iso/d.pkl')
counts = pd.read_pickle('/home/claude/kv4/counts.pkl')
anno = pd.read_csv('/home/claude/kv4/anno.csv')

final_map = {
 'Kcnd1-201': 'Kv4.1', 'Kcnd2-201': 'Kv4.2',
 'Kcnd3-202': 'Kv4.3S', 'Kcnd3-201': 'Kv4.3L', 'Kcnd3-203': 'Kv4.3L',
}
VARIANTS = ['Kv4.1', 'Kv4.2', 'Kv4.3S', 'Kv4.3L']
COLOR = {'Kv4.1': '#4FA8D8', 'Kv4.2': '#E8822D', 'Kv4.3S': '#B9A3D9', 'Kv4.3L': '#7B4FA3'}
OFFSET = {'Kv4.1': -0.27, 'Kv4.2': -0.09, 'Kv4.3S': 0.09, 'Kv4.3L': 0.27}

dd = d[d['gene'].isin(['Kcnd1', 'Kcnd2', 'Kcnd3'])].copy()
dd['vclass'] = dd['transcript'].map(final_map)
dd = dd.dropna(subset=['vclass'])

percell = dd.groupby(['sample_name', 'vclass'])['reads'].sum().reset_index()
wide = percell.pivot(index='sample_name', columns='vclass', values='reads').reindex(columns=VARIANTS).fillna(0.0)

meta = dd[['sample_name', 'type', 'primary_type_id', 'cls']].drop_duplicates().set_index('sample_name')
wide = wide.join(meta, how='left')

depth = counts.sum(axis=0)
depth.index.name = 'sample_name'
wide['depth'] = depth.reindex(wide.index)
wide = wide.dropna(subset=['depth'])

for v in VARIANTS:
    wide[v + '_cpm'] = wide[v] / wide['depth'] * 1e6
VCPM = [v + '_cpm' for v in VARIANTS]
wide['total_cpm'] = wide[VCPM].sum(axis=1)

all_types = sorted(wide['type'].unique())
palette = [cm.tab20(i % 20) if i < 20 else cm.tab20b((i - 20) % 20) for i in range(len(all_types))]
TYPE_COLOR = dict(zip(all_types, palette))

CLS_ROWS = [('GABA', 'Inhibitory (GABA)'), ('Glut', 'Excitatory (Glut)'), ('NonN', 'Non-neuronal')]

fig = plt.figure(figsize=(20, 17))
outer = gridspec.GridSpec(3, 1, hspace=0.6, top=0.90, bottom=0.07, left=0.06, right=0.98, figure=fig)
axes = []
for block in range(3):
    inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[block],
                                              height_ratios=[4, 0.35], hspace=0.03)
    axes.append(fig.add_subplot(inner[0]))
    axes.append(fig.add_subplot(inner[1]))

for i, (cls_key, cls_label) in enumerate(CLS_ROWS):
    ax, strip_ax = axes[2 * i], axes[2 * i + 1]
    sub = wide[wide['cls'] == cls_key].copy()
    if sub.empty:
        ax.set_visible(False); strip_ax.set_visible(False)
        continue

    # cluster by type, ordered by that type's mean Kv4 signal (high to low);
    # within each type block, order cells by their own total (high to low, left to right)
    type_mean = sub.groupby('type')['total_cpm'].mean().sort_values(ascending=False)
    type_order = type_mean.index.tolist()

    x = 0.0
    block_bounds = []  # (type, start_x, end_x, mean)
    pieces = []
    for t in type_order:
        block = sub[sub['type'] == t].sort_values('total_cpm', ascending=False).copy()
        block['x'] = np.arange(len(block)) + x
        pieces.append(block)
        block_bounds.append((t, x, x + len(block) - 1, type_mean[t]))
        x += len(block)
    sub = pd.concat(pieces)
    n_total = len(sub)

    for v in VARIANTS:
        vals = sub[v + '_cpm']
        mask = vals > 0
        ax.scatter(sub.loc[mask, 'x'] + OFFSET[v], vals[mask], s=50, color=COLOR[v],
                   edgecolor='black', linewidth=0.3, alpha=0.75, label=v, zorder=3)

    # explicit "no Kv4 detected" marker, instead of leaving a blank gap
    miss = sub[sub['total_cpm'] == 0]
    if len(miss):
        ax.scatter(miss['x'], np.full(len(miss), 0.42), marker='x', color='#888888', s=28,
                   linewidth=1.0, zorder=3, label='no Kv4 detected')

    # per-type mean, drawn as a horizontal segment spanning its block (the "staircase" of block averages)
    for t, x0, x1, m in block_bounds:
        ax.plot([x0 - 0.5, x1 + 0.5], [m, m], color='black', linewidth=1.6, alpha=0.6, zorder=2)
        ax.axvline(x0 - 0.5, color='#CCCCCC', linewidth=0.5, zorder=0.5)

    ax.set_yscale('log')
    ax.set_ylim(0.3, None)
    ax.set_xlim(-1, n_total)
    ax.set_xticks([])
    ax.set_title(f'{cls_label}  (n = {n_total} cells, {len(type_order)} types, '
                 f'blocks ordered by type-mean Kv4 signal)', fontsize=11, fontweight='bold', loc='left')
    ax.set_ylabel('Kv4 reads per\nmillion total reads', fontsize=9)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    # type-identity strip + type name label at block center
    strip_colors = [TYPE_COLOR[t] for t in sub['type']]
    strip_ax.bar(sub['x'], height=1.0, width=1.0, color=strip_colors, edgecolor='none')
    strip_ax.set_xlim(-1, n_total)
    strip_ax.set_ylim(0, 1)
    strip_ax.set_xticks([]); strip_ax.set_yticks([])
    for spine in strip_ax.spines.values():
        spine.set_visible(False)
    for t, x0, x1, m in block_bounds:
        cx = (x0 + x1) / 2
        strip_ax.text(cx, -1.6, t, rotation=90, fontsize=9, ha='center', va='top',
                       color=TYPE_COLOR[t], fontweight='bold')

axes[-2].set_xlabel('individual cells, grouped into cell-type blocks (block = type-mean Kv4, high to low); '
                     'within block, high to low Kv4 left to right', fontsize=10)
handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, loc='upper center', fontsize=13, framealpha=0.9, ncol=5)

fig.suptitle('Kv4 pore-subunit signal per cell, clustered by cell type\n'
             'blocks ordered by type-mean Kv4 (high → low); black step = per-type mean; '
             'gray X = cell with no detectable Kv4\n460-cell preview, Tasic 2016 (GSE71585)', fontsize=11)
fig.savefig('/tmp/claude-0/-home-claude/85b0e16d-66ee-57f3-a0b0-9b216861ee30/scratchpad/kv4_percell_reads_v5_bytype.png', dpi=200)
print('done')
