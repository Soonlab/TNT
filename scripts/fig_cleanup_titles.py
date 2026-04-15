"""
Cleanup pass:
  - Fix Fig 2B: push good/bad x-tick labels lower (pad)
  - Fix Fig 2C: move legend further right (bbox_to_anchor outside axis)
  - Remove titles from Fig 1 and Fig 5 panels
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch
from scipy.stats import mannwhitneyu

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'; GOLD_DEEP = '#d4a300'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}
PAL_STAGE_DEEP = {'T2':'#7fb0c4', 'T2/T3':'#1c5d7e', 'T3':'#0e2a47', 'T4':'#a01b2b'}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'

# ============================================================
# Load data for Fig 2 B, C
# ============================================================
clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(ROOT/'02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
integ = pd.read_csv(ROOT/'tables/integrated_subject_master.tsv', sep='\t')
rfeat = pd.read_csv(ROOT/'tables/response_feature_stats.tsv', sep='\t')
rf_imp = pd.read_csv(ROOT/'10_ml_predictor/rf_feature_importance.tsv', sep='\t')
UNMATCHED = [13,15,16,17,18,19,33]

def raincloud_deep(ax, data, x, y, order, palette, width=0.4, alpha=0.75):
    for i, g in enumerate(order):
        vals = pd.to_numeric(data[data[x]==g][y], errors='coerce').dropna().values
        if len(vals) < 2: continue
        color = palette[g]
        parts = ax.violinplot(vals, positions=[i-0.15], widths=width, showmeans=False,
                              showextrema=False, showmedians=False)
        for pc in parts['bodies']:
            pc.set_facecolor(color); pc.set_alpha(0.55); pc.set_edgecolor(color); pc.set_linewidth(0.9)
            m = np.mean(pc.get_paths()[0].vertices[:,0])
            pc.get_paths()[0].vertices[:,0] = np.clip(pc.get_paths()[0].vertices[:,0], -np.inf, m)
        bp = ax.boxplot([vals], positions=[i+0.08], widths=0.15, patch_artist=True, showfliers=False,
                        medianprops=dict(color='white', linewidth=1.6),
                        boxprops=dict(facecolor=color, alpha=0.95, edgecolor=color, linewidth=0.9),
                        whiskerprops=dict(color=color, linewidth=0.9),
                        capprops=dict(color=color, linewidth=0.9))
        jitter_x = i + 0.28 + np.random.uniform(-0.07, 0.07, size=len(vals))
        ax.scatter(jitter_x, vals, s=18, color=color, alpha=alpha, edgecolor='white', linewidth=0.5, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=11, color='#0e2a47')
    ax.tick_params(axis='x', pad=10)  # push labels DOWN

# ============================================================
# Fig 2B — TMB raincloud with lower x-labels
# ============================================================
fig, ax = plt.subplots(figsize=(4.8, 4.5))
pre_m = tmb[(tmb.timepoint=='pre') & (~tmb.subject_id.isin(UNMATCHED))]
raincloud_deep(ax, pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ['good','bad'], PAL)
ax.axhline(10, color='#c01a1a', ls='--', lw=1, alpha=0.85)
ax.text(1.5, 10, 'TMB-high (10/Mb)', fontsize=8.5, color='#c01a1a', va='bottom', ha='right', fontweight='bold')
g = pre_m[pre_m.response_bin=='good'].TMB_nonsyn_per_Mb
b = pre_m[pre_m.response_bin=='bad'].TMB_nonsyn_per_Mb
p = mannwhitneyu(g, b).pvalue
stat_bracket(ax, 0, 1, max(g.max(), b.max())+0.3, p)
ax.set_ylabel('Nonsynonymous TMB (/Mb)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlabel('')
# n annotation with more vertical space
ymin = ax.get_ylim()[0]
for i, gname in enumerate(['good','bad']):
    n = pre_m[pre_m['response_bin']==gname]['TMB_nonsyn_per_Mb'].dropna().shape[0]
    ax.text(i, ymin - 0.45, f'n={n}', ha='center', va='top', fontsize=8.5, color='#5a6772')
save_panel(fig, 'Fig2B_TMB_raincloud', OUT)

# ============================================================
# Fig 2C — MSI × TMB scatter, legend moved right outside axis
# ============================================================
fig, ax = plt.subplots(figsize=(6, 4.8))
mdf = msi.merge(tmb[['sample_id','TMB_nonsyn_per_Mb']], on='sample_id', how='left')
mdf = mdf.dropna(subset=['TMB_nonsyn_per_Mb','MSI_pct'])
for resp in ['good','bad']:
    sub = mdf[mdf.response_bin==resp]
    ax.scatter(sub.TMB_nonsyn_per_Mb, sub.MSI_pct, color=PAL[resp], s=85, alpha=0.85,
               edgecolor='white', linewidth=1.0, label=f'{resp} (n={len(sub)})', zorder=3)
ax.axhspan(20, 100, alpha=0.07, color=GOOD_DEEP)
ax.axvspan(10, 100, alpha=0.07, color=GOLD_DEEP)
ax.axhline(20, color=GOOD_DEEP, ls='--', lw=0.9, alpha=0.6)
ax.axvline(10, color=GOLD_DEEP, ls='--', lw=0.9, alpha=0.6)
ax.text(9, 0.01, 'MSS zone\n(TMB-low)', fontsize=9.5, color='#0e2a47', ha='right', va='bottom',
        style='italic', fontweight='bold')
ax.set_xlabel('Nonsynonymous TMB (/Mb)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('MSI percentage (%)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(-0.3, max(tmb.TMB_nonsyn_per_Mb.max()+1, 13))
ax.set_ylim(-0.02, 0.8)
# Legend OUTSIDE axis to the right
ax.legend(fontsize=10, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False)
save_panel(fig, 'Fig2C_MSI_TMB_scatter', OUT)

# ============================================================
# Fig 1 re-regenerate without titles
# ============================================================
# 1A — Sankey
fig, ax = plt.subplots(figsize=(7, 4.5))
from matplotlib.patches import Rectangle as Rect, Polygon
col_x = {0: 0.05, 1: 0.35, 2: 0.65, 3: 0.95}

def box(ax, x, y_bot, y_top, color, label, count, textcolor='white'):
    w = 0.08
    ax.add_patch(Rect((x-w/2, y_bot), w, y_top-y_bot,
                       facecolor=color, edgecolor='#0e2a47', linewidth=1.2))
    ax.text(x, (y_bot+y_top)/2, f'{label}\n(n={count})',
            ha='center', va='center', fontsize=9, color=textcolor, fontweight='bold')

def sankey_flow(ax, x1, y1_bot, y1_top, x2, y2_bot, y2_top, color, alpha=0.4):
    w = 0.08
    x1r = x1 + w/2; x2l = x2 - w/2
    path_vals = np.array([
        [x1r, y1_bot], [x1r, y1_top],
        [(x1r+x2l)/2, y1_top], [(x1r+x2l)/2, y2_top],
        [x2l, y2_top], [x2l, y2_bot],
        [(x1r+x2l)/2, y2_bot], [(x1r+x2l)/2, y1_bot]
    ])
    cp = [MPath.MOVETO, MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
          MPath.LINETO, MPath.CURVE4, MPath.CURVE4]
    pp = PathPatch(MPath(path_vals, cp), facecolor=color, alpha=alpha, edgecolor='none')
    ax.add_patch(pp)

n_all = 35
box(ax, col_x[0], 0.15, 0.85, BLACK_DEEP, 'LARC\nenrolled', n_all)
n_good, n_bad = 18, 17
top = 0.85; good_frac = n_good/n_all
good_top = top; good_bot = top - good_frac*0.7
bad_top = good_bot; bad_bot = 0.15
box(ax, col_x[1], good_bot, good_top, GOOD_DEEP, 'Good\n(TRG 0-1)', n_good)
box(ax, col_x[1], bad_bot, bad_top, BAD_DEEP, 'Poor\n(TRG 2-3)', n_bad)
sankey_flow(ax, col_x[0], good_bot, good_top, col_x[1], good_bot, good_top, GOOD_DEEP, 0.45)
sankey_flow(ax, col_x[0], 0.15, good_bot, col_x[1], bad_bot, bad_top, BAD_DEEP, 0.45)
paired_good = 7; paired_bad = 7
unp_good = 11; unp_bad = 10
paired_top = 0.85; paired_bot = 0.85 - (paired_good+paired_bad)/n_all*0.7
unp_top = paired_bot; unp_bot = 0.15
box(ax, col_x[2], paired_bot, paired_top, BLACK_DEEP, 'Paired\npre+post', paired_good+paired_bad)
box(ax, col_x[2], unp_bot, unp_top, '#5a6772', 'Unpaired\npre only', unp_good+unp_bad)
g_paired_frac = paired_good/n_good
g_mid = good_bot + (good_top-good_bot)*(1-g_paired_frac)
sankey_flow(ax, col_x[1], g_mid, good_top, col_x[2], paired_bot+(paired_top-paired_bot)*(paired_bad/(paired_good+paired_bad)), paired_top, GOOD_DEEP, 0.35)
sankey_flow(ax, col_x[1], good_bot, g_mid, col_x[2], unp_bot+(unp_top-unp_bot)*(unp_bad/(unp_good+unp_bad)), unp_top, GOOD_DEEP, 0.3)
b_paired_frac = paired_bad/n_bad
b_mid = bad_bot + (bad_top-bad_bot)*(1-b_paired_frac)
sankey_flow(ax, col_x[1], b_mid, bad_top, col_x[2], paired_bot, paired_bot+(paired_top-paired_bot)*(paired_bad/(paired_good+paired_bad)), BAD_DEEP, 0.35)
sankey_flow(ax, col_x[1], bad_bot, b_mid, col_x[2], unp_bot, unp_bot+(unp_top-unp_bot)*(unp_bad/(unp_good+unp_bad)), BAD_DEEP, 0.3)
box(ax, col_x[3], 0.55, 0.85, '#00567d', 'WES\n77 samples', 77)
box(ax, col_x[3], 0.15, 0.45, '#c96806', 'RNA-seq\n56 samples', 56)
sankey_flow(ax, col_x[2], 0.15, 0.85, col_x[3], 0.55, 0.85, '#00567d', 0.3)
sankey_flow(ax, col_x[2], 0.15, 0.85, col_x[3], 0.15, 0.45, '#c96806', 0.3)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks([]); ax.set_yticks([])
for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
save_panel(fig, 'Fig1A_sankey', OUT)

# 1B — Waterfall (no titles)
fig = plt.figure(figsize=(13, 5))
gs = fig.add_gridspec(3, 1, height_ratios=[4, 0.6, 1.3], hspace=0.18)
ax_main = fig.add_subplot(gs[0]); ax_anno = fig.add_subplot(gs[1]); ax_legend = fig.add_subplot(gs[2])
order = clin.sort_values(['response_num','subject_id']).reset_index(drop=True)
n = len(order); x = np.arange(n)
trg_colors = {0: '#0f8b78', 1: '#4fb3a3', 2: '#f4a261', 3: '#c1272d'}
bar_colors = [trg_colors[r] for r in order['response_num']]
bars = ax_main.bar(x, order['response_num'], color=bar_colors, edgecolor='white', linewidth=1.5, width=0.85)
for i, (_, r) in enumerate(order.iterrows()):
    ax_main.text(i, r['response_num']+0.05, f"TRG{r['response_num']}", ha='center', va='bottom', fontsize=7, color='#0e2a47')
ax_main.axhline(1.5, color='#0e2a47', ls='--', lw=0.9, alpha=0.6)
ax_main.text(n-0.5, 1.55, 'good / poor boundary', fontsize=8.5, color='#0e2a47', ha='right', va='bottom', fontweight='bold')
ax_main.set_ylim(0, 3.5)
ax_main.set_ylabel('Tumor regression grade (TRG)', fontsize=11, fontweight='bold', color='#0e2a47')
ax_main.set_xlim(-0.6, n-0.4); ax_main.set_xticks([])
add_axis_spines(ax_main)

cT_map = {'T2':0, 'T2/T3':1, 'T3':2, 'T4':3}
ages = order['age'].values; age_norm = (ages - ages.min())/(ages.max()-ages.min())
anno_mat = np.zeros((3, n, 3))
age_cmap = plt.cm.YlOrRd
for j, (_, r) in enumerate(order.iterrows()):
    anno_mat[0, j] = matplotlib.colors.to_rgb(PAL_STAGE_DEEP.get(r['cT'], '#5a6772'))
    anno_mat[1, j] = matplotlib.colors.to_rgb('#0e2a47' if r['sex']=='M' else '#a01b2b')
    anno_mat[2, j] = age_cmap(age_norm[j])[:3]
ax_anno.imshow(anno_mat, aspect='auto', interpolation='nearest')
ax_anno.set_yticks([0,1,2]); ax_anno.set_yticklabels(['cT stage','Sex','Age'], fontsize=9, color='#0e2a47', fontweight='bold')
ax_anno.set_xticks(x); ax_anno.set_xticklabels([f"S{r['subject_id']}" for _, r in order.iterrows()], fontsize=6, rotation=90)
ax_anno.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_anno.spines[s].set_visible(False)

ax_legend.axis('off')
leg_items = []
for t, c in [('TRG 0 (CR)', trg_colors[0]), ('TRG 1 (near-CR)', trg_colors[1]),
             ('TRG 2 (PR)', trg_colors[2]), ('TRG 3 (poor)', trg_colors[3])]:
    leg_items.append(mpatches.Patch(color=c, label=t))
for t, c in PAL_STAGE_DEEP.items():
    leg_items.append(mpatches.Patch(color=c, label=t))
leg_items.append(mpatches.Patch(color='#0e2a47', label='Male'))
leg_items.append(mpatches.Patch(color='#a01b2b', label='Female'))
ax_legend.legend(handles=leg_items, ncol=5, loc='center', fontsize=9, frameon=False,
                 bbox_to_anchor=(0.5, 0.5))
save_panel(fig, 'Fig1B_waterfall', OUT)

# 1C — Clinical comparison
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8),
                         gridspec_kw=dict(width_ratios=[1.4, 1, 1], wspace=0.45))
ax = axes[0]
raincloud_deep(ax, clin, 'response_bin', 'age', ['good','bad'], PAL)
p_age = mannwhitneyu(clin[clin.response_bin=='good'].age, clin[clin.response_bin=='bad'].age).pvalue
stat_bracket(ax, 0, 1, clin.age.max()+2, p_age)
ax.set_ylabel('Age (years)', fontsize=11, fontweight='bold', color='#0e2a47')

ax = axes[1]
ct_tab = pd.crosstab(clin.cT, clin.response_bin, normalize='columns')*100
ct_tab = ct_tab.reindex(['T2','T2/T3','T3','T4']).fillna(0)[['good','bad']]
bottom = np.zeros(2)
for i, t in enumerate(['T2','T2/T3','T3','T4']):
    ax.bar(range(2), ct_tab.loc[t].values, bottom=bottom, label=t,
           color=PAL_STAGE_DEEP[t], edgecolor='white', linewidth=1.5, width=0.6)
    for j in range(2):
        if ct_tab.loc[t].values[j] > 8:
            ax.text(j, bottom[j]+ct_tab.loc[t].values[j]/2, f"{ct_tab.loc[t].values[j]:.0f}%",
                    ha='center', va='center', fontsize=8.5, color='white', fontweight='bold')
    bottom += ct_tab.loc[t].values
ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'])
ax.tick_params(axis='x', pad=10)
ax.set_ylabel('Percentage (%)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylim(0,102)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=9, title='cT stage',
          title_fontsize=9.5, frameon=False)

ax = axes[2]
sex_tab = pd.crosstab(clin.sex, clin.response_bin)
for i, sex in enumerate(['M','F']):
    vals = sex_tab.loc[sex].values if sex in sex_tab.index else [0,0]
    ax.bar([j + (i-0.5)*0.35 for j in range(2)], vals, width=0.32,
           color='#0e2a47' if sex=='M' else '#a01b2b',
           edgecolor='white', linewidth=1.2, label=sex)
    for j, v in enumerate(vals):
        ax.text(j + (i-0.5)*0.35, v+0.3, str(v), ha='center', fontsize=9, fontweight='bold')
ax.set_xticks(range(2)); ax.set_xticklabels(['good','bad'])
ax.tick_params(axis='x', pad=10)
ax.set_ylabel('Patient count', fontsize=11, fontweight='bold', color='#0e2a47')
ax.legend(title='Sex', fontsize=9, loc='upper right', frameon=False)
save_panel(fig, 'Fig1C_clinical', OUT)

# 1D — Sample matrix (retain current layout but no title)
fig = plt.figure(figsize=(14.5, 6.5))
gs = fig.add_gridspec(8, 1, height_ratios=[0.35,0.35,0.35,0.35,0.35,0.35,0.2,5], hspace=0.15)
subs = sorted(clin.subject_id)
order_subj = clin.sort_values(['response_num','subject_id']).subject_id.tolist()

def ann_row(ax, values, label, palette_map=None, cmap=None, vmin=None, vmax=None, is_cat=True, na_color='#ecf0f1'):
    if is_cat:
        colors_arr = np.array([[matplotlib.colors.to_rgb(palette_map.get(v, na_color)) for v in values]])
    else:
        values_num = pd.to_numeric(pd.Series(values), errors='coerce')
        if vmin is None: vmin = np.nanmin(values_num)
        if vmax is None: vmax = np.nanmax(values_num)
        values_num = values_num.fillna(vmin)
        norm = (values_num - vmin)/(vmax-vmin+1e-9)
        colors_arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(colors_arr, aspect='auto', interpolation='nearest')
    ax.set_yticks([0]); ax.set_yticklabels([label], fontsize=9, color='#0e2a47', fontweight='bold')
    ax.set_xticks([]); ax.tick_params(length=0)
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)

ann_data = pd.DataFrame({'subject_id': order_subj})
ann_data = ann_data.merge(clin[['subject_id','response_bin','cT','sex','age']], on='subject_id')
tmb_pre = tmb[tmb.timepoint=='pre'][['subject_id','TMB_nonsyn_per_Mb']].drop_duplicates('subject_id')
ann_data = ann_data.merge(tmb_pre, on='subject_id', how='left')
msi_pre = msi[msi.timepoint=='pre'][['subject_id','MSI_pct']].drop_duplicates('subject_id')
ann_data = ann_data.merge(msi_pre, on='subject_id', how='left')

ax0 = fig.add_subplot(gs[0]); ann_row(ax0, ann_data.response_bin.tolist(), 'Response', PAL)
ax1 = fig.add_subplot(gs[1]); ann_row(ax1, ann_data.cT.tolist(), 'cT stage', PAL_STAGE_DEEP)
ax2 = fig.add_subplot(gs[2]); ann_row(ax2, ann_data.sex.tolist(), 'Sex', {'M':'#0e2a47','F':'#a01b2b'})
ax3 = fig.add_subplot(gs[3]); ann_row(ax3, ann_data.age.tolist(), 'Age', cmap=plt.cm.YlOrRd, is_cat=False, vmin=30, vmax=80)
ax4 = fig.add_subplot(gs[4]); ann_row(ax4, ann_data.TMB_nonsyn_per_Mb.tolist(), 'TMB/Mb', cmap=plt.cm.Purples, is_cat=False, vmin=0, vmax=3)
ax5 = fig.add_subplot(gs[5]); ann_row(ax5, ann_data.MSI_pct.tolist(), 'MSI %', cmap=plt.cm.Greens, is_cat=False, vmin=0, vmax=0.3)

ax_main = fig.add_subplot(gs[7])
mat = np.zeros((6, len(order_subj)))
for j, s in enumerate(order_subj):
    w = wes_inv[wes_inv.subject_id==s]; r = rna_inv[rna_inv.subject_id==s]
    mat[0,j] = 1 if (w.timepoint=='normal').any() else 0
    mat[1,j] = 1 if (w.timepoint=='pre').any() else 0
    mat[2,j] = 1 if (w.timepoint=='post').any() else 0
    mat[3,j] = 1 if (r.timepoint=='normal').any() else 0
    mat[4,j] = 1 if (r.timepoint=='pre').any() else 0
    mat[5,j] = 1 if (r.timepoint=='post').any() else 0
color_mat = np.ones((6, len(order_subj), 3))
sample_colors = ['#00567d','#0e2a47','#7fb0c4','#c96806','#c53e1f','#f4a261']
for i in range(6):
    for j in range(len(order_subj)):
        if mat[i,j]:
            color_mat[i,j] = matplotlib.colors.to_rgb(sample_colors[i])
        else:
            color_mat[i,j] = matplotlib.colors.to_rgb('#ecf0f1')
ax_main.imshow(color_mat, aspect='auto', interpolation='nearest')
ax_main.set_yticks(range(6))
ax_main.set_yticklabels(['WES normal','WES pre','WES post','RNA normal','RNA pre','RNA post'],
                        fontsize=9.5, color='#0e2a47', fontweight='bold')
ax_main.set_xticks(range(len(order_subj)))
ax_main.set_xticklabels([f'S{s}' for s in order_subj], fontsize=6.5, rotation=90, color='#0e2a47')
ax_main.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_main.spines[s].set_visible(False)
save_panel(fig, 'Fig1D_sample_matrix', OUT)

# 1E — Study design (no title)
fig, ax = plt.subplots(figsize=(7, 6))
ax.axis('off')

def flow_box(x, y, w, h, text, color, textcolor='white', bold=True):
    ax.add_patch(FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.02",
                 linewidth=1.3, edgecolor='#0e2a47', facecolor=color))
    ax.text(x, y, text, ha='center', va='center', fontsize=9, color=textcolor,
            fontweight='bold' if bold else 'normal')

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=1.3, color='#0e2a47'))

flow_box(0.5, 0.95, 0.42, 0.08, 'LARC patients enrolled (N = 35)', BLACK_DEEP)
flow_box(0.22, 0.78, 0.30, 0.08, 'Good responders\nTRG 0-1  (n = 18)', GOOD_DEEP)
flow_box(0.78, 0.78, 0.30, 0.08, 'Poor responders\nTRG 2-3  (n = 17)', BAD_DEEP)
arrow(0.45, 0.91, 0.30, 0.82); arrow(0.55, 0.91, 0.70, 0.82)
flow_box(0.22, 0.60, 0.30, 0.09, 'WES 38 samples\nRNA-seq 33 samples', '#00567d')
flow_box(0.78, 0.60, 0.30, 0.09, 'WES 39 samples\nRNA-seq 23 samples', '#00567d')
arrow(0.22, 0.74, 0.22, 0.65); arrow(0.78, 0.74, 0.78, 0.65)
flow_box(0.5, 0.42, 0.55, 0.08, 'Matched pre/post pairs (n = 14)\nUnmatched (n = 21)', '#5a6772', 'white')
arrow(0.22, 0.55, 0.42, 0.46); arrow(0.78, 0.55, 0.58, 0.46)
flow_box(0.5, 0.22, 0.9, 0.16,
         'Integrated analyses:\nWES somatic · SBS signatures · MSI · CNV · HLA class I · HLA LOH · neoantigens\nRNA DEG · GSEA · ssGSEA · 22 immune signatures · TCR/BCR · CMS\nClonal evolution · ML LOOCV predictor · 7-cohort external meta-analysis',
         '#0e2a47', 'white', bold=False)
arrow(0.5, 0.38, 0.5, 0.30)
flow_box(0.5, 0.04, 0.9, 0.06,
         'Primary finding: DNA-repair / cell-cycle programs predict response; treatment cascade in responders only.',
         GOLD_DEEP, '#0e2a47')
arrow(0.5, 0.14, 0.5, 0.07)
ax.set_xlim(0,1); ax.set_ylim(0,1)
save_panel(fig, 'Fig1E_design', OUT)

# ============================================================
# Fig 5 — remove titles from panels B, C, D, E, F
# (regenerate relevant pieces without titles)
# ============================================================
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc
import warnings; warnings.filterwarnings('ignore')
import umap
from sklearn.decomposition import PCA
import shap

def assign_omic(f):
    if f in ['TMB_nonsyn_per_Mb','n_nonsyn','MMR_prop','SBS5','SBS3','MSI_pct']: return 'WES_genomic'
    if f in ['CIN','frac_amp','frac_del']: return 'WES_CNV'
    if 'CD8' in f or 'MHC' in f or 'IFN' in f or 'TLS' in f or 'NLRC' in f or 'NK_' in f or 'B_' in f or 'Mac_' in f or 'Treg' in f: return 'RNA_immune'
    if 'EMT' in f or 'TGF' in f or 'CAF' in f or 'Hypox' in f or 'Stem' in f or 'Epith' in f: return 'RNA_stromal'
    if 'Repair' in f or 'HDR' in f or 'E2F' in f or 'G2-M' in f or 'Myc' in f or 'Hallmark' in f: return 'RNA_pathway'
    if f in ['age']: return 'Clinical'
    return 'Other'
OMIC_COLORS = {
    'Clinical':'#0e2a47','WES_genomic':'#118ab2','WES_CNV':'#06aed5',
    'RNA_pathway':'#057a64','RNA_immune':'#c11456','RNA_stromal':'#b03219','Other':'#5a6772',
}
feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
omic_map = {f: assign_omic(f) for f in feats_num}
X_raw = integ[feats_num].apply(pd.to_numeric, errors='coerce')
X_imputed = X_raw.fillna(X_raw.median())
y = (integ['response_bin']=='good').astype(int).values
X = X_raw.values

# 5B ROC
models = {
    'LASSO LR':  Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='l1', solver='saga', C=0.5, max_iter=10000))]),
    'Elastic Net LR': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, C=0.5, max_iter=10000))]),
    'Random Forest': Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
                           ('clf',RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))]),
}
model_colors = {'LASSO LR':'#118ab2','Elastic Net LR':'#06aed5','Random Forest':'#7a3aad'}
loo = LeaveOneOut()

fig, ax = plt.subplots(figsize=(6.5, 6))
np.random.seed(42); mean_fpr = np.linspace(0, 1, 100)
for name, mdl in models.items():
    proba = cross_val_predict(mdl, X, y, cv=loo, method='predict_proba')[:, 1]
    fpr, tpr, _ = roc_curve(y, proba); auc_val = auc(fpr, tpr)
    bs_tprs = []; bs_aucs = []
    for _ in range(200):
        idx = np.random.choice(len(y), len(y), replace=True)
        if len(set(y[idx])) < 2: continue
        try:
            f_b, t_b, _ = roc_curve(y[idx], proba[idx])
            interp_t = np.interp(mean_fpr, f_b, t_b); interp_t[0] = 0
            bs_tprs.append(interp_t); bs_aucs.append(auc(f_b, t_b))
        except: pass
    bs_tprs = np.array(bs_tprs)
    if len(bs_tprs) > 0:
        tpr_lo = np.percentile(bs_tprs, 2.5, axis=0)
        tpr_hi = np.percentile(bs_tprs, 97.5, axis=0)
        ax.fill_between(mean_fpr, tpr_lo, tpr_hi, color=model_colors[name], alpha=0.18)
    auc_lo, auc_hi = np.percentile(bs_aucs, [2.5, 97.5])
    ax.plot(fpr, tpr, color=model_colors[name], lw=2.6,
            label=f'{name}\n  AUC = {auc_val:.3f}  [95% CI: {auc_lo:.2f}–{auc_hi:.2f}]')
ax.plot([0,1],[0,1], ls='--', color='#5a6772', lw=1.0, label='Random (AUC = 0.5)')
ax.set_xlabel('False positive rate (1 − specificity)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('True positive rate (sensitivity)', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.legend(loc='lower right', fontsize=9, frameon=True, framealpha=0.95, edgecolor='#0e2a47')
add_axis_spines(ax)
save_panel(fig, 'Fig5B_ROC_CI', OUT)

# 5C Forest
fig, ax = plt.subplots(figsize=(9, 8))
top = rfeat.head(20).iloc[::-1].reset_index(drop=True)
np.random.seed(42)
def bootstrap_delta(g_vals, b_vals, n=1000):
    d = []
    for _ in range(n):
        gi = np.random.choice(g_vals, len(g_vals), replace=True)
        bi = np.random.choice(b_vals, len(b_vals), replace=True)
        d.append(np.median(gi) - np.median(bi))
    return np.percentile(d, [2.5, 97.5])
cis = []
for _, r in top.iterrows():
    if r.feature not in integ.columns: cis.append((0,0)); continue
    vals = pd.to_numeric(integ[r.feature], errors='coerce')
    g_vals = vals[integ.response_bin=='good'].dropna().values
    b_vals = vals[integ.response_bin=='bad'].dropna().values
    if len(g_vals)>=3 and len(b_vals)>=3:
        cis.append(bootstrap_delta(g_vals, b_vals))
    else: cis.append((0,0))
top['ci_low'] = [c[0] for c in cis]; top['ci_high'] = [c[1] for c in cis]
y_pos = np.arange(len(top))
xmax_data = max(top.ci_high.max(), top.delta_med.max())
xmin_data = min(top.ci_low.min(), top.delta_med.min())
right_margin = (xmax_data - xmin_data) * 0.55
for i, (_, r) in enumerate(top.iterrows()):
    color = GOOD_DEEP if r.delta_med>0 else BAD_DEEP
    ax.plot([r.ci_low, r.ci_high], [i, i], color=color, lw=2.2, alpha=0.85, solid_capstyle='round')
    ax.scatter(r.delta_med, i, s=180 if r.pvalue<0.05 else 100, color=color, edgecolor='white', linewidth=1.4, zorder=3)
    star = sig_symbol(r.pvalue); star = star if star != 'ns' else ''
    label = f'p = {r.pvalue:.3g} {star}'
    ax.text(xmax_data + right_margin*0.04, i, label, va='center', ha='left', fontsize=9, color='#0e2a47')
ax.axvline(0, color='#0e2a47', lw=1.0)
ax.axvspan(-0.03, 0.03, color='#dee2e6', alpha=0.4)
ax.set_yticks(y_pos)
ax.set_yticklabels([f[:50] for f in top.feature], fontsize=9.5, color='#0e2a47')
ax.set_xlabel('Median Δ (good − poor),  95% bootstrap CI', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_xlim(xmin_data*1.1, xmax_data + right_margin)
add_axis_spines(ax)
ax.text(0.98, 0.98, '↑ Good      ↑ Poor →', transform=ax.transAxes, ha='right', va='top',
        fontsize=9, color='#0e2a47', fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.4'))
save_panel(fig, 'Fig5C_forest_CI', OUT)

# 5D UMAP
X_scaled = StandardScaler().fit_transform(X_imputed.values)
try:
    reducer = umap.UMAP(n_neighbors=8, min_dist=0.3, random_state=42)
    embedding = reducer.fit_transform(X_scaled)
    method_label = 'UMAP'
except Exception:
    pca = PCA(n_components=2); embedding = pca.fit_transform(X_scaled)
    method_label = 'PCA'

fig, ax = plt.subplots(figsize=(7.5, 6.5))
for resp in ['good','bad']:
    mask = (integ.response_bin == resp).values
    try:
        sns.kdeplot(x=embedding[mask, 0], y=embedding[mask, 1], ax=ax,
                    color=PAL[resp], levels=4, alpha=0.4, linewidths=1.2)
    except: pass
for resp in ['good','bad']:
    mask = (integ.response_bin == resp).values
    ax.scatter(embedding[mask, 0], embedding[mask, 1], color=PAL[resp], s=140, alpha=0.92,
               edgecolor='white', linewidth=1.4, label=f'{resp} (n={mask.sum()})', zorder=4)
for i in range(len(integ)):
    sid = integ.iloc[i].subject_id
    ax.annotate(f'S{int(sid)}', (embedding[i, 0], embedding[i, 1]), xytext=(5, 5),
                textcoords='offset points', fontsize=7, color='#0e2a47', alpha=0.6)
ax.set_xlabel(f'{method_label} dim 1', fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel(f'{method_label} dim 2', fontsize=11, fontweight='bold', color='#0e2a47')
ax.legend(loc='upper right', fontsize=10, frameon=True, framealpha=0.92, edgecolor='#0e2a47')
add_axis_spines(ax)
save_panel(fig, 'Fig5D_UMAP', OUT)

# 5E SHAP beeswarm
rf = Pipeline([('imp',SimpleImputer(strategy='median')),('sc',StandardScaler()),
               ('clf',RandomForestClassifier(n_estimators=500, random_state=0, n_jobs=4))])
rf.fit(X, y)
X_proc = rf.named_steps['sc'].transform(rf.named_steps['imp'].transform(X))
explainer = shap.TreeExplainer(rf.named_steps['clf'])
shap_values = explainer.shap_values(X_proc)
if isinstance(shap_values, list): sv = shap_values[1]
elif shap_values.ndim == 3: sv = shap_values[:, :, 1]
else: sv = shap_values
mean_abs = np.abs(sv).mean(axis=0)
top_idx = np.argsort(mean_abs)[-15:]
top_feat_names = [feats_num[i] for i in top_idx]

fig, ax = plt.subplots(figsize=(9, 7))
np.random.seed(0)
cmap_shap = LinearSegmentedColormap.from_list('shap', ['#118ab2','#dde2e8','#c11456'])
for plot_i, fi in enumerate(top_idx):
    sv_f = sv[:, fi]; val_f = X_proc[:, fi]
    val_norm = (val_f - val_f.min()) / (val_f.max() - val_f.min() + 1e-9)
    colors = cmap_shap(val_norm)
    jitter_y = plot_i + np.random.uniform(-0.3, 0.3, len(sv_f))
    ax.scatter(sv_f, jitter_y, s=45, c=colors, alpha=0.85, edgecolor='#0e2a47', linewidth=0.4, zorder=3)
ax.axvline(0, color='#0e2a47', lw=0.9)
ax.set_yticks(range(len(top_idx)))
ax.set_yticklabels([f[:48] for f in top_feat_names], fontsize=9.5, color='#0e2a47')
ax.set_xlabel('SHAP value (impact on predicted "good responder" probability)',
              fontsize=10.5, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)
import matplotlib.colors as mcolors
norm = mcolors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap_shap, norm=norm); sm.set_array([])
cax_sb = fig.add_axes([0.97, 0.3, 0.012, 0.4])
cb = fig.colorbar(sm, cax=cax_sb, ticks=[0, 1])
cb.ax.set_yticklabels(['low', 'high'], fontsize=8.5, color='#0e2a47')
cb.set_label('Feature value\n(scaled)', fontsize=9.5, color='#0e2a47', fontweight='bold')
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)
save_panel(fig, 'Fig5E_SHAP_beeswarm', OUT)

# 5F per-subject prediction waterfall
mdl_best = models['LASSO LR']
proba_best = cross_val_predict(mdl_best, X, y, cv=loo, method='predict_proba')[:, 1]
pred_df = pd.DataFrame({'subject_id': integ.subject_id, 'true':y, 'true_label': integ.response_bin,
                         'prob_good': proba_best})
pred_df['correct'] = (pred_df.prob_good >= 0.5) == pred_df.true.astype(bool)
pred_df = pred_df.sort_values('prob_good', ascending=False).reset_index(drop=True)
pred_df['height'] = pred_df.prob_good - 0.5

fig, ax = plt.subplots(figsize=(13, 5))
x_pos = np.arange(len(pred_df))
colors = [PAL[lbl] for lbl in pred_df.true_label]
bars = ax.bar(x_pos, pred_df.height, color=colors, edgecolor='white', linewidth=0.8, width=0.85)
ax.axhline(0, color='#0e2a47', lw=1.0); ax.set_ylim(-0.55, 0.55)
ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
ax.set_yticklabels(['0\n(predicted poor)', '0.25', '0.5\n(decision)', '0.75', '1.0\n(predicted good)'],
                   fontsize=8.5, color='#0e2a47')
ax.set_xticks(x_pos)
ax.set_xticklabels([f'S{int(s)}' for s in pred_df.subject_id], fontsize=7, rotation=90, color='#0e2a47')
ax.set_xlim(-0.6, len(pred_df)-0.4)
mis = pred_df[~pred_df.correct]
for _, r in mis.iterrows():
    i = pred_df.index[pred_df.subject_id==r.subject_id][0]
    h = r.height; y_marker = h + (0.04 if h>=0 else -0.04)
    ax.scatter(i, y_marker, marker='x', s=70, color='#c01a1a', linewidth=2.2, zorder=5)
handles = [
    mpatches.Patch(color=GOOD_DEEP, label=f'True good (n={int(pred_df.true.sum())})'),
    mpatches.Patch(color=BAD_DEEP, label=f'True poor (n={int((1-pred_df.true).sum())})'),
    matplotlib.lines.Line2D([0],[0], marker='x', color='#c01a1a', markersize=10, lw=0,
                             label=f'Misclassified (n={int((~pred_df.correct).sum())})'),
]
ax.legend(handles=handles, loc='upper right', fontsize=9.5, frameon=True, framealpha=0.95, edgecolor='#0e2a47')
ax.set_xlabel('Subject (ranked by predicted good-responder probability)',
              fontsize=10.5, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)
acc = pred_df.correct.mean()
ax.text(0.01, 0.96, f'LASSO LR LOOCV accuracy = {acc:.1%}', transform=ax.transAxes,
        ha='left', va='top', fontsize=10, color='#0e2a47', fontweight='bold',
        bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.4'))
save_panel(fig, 'Fig5F_per_subject_prediction', OUT)

print('\n=== Cleanup pass saved (Fig 1 + Fig 2B/C + Fig 5 B-F, titles removed) ===')
