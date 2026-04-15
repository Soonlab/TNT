"""
Figure 1 — Cohort & study design
  1A: Sankey-style cohort flow (matplotlib custom drawing)
  1B: Waterfall / swimmer plot (35 subjects × TRG + clinical annotation tiers)
  1C: Clinical comparison (raincloud age + stacked cT + sex dot)
  1D: Sample × feature heatmap with rich multi-tier column annotation
  1E: CONSORT-style study design
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v2'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(ROOT/'02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')

# ==========================================================
# 1A — Sankey flow
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 4.5))
from matplotlib.patches import Rectangle, Polygon

# Columns: x positions
col_x = {0: 0.05, 1: 0.35, 2: 0.65, 3: 0.95}

def box(ax, x, y_bot, y_top, color, label, count, textcolor='white'):
    w = 0.08
    ax.add_patch(Rectangle((x-w/2, y_bot), w, y_top-y_bot,
                            facecolor=color, edgecolor='#1d3557', linewidth=1.2))
    ax.text(x, (y_bot+y_top)/2, f'{label}\n(n={count})',
            ha='center', va='center', fontsize=9, color=textcolor, fontweight='bold')

def sankey_flow(ax, x1, y1_bot, y1_top, x2, y2_bot, y2_top, color, alpha=0.4):
    w = 0.08
    x1r = x1 + w/2; x2l = x2 - w/2
    path = np.array([
        [x1r, y1_bot], [x1r, y1_top],
        [(x1r+x2l)/2, y1_top], [(x1r+x2l)/2, y2_top],
        [x2l, y2_top], [x2l, y2_bot],
        [(x1r+x2l)/2, y2_bot], [(x1r+x2l)/2, y1_bot]
    ])
    # Use Bezier-ish smoothing
    from matplotlib.path import Path as MPath
    from matplotlib.patches import PathPatch
    cp = [MPath.MOVETO, MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
          MPath.LINETO, MPath.CURVE4, MPath.CURVE4]
    path_patch = PathPatch(MPath(path, cp), facecolor=color, alpha=alpha, edgecolor='none')
    ax.add_patch(path_patch)

# Column 0: all patients
n_all = 35
box(ax, col_x[0], 0.15, 0.85, BLACK, 'LARC\nenrolled', n_all)

# Column 1: good / bad
n_good, n_bad = 18, 17
top = 0.85; good_frac = n_good/n_all
good_top = top; good_bot = top - good_frac*0.7
bad_top = good_bot; bad_bot = 0.15
box(ax, col_x[1], good_bot, good_top, GOOD, 'Good\n(TRG 0-1)', n_good)
box(ax, col_x[1], bad_bot, bad_top, BAD, 'Poor\n(TRG 2-3)', n_bad)

# Flows 0→1
sankey_flow(ax, col_x[0], good_bot, good_top, col_x[1], good_bot, good_top, GOOD, 0.45)
sankey_flow(ax, col_x[0], 0.15, good_bot, col_x[1], bad_bot, bad_top, BAD, 0.45)

# Column 2: paired / unpaired
paired_good = 7; paired_bad = 7
unp_good = 11; unp_bad = 10
paired_top = 0.85; paired_bot = 0.85 - (paired_good+paired_bad)/n_all*0.7
unp_top = paired_bot; unp_bot = 0.15
box(ax, col_x[2], paired_bot, paired_top, BLACK, 'Paired\npre+post', paired_good+paired_bad)
box(ax, col_x[2], unp_bot, unp_top, NORMAL, 'Unpaired\npre only', unp_good+unp_bad)

# flows 1→2
# good → paired/unpaired proportions
g_paired_frac = paired_good/n_good
g_mid = good_bot + (good_top-good_bot)*(1-g_paired_frac)
sankey_flow(ax, col_x[1], g_mid, good_top, col_x[2], paired_bot+(paired_top-paired_bot)*(paired_bad/(paired_good+paired_bad)), paired_top, GOOD, 0.35)
sankey_flow(ax, col_x[1], good_bot, g_mid, col_x[2], unp_bot+(unp_top-unp_bot)*(unp_bad/(unp_good+unp_bad)), unp_top, GOOD, 0.3)

b_paired_frac = paired_bad/n_bad
b_mid = bad_bot + (bad_top-bad_bot)*(1-b_paired_frac)
sankey_flow(ax, col_x[1], b_mid, bad_top, col_x[2], paired_bot, paired_bot+(paired_top-paired_bot)*(paired_bad/(paired_good+paired_bad)), BAD, 0.35)
sankey_flow(ax, col_x[1], bad_bot, b_mid, col_x[2], unp_bot, unp_bot+(unp_top-unp_bot)*(unp_bad/(unp_good+unp_bad)), BAD, 0.3)

# Column 3: sample counts
box(ax, col_x[3], 0.55, 0.85, PRE, 'WES\n77 samples', 77)
box(ax, col_x[3], 0.15, 0.45, POST, 'RNA-seq\n56 samples', 56)
sankey_flow(ax, col_x[2], 0.15, 0.85, col_x[3], 0.55, 0.85, PRE, 0.3)
sankey_flow(ax, col_x[2], 0.15, 0.85, col_x[3], 0.15, 0.45, POST, 0.3)

ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks([]); ax.set_yticks([])
for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
ax.set_title('Patient enrollment, response stratification, and multi-omic sampling', pad=15)
save_panel(fig, 'Fig1A_sankey', OUT)

# ==========================================================
# 1B — Waterfall plot (35 patients sorted by TRG)
# ==========================================================
fig = plt.figure(figsize=(13, 5))
gs = fig.add_gridspec(3, 1, height_ratios=[4, 0.6, 1.3], hspace=0.18)
ax_main = fig.add_subplot(gs[0])
ax_anno = fig.add_subplot(gs[1])
ax_legend = fig.add_subplot(gs[2])

# Sort by response then TRG ascending (best first), so bars go low→high
order = clin.sort_values(['response_num','subject_id']).reset_index(drop=True)
n = len(order)
x = np.arange(n)
# TRG colormap: 0=best green, 3=worst red
trg_colors = {0: '#0f8b78', 1: '#4fb3a3', 2: '#f4a261', 3: '#c1272d'}
bar_colors = [trg_colors[r] for r in order['response_num']]
bars = ax_main.bar(x, order['response_num'], color=bar_colors, edgecolor='white', linewidth=1.5, width=0.85)
# Add TRG labels on top of bars
for i, (_, r) in enumerate(order.iterrows()):
    ax_main.text(i, r['response_num']+0.05, f"TRG{r['response_num']}", ha='center', va='bottom', fontsize=7, color='#1d3557')
ax_main.axhline(1.5, color='#1d3557', ls='--', lw=0.9, alpha=0.6)
ax_main.text(n-0.5, 1.55, 'good / poor boundary', fontsize=8, color='#1d3557', ha='right', va='bottom')
ax_main.set_ylim(0, 3.5)
ax_main.set_ylabel('Tumor regression grade (TRG)')
ax_main.set_xlim(-0.6, n-0.4)
ax_main.set_xticks([])
ax_main.set_title('Individual response distribution (n=35; TRG0=complete, TRG3=poor)')
add_axis_spines(ax_main)

# Annotation strip: cT, sex, age group
cT_map = {'T2':0, 'T2/T3':1, 'T3':2, 'T4':3}
ages = order['age'].values; age_norm = (ages - ages.min())/(ages.max()-ages.min())
anno_rows = 3
anno_mat = np.zeros((anno_rows, n, 3))
age_cmap = plt.cm.YlOrRd
for j, (_, r) in enumerate(order.iterrows()):
    # row 0: cT stage
    anno_mat[0, j] = matplotlib.colors.to_rgb(PAL_STAGE.get(r['cT'], NORMAL))
    # row 1: sex
    anno_mat[1, j] = matplotlib.colors.to_rgb('#264653' if r['sex']=='M' else '#e63946')
    # row 2: age (gradient)
    anno_mat[2, j] = age_cmap(age_norm[j])[:3]
ax_anno.imshow(anno_mat, aspect='auto', interpolation='nearest')
ax_anno.set_yticks([0,1,2])
ax_anno.set_yticklabels(['cT stage','Sex','Age'], fontsize=8)
ax_anno.set_xticks(x)
ax_anno.set_xticklabels([f"S{r['subject_id']}" for _, r in order.iterrows()], fontsize=6, rotation=90)
ax_anno.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_anno.spines[s].set_visible(False)

# Legend
ax_legend.axis('off')
leg_items = []
for t, c in [('TRG 0 (CR)', trg_colors[0]), ('TRG 1 (near-CR)', trg_colors[1]), ('TRG 2 (PR)', trg_colors[2]), ('TRG 3 (poor)', trg_colors[3])]:
    leg_items.append(mpatches.Patch(color=c, label=t))
for t, c in [('cT2', PAL_STAGE['T2']), ('cT2/3', PAL_STAGE['T2/T3']), ('cT3', PAL_STAGE['T3']), ('cT4', PAL_STAGE['T4'])]:
    leg_items.append(mpatches.Patch(color=c, label=t))
for t, c in [('Male', '#264653'), ('Female', '#e63946')]:
    leg_items.append(mpatches.Patch(color=c, label=t))
ax_legend.legend(handles=leg_items, ncol=5, loc='center', fontsize=8.5, frameon=False,
                 bbox_to_anchor=(0.5, 0.5))

save_panel(fig, 'Fig1B_waterfall', OUT)

# ==========================================================
# 1C — Clinical summary (age raincloud + cT stacked + sex dot)
# ==========================================================
fig, axes = plt.subplots(1, 3, figsize=(11, 3.8),
                         gridspec_kw=dict(width_ratios=[1.4, 1, 1], wspace=0.45))

# Age raincloud
ax = axes[0]
raincloud(ax, clin, 'response_bin', 'age', order=['good','bad'], palette=PAL_RESP)
p_age = stats.mannwhitneyu(clin[clin.response_bin=='good'].age, clin[clin.response_bin=='bad'].age).pvalue
stat_bracket(ax, 0, 1, clin.age.max()+2, p_age)
ax.set_ylabel('Age (years)'); ax.set_xlabel('')
ax.set_title('Age distribution')
annotate_count(ax, clin, 'response_bin', 'age', ['good','bad'])

# cT stacked
ax = axes[1]
ct_tab = pd.crosstab(clin.cT, clin.response_bin, normalize='columns')*100
ct_tab = ct_tab.reindex(['T2','T2/T3','T3','T4']).fillna(0)[['good','bad']]
bottom = np.zeros(2)
for i, t in enumerate(['T2','T2/T3','T3','T4']):
    ax.bar(range(2), ct_tab.loc[t].values, bottom=bottom, label=t,
           color=PAL_STAGE[t], edgecolor='white', linewidth=1.5, width=0.6)
    # label inside bar if >10%
    for j in range(2):
        if ct_tab.loc[t].values[j] > 8:
            ax.text(j, bottom[j]+ct_tab.loc[t].values[j]/2, f"{ct_tab.loc[t].values[j]:.0f}%",
                    ha='center', va='center', fontsize=8.5, color='white', fontweight='bold')
    bottom += ct_tab.loc[t].values
ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'])
ax.set_ylabel('Percentage (%)'); ax.set_ylim(0,102)
ax.legend(title='cT stage', loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=8, title_fontsize=8)
p_ct = stats.chi2_contingency(pd.crosstab(clin.response_bin, clin.cT)).pvalue
ax.set_title(f'Clinical T-stage  (χ² p = {p_ct:.3f})')

# Sex breakdown
ax = axes[2]
sex_tab = pd.crosstab(clin.sex, clin.response_bin)
x_labels = ['good','bad']
for i, sex in enumerate(['M','F']):
    vals = sex_tab.loc[sex].values if sex in sex_tab.index else [0,0]
    ax.bar([j + (i-0.5)*0.35 for j in range(2)], vals, width=0.32,
           color='#264653' if sex=='M' else '#e76f51',
           edgecolor='white', linewidth=1.2, label=sex)
    for j, v in enumerate(vals):
        ax.text(j + (i-0.5)*0.35, v+0.3, str(v), ha='center', fontsize=9, fontweight='bold')
ax.set_xticks(range(2)); ax.set_xticklabels(x_labels)
ax.set_ylabel('Patient count')
ax.legend(title='Sex', fontsize=8, loc='upper right')
p_sex = stats.fisher_exact(sex_tab.values).pvalue
ax.set_title(f'Sex  (Fisher p = {p_sex:.3f})')

save_panel(fig, 'Fig1C_clinical', OUT)

# ==========================================================
# 1D — Sample availability heatmap with rich top annotation
# ==========================================================
fig = plt.figure(figsize=(14.5, 6.5))
gs = fig.add_gridspec(8, 1, height_ratios=[0.35,0.35,0.35,0.35,0.35,0.35,0.2,5], hspace=0.15)
# Top 6 annotation tracks + main matrix
subs = sorted(clin.subject_id)
order_subj = clin.sort_values(['response_num','subject_id']).subject_id.tolist()

def ann_row(ax, values, label, palette_map=None, cmap=None, vmin=None, vmax=None,
            is_cat=True, na_color='#ecf0f1'):
    """Draw one annotation track."""
    if is_cat:
        colors_arr = np.array([[matplotlib.colors.to_rgb(palette_map.get(v, na_color)) for v in values]])
    else:
        # continuous
        values_num = pd.to_numeric(pd.Series(values), errors='coerce')
        if vmin is None: vmin = np.nanmin(values_num)
        if vmax is None: vmax = np.nanmax(values_num)
        values_num = values_num.fillna(vmin)
        norm = (values_num - vmin)/(vmax-vmin+1e-9)
        colors_arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(colors_arr, aspect='auto', interpolation='nearest')
    ax.set_yticks([0]); ax.set_yticklabels([label], fontsize=8.5)
    ax.set_xticks([]); ax.tick_params(length=0)
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)

# Gather annotation data
ann_data = pd.DataFrame({'subject_id': order_subj})
ann_data = ann_data.merge(clin[['subject_id','response_bin','cT','sex','age']], on='subject_id')
# TMB per subject (pre-treatment)
tmb_pre = tmb[tmb.timepoint=='pre'][['subject_id','TMB_nonsyn_per_Mb']].drop_duplicates('subject_id')
ann_data = ann_data.merge(tmb_pre, on='subject_id', how='left')
# MSI per subject
msi_pre = msi[msi.timepoint=='pre'][['subject_id','MSI_pct']].drop_duplicates('subject_id')
ann_data = ann_data.merge(msi_pre, on='subject_id', how='left')

ax0 = fig.add_subplot(gs[0]); ann_row(ax0, ann_data.response_bin.tolist(), 'Response', PAL_RESP)
ax1 = fig.add_subplot(gs[1]); ann_row(ax1, ann_data.cT.tolist(), 'cT stage', PAL_STAGE)
ax2 = fig.add_subplot(gs[2]); ann_row(ax2, ann_data.sex.tolist(), 'Sex', {'M':'#264653','F':'#e63946'})
ax3 = fig.add_subplot(gs[3]); ann_row(ax3, ann_data.age.tolist(), 'Age', cmap=plt.cm.YlOrRd, is_cat=False, vmin=30, vmax=80)
ax4 = fig.add_subplot(gs[4]); ann_row(ax4, ann_data.TMB_nonsyn_per_Mb.tolist(), 'TMB/Mb', cmap=plt.cm.Purples, is_cat=False, vmin=0, vmax=3)
ax5 = fig.add_subplot(gs[5]); ann_row(ax5, ann_data.MSI_pct.tolist(), 'MSI %', cmap=plt.cm.Greens, is_cat=False, vmin=0, vmax=0.3)

# spacer
# main sample matrix
ax_main = fig.add_subplot(gs[7])
mat = np.zeros((6, len(order_subj)))
for j, s in enumerate(order_subj):
    w = wes_inv[wes_inv.subject_id==s]
    r = rna_inv[rna_inv.subject_id==s]
    mat[0,j] = 1 if (w.timepoint=='normal').any() else 0
    mat[1,j] = 1 if (w.timepoint=='pre').any() else 0
    mat[2,j] = 1 if (w.timepoint=='post').any() else 0
    mat[3,j] = 1 if (r.timepoint=='normal').any() else 0
    mat[4,j] = 1 if (r.timepoint=='pre').any() else 0
    mat[5,j] = 1 if (r.timepoint=='post').any() else 0
# colored by sample type
color_mat = np.ones((6, len(order_subj), 3))
sample_colors = ['#457b9d','#1d3557','#a8dadc','#f4a261','#e76f51','#ffd6a5']
for i in range(6):
    for j in range(len(order_subj)):
        if mat[i,j]:
            color_mat[i,j] = matplotlib.colors.to_rgb(sample_colors[i])
        else:
            color_mat[i,j] = matplotlib.colors.to_rgb('#ecf0f1')
ax_main.imshow(color_mat, aspect='auto', interpolation='nearest')
ax_main.set_yticks(range(6))
ax_main.set_yticklabels(['WES normal','WES pre','WES post','RNA normal','RNA pre','RNA post'], fontsize=8.5)
ax_main.set_xticks(range(len(order_subj)))
ax_main.set_xticklabels([f'S{s}' for s in order_subj], fontsize=6.5, rotation=90)
ax_main.tick_params(length=0)
for s in ['top','right','left','bottom']: ax_main.spines[s].set_visible(False)

# Title across the whole figure
fig.suptitle('Sample availability per subject with clinical & molecular annotation tiers',
             fontsize=12, fontweight='bold', y=0.94, color='#1d3557')
save_panel(fig, 'Fig1D_sample_matrix', OUT)

# ==========================================================
# 1E — Study design (CONSORT-style, redesigned)
# ==========================================================
fig, ax = plt.subplots(figsize=(7, 6))
ax.axis('off')

def flow_box(x, y, w, h, text, color, textcolor='white', bold=True):
    ax.add_patch(mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.02",
                 linewidth=1.3, edgecolor='#1d3557', facecolor=color))
    ax.text(x, y, text, ha='center', va='center', fontsize=9, color=textcolor,
            fontweight='bold' if bold else 'normal')

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', lw=1.3, color='#1d3557'))

# Step 1: enrolled
flow_box(0.5, 0.95, 0.42, 0.08, 'LARC patients enrolled (N = 35)', BLACK)
# Step 2: response split
flow_box(0.22, 0.78, 0.30, 0.08, 'Good responders\nTRG 0-1  (n = 18)', GOOD)
flow_box(0.78, 0.78, 0.30, 0.08, 'Poor responders\nTRG 2-3  (n = 17)', BAD)
arrow(0.45, 0.91, 0.30, 0.82); arrow(0.55, 0.91, 0.70, 0.82)
# Step 3: sampling
flow_box(0.22, 0.60, 0.30, 0.09, 'WES 38 samples\nRNA-seq 33 samples', '#457b9d')
flow_box(0.78, 0.60, 0.30, 0.09, 'WES 39 samples\nRNA-seq 23 samples', '#457b9d')
arrow(0.22, 0.74, 0.22, 0.65); arrow(0.78, 0.74, 0.78, 0.65)
# Step 4: matched structure
flow_box(0.5, 0.42, 0.55, 0.08, 'Matched pre/post pairs (n = 14)\nUnmatched (n = 21)', '#8d99ae', 'white')
arrow(0.22, 0.55, 0.42, 0.46); arrow(0.78, 0.55, 0.58, 0.46)
# Step 5: analyses
flow_box(0.5, 0.22, 0.9, 0.16,
         'Integrated analyses:\nWES somatic · SBS signatures · MSI · CNV · HLA class I · HLA LOH · neoantigens\nRNA DEG · GSEA · ssGSEA · 22 immune signatures · TCR/BCR · CMS\nClonal evolution · ML LOOCV predictor · 7-cohort external meta-analysis',
         '#1d3557', 'white', bold=False)
arrow(0.5, 0.38, 0.5, 0.30)
# Step 6: findings
flow_box(0.5, 0.04, 0.9, 0.06,
         'Primary finding: DNA-repair / cell-cycle programs predict response; treatment cascade in responders only.',
         HIGHLIGHT, '#1d3557')
arrow(0.5, 0.14, 0.5, 0.07)

ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title('Study design and analytical workflow', fontsize=12, fontweight='bold', pad=8, color='#1d3557')
save_panel(fig, 'Fig1E_design', OUT)

print('\n=== Fig 1 (5 panels) complete in', OUT, '===')
