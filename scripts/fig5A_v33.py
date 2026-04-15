"""
Fig 5A v3.3 — omic group legend moved BELOW colorbar (same column, stacked).
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb
from scipy.cluster.hierarchy import linkage, leaves_list
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
OMIC_COLORS = {
    'Clinical':'#0e2a47','WES_genomic':'#118ab2','WES_CNV':'#06aed5',
    'RNA_pathway':'#057a64','RNA_immune':'#c11456','RNA_stromal':'#b03219','Other':'#5a6772',
}
def assign_omic(f):
    if f in ['TMB_nonsyn_per_Mb','n_nonsyn','MMR_prop','SBS5','SBS3','MSI_pct']: return 'WES_genomic'
    if f in ['CIN','frac_amp','frac_del']: return 'WES_CNV'
    if 'CD8' in f or 'MHC' in f or 'IFN' in f or 'TLS' in f or 'NLRC' in f or 'NK_' in f or 'B_' in f or 'Mac_' in f or 'Treg' in f: return 'RNA_immune'
    if 'EMT' in f or 'TGF' in f or 'CAF' in f or 'Hypox' in f or 'Stem' in f or 'Epith' in f: return 'RNA_stromal'
    if 'Repair' in f or 'HDR' in f or 'E2F' in f or 'G2-M' in f or 'Myc' in f or 'Hallmark' in f: return 'RNA_pathway'
    if f in ['age']: return 'Clinical'
    return 'Other'

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'
integ = pd.read_csv(ROOT/'tables/integrated_subject_master.tsv', sep='\t')
feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
omic_map = {f: assign_omic(f) for f in feats_num}

X_raw = integ[feats_num].apply(pd.to_numeric, errors='coerce')
X_imputed = X_raw.fillna(X_raw.median())
corr = X_imputed.corr(method='spearman')

omic_order = ['Clinical','WES_genomic','WES_CNV','RNA_pathway','RNA_immune','RNA_stromal','Other']
ordered_features = []
for og in omic_order:
    in_og = [f for f in feats_num if omic_map[f]==og]
    if len(in_og) > 1:
        sub_corr = corr.loc[in_og, in_og].fillna(0)
        Z = linkage(sub_corr.values, method='average')
        idx = leaves_list(Z)
        ordered_features.extend([in_og[i] for i in idx])
    else:
        ordered_features.extend(in_og)
corr_o = corr.loc[ordered_features, ordered_features]
N = len(ordered_features)

# ============================================================
# Layout: cols = [y-label 2.0, left bar 0.14, heatmap 8, gap 0.4, right column 1.9]
# Right column: sub-gridspec (colorbar top + legend below)
# ============================================================
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(4, 5,
    height_ratios=[0.12, 0.14, 8, 1.8],
    width_ratios=[2.0, 0.14, 8, 0.4, 1.9],
    hspace=0.04, wspace=0.04)

top_colors = [OMIC_COLORS[omic_map[f]] for f in ordered_features]

# Top color bar
ax_top = fig.add_subplot(gs[1, 2])
arr_top = np.array([[to_rgb(c) for c in top_colors]])
ax_top.imshow(arr_top, aspect='auto', interpolation='nearest', extent=[0, N, 0, 1])
ax_top.set_xticks([]); ax_top.set_yticks([])
for s in ['top','right','left','bottom']: ax_top.spines[s].set_visible(False)
ax_top.set_xlim(0, N)

# Left color bar
ax_left = fig.add_subplot(gs[2, 1])
left_arr = np.array([[to_rgb(c) for c in top_colors]])
ax_left.imshow(left_arr.transpose(1,0,2), aspect='auto', interpolation='nearest', extent=[0, 1, 0, N])
ax_left.set_xticks([]); ax_left.set_yticks([])
for s in ['top','right','left','bottom']: ax_left.spines[s].set_visible(False)
ax_left.invert_yaxis()

# Y-label column
ax_ylbl = fig.add_subplot(gs[2, 0])
ax_ylbl.axis('off')
ax_ylbl.set_xlim(0, 1); ax_ylbl.set_ylim(0, N)
ax_ylbl.invert_yaxis()
for i, f in enumerate(ordered_features):
    ax_ylbl.text(0.99, i+0.5, f[:28], ha='right', va='center', fontsize=7.5, color='#0e2a47')

# Main heatmap
ax_h = fig.add_subplot(gs[2, 2])
im = ax_h.imshow(corr_o.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto',
                 interpolation='nearest', extent=[0, N, 0, N])
ax_h.set_xticks([]); ax_h.set_yticks([])
ax_h.invert_yaxis()
ax_h.set_xlim(0, N); ax_h.set_ylim(N, 0)
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# X-label row
ax_xlbl = fig.add_subplot(gs[3, 2])
ax_xlbl.axis('off')
ax_xlbl.set_xlim(0, N); ax_xlbl.set_ylim(0, 1)
for i, f in enumerate(ordered_features):
    ax_xlbl.text(i+0.5, 0.98, f[:28], ha='right', va='top', rotation=90,
                 fontsize=7.5, color='#0e2a47')

# Group dividers
cur = 0
group_starts = []
for og in omic_order:
    in_og_count = sum(1 for f in ordered_features if omic_map[f]==og)
    if in_og_count == 0: continue
    group_starts.append(cur)
    cur += in_og_count
for gs_pos in group_starts[1:]:
    ax_h.axvline(gs_pos, color='white', lw=2)
    ax_h.axhline(gs_pos, color='white', lw=2)

# ============================================================
# Right column (col 4, row 2): colorbar ON TOP + omic legend BELOW
# Use sub-gridspec to split into 2 vertical zones
# ============================================================
right_gs = gs[2, 4].subgridspec(2, 1, height_ratios=[1.5, 3], hspace=0.05)

# Colorbar zone (top)
ax_cb_holder = fig.add_subplot(right_gs[0])
ax_cb_holder.axis('off')
cbar_ax = inset_axes(ax_cb_holder, width="22%", height="90%", loc='upper left',
                     bbox_to_anchor=(0.0, 0, 1, 1), bbox_transform=ax_cb_holder.transAxes)
cb = fig.colorbar(im, cax=cbar_ax, orientation='vertical')
cb.set_label('Spearman ρ', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)

# Legend zone (below colorbar)
ax_leg = fig.add_subplot(right_gs[1])
ax_leg.axis('off')
ax_leg.text(0.0, 0.98, 'Omic group', fontsize=10.5, color='#0e2a47', fontweight='bold',
            ha='left', va='top', transform=ax_leg.transAxes)
omic_present = [og for og in omic_order if any(omic_map[f]==og for f in ordered_features)]
y_step = 0.11
for i, og in enumerate(omic_present):
    y_pos = 0.88 - i*y_step
    ax_leg.add_patch(Rectangle((0.02, y_pos-0.030), 0.18, 0.060, color=OMIC_COLORS[og],
                                transform=ax_leg.transAxes))
    ax_leg.text(0.24, y_pos, og.replace('_',' '), transform=ax_leg.transAxes,
                fontsize=9.5, va='center', color='#0e2a47')

save_panel(fig, 'Fig5A_correlation', OUT)
print('=== Fig 5A v3.3 saved (legend below colorbar) ===')
