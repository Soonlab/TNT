"""
Fig 3 v3.3 — fix specific overlap issues:
  3B: y-tick gene labels manually placed left of module strip (no overlap)
      bottom legend properly positioned (no bbox_to_anchor pushing into heatmap)
      colorbar / module legend with clear gaps
  3E: revert to v3.1 (smaller padding, lower-left legend, quadrant labels at 0.85 max)
  3F: TLS sig label brought close to bar (use ylabel with labelpad)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap, to_rgb

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}
PAL_STAGE_DEEP = {'T2':'#7fb0c4', 'T2/T3':'#1c5d7e', 'T3':'#0e2a47', 'T4':'#a01b2b'}
CAT_COLORS_DEEP = {
    'Cytotoxic T-cell':'#057a64','Antigen presentation':'#00567d','IFN response':'#0099b8',
    'B-cell / TLS':'#d4a300','Innate':'#d96125','Regulatory':'#7a3aad',
    'Stromal/EMT':'#b03219','Hypoxia':'#0e2a47','Other':'#5a6772',
}
SIG_CATEGORY = {
    'CD8_proliferation':'Cytotoxic T-cell','CD8_activation':'Cytotoxic T-cell',
    'CD8_exhaustion':'Cytotoxic T-cell','Cytolytic_activity':'Cytotoxic T-cell',
    'Antigen_presentation':'Antigen presentation','MHC_II':'Antigen presentation',
    'NLRC5_HLA_IFNG':'Antigen presentation','IFNg_Ayers_18':'IFN response',
    'TLS_Cabrita':'B-cell / TLS','B_cell':'B-cell / TLS',
    'NK_cell':'Innate','Mac_M1':'Innate','Mac_M2':'Innate',
    'Treg':'Regulatory','Checkpoint_inhibitory':'Regulatory',
    'TGFb_Mariathasan':'Stromal/EMT','EMT_Mak':'Stromal/EMT',
    'CAF_iCAF':'Stromal/EMT','CAF_myCAF':'Stromal/EMT',
    'Hypoxia_Buffa':'Hypoxia','Stemness_mRNAsi_proxy':'Other','Epithelial':'Other',
}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'
clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
tpm = pd.read_csv(ROOT/'06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index(); sigs_m.columns = ['sample_id'] + list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

# ============================================================
# FIG 3B — clean layout, no overlaps
# ============================================================
sig_cols = [c for c in SIG_CATEGORY if c in sigs_m.columns]
pre = sigs_m[sigs_m.timepoint=='pre']
mat = pre.set_index('sample_id')[sig_cols]
mat_z = mat.sub(mat.mean()).div(mat.std())
order_resp = pre.sort_values('response_bin').sample_id.tolist()
mat_z = mat_z.loc[order_resp]
sig_cols_sorted = sorted(sig_cols, key=lambda s: (SIG_CATEGORY.get(s,'Other'), s))

def s2subj(s): return rna_inv[rna_inv.sample_id==s].subject_id.iloc[0]
resp_a = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_resp], index=order_resp)
ct_a = pd.Series([clin[clin.subject_id==s2subj(s)].cT.iloc[0] for s in order_resp], index=order_resp)
sex_a = pd.Series([clin[clin.subject_id==s2subj(s)].sex.iloc[0] for s in order_resp], index=order_resp)
age_a = pd.Series([clin[clin.subject_id==s2subj(s)].age.iloc[0] for s in order_resp], index=order_resp)
tmb_map = tmb.set_index('sample_id')['TMB_nonsyn_per_Mb']
tmb_a = pd.Series([tmb_map.get(f'{s2subj(s)}-PR', np.nan) for s in order_resp], index=order_resp)

fig = plt.figure(figsize=(17, 9))
# Layout columns:
#   col 0: blank pad (for ann labels + sig labels)
#   col 1: ann labels (text only, no bg)
#   col 2: module color strip (only at row 5)
#   col 3: HEATMAP (also annotation strips above)
#   col 4: gap
#   col 5: colorbar
#   col 6: gap
#   col 7: module legend
gs = fig.add_gridspec(7, 8,
    height_ratios=[0.16, 0.16, 0.16, 0.16, 0.16, 5.4, 1.0],
    width_ratios=[0.05, 1.6, 0.18, 6.0, 0.4, 0.20, 0.4, 1.7],
    hspace=0.06, wspace=0.0)

heat = mat_z[sig_cols_sorted].T
n_samples = heat.shape[1]; n_sigs = heat.shape[0]

def ann_row_strip(ax, vals, palette_map=None, cmap=None, vmin=None, vmax=None):
    if palette_map:
        arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        v_n = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(v_n)
        if vmax is None: vmax = np.nanmax(v_n)
        norm = ((v_n-vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[0, n_samples, 0, 1])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, n_samples)

ann_data_list = [
    ('Response',  resp_a.values.tolist(), {'palette_map':PAL}),
    ('cT stage',  ct_a.values.tolist(),   {'palette_map':PAL_STAGE_DEEP}),
    ('Sex',       sex_a.values.tolist(),  {'palette_map':{'M':'#0e2a47','F':'#a01b2b'}}),
    ('Age',       age_a.values.tolist(),  {'cmap':plt.cm.YlOrRd, 'vmin':30, 'vmax':80}),
    ('TMB / Mb',  tmb_a.values.tolist(),  {'cmap':plt.cm.Purples, 'vmin':0, 'vmax':3}),
]
for i, (lbl, vals, kwargs) in enumerate(ann_data_list):
    # Annotation strip (col 3, the heatmap col)
    ax_a = fig.add_subplot(gs[i, 3])
    ann_row_strip(ax_a, vals, **kwargs)
    # Annotation label (col 1, text only, right-aligned)
    ax_lbl = fig.add_subplot(gs[i, 1])
    ax_lbl.axis('off')
    ax_lbl.text(1.0, 0.5, lbl, ha='right', va='center', fontsize=10.5,
                color='#0e2a47', fontweight='bold', transform=ax_lbl.transAxes)

# Module color strip (col 2, row 5) — directly left of heatmap
ax_mod = fig.add_subplot(gs[5, 2])
row_colors = [CAT_COLORS_DEEP[SIG_CATEGORY[s]] for s in sig_cols_sorted]
arr_rc = np.array([[to_rgb(c) for c in row_colors]])
ax_mod.imshow(arr_rc.transpose(1,0,2), aspect='auto', interpolation='nearest',
              extent=[0, 1, 0, n_sigs])
ax_mod.set_xticks([]); ax_mod.set_yticks([])
for s in ['top','right','left','bottom']: ax_mod.spines[s].set_visible(False)
ax_mod.invert_yaxis()

# Signature labels in col 1 (LEFT of module strip — no overlap with module color)
ax_sig_lbl = fig.add_subplot(gs[5, 1])
ax_sig_lbl.axis('off')
ax_sig_lbl.set_xlim(0, 1); ax_sig_lbl.set_ylim(0, n_sigs)
ax_sig_lbl.invert_yaxis()
for i, s in enumerate(sig_cols_sorted):
    ax_sig_lbl.text(1.0, i+0.5, s.replace('_',' '), ha='right', va='center',
                    fontsize=10, color='#0e2a47')

# Main heatmap (col 3, row 5)
ax_h = fig.add_subplot(gs[5, 3])
im = ax_h.imshow(heat.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                 aspect='auto', interpolation='nearest',
                 extent=[0, n_samples, 0, n_sigs])
ax_h.set_xticks([]); ax_h.set_yticks([])
ax_h.invert_yaxis()
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# Colorbar (col 5, rows 0-2 i.e., top-right area)
cax = fig.add_subplot(gs[0:3, 5])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('z-score', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.8)

# Module legend (col 7, top portion)
ax_modleg = fig.add_subplot(gs[0:5, 7])
ax_modleg.axis('off')
ax_modleg.text(0.0, 0.99, 'Signature module', fontsize=10.5, color='#0e2a47',
               fontweight='bold', ha='left', va='top', transform=ax_modleg.transAxes)
mod_present = []
for s in sig_cols_sorted:
    cat = SIG_CATEGORY[s]
    if cat not in mod_present: mod_present.append(cat)
y_step = 0.085
for i, mod in enumerate(mod_present):
    y = 0.92 - i*y_step
    ax_modleg.add_patch(Rectangle((0.0, y-0.025), 0.18, 0.05, color=CAT_COLORS_DEEP[mod],
                                   transform=ax_modleg.transAxes))
    ax_modleg.text(0.22, y, mod, fontsize=9.5, va='center',
                   transform=ax_modleg.transAxes, color='#0e2a47')

# Bottom clinical legend (row 6, full width col 1-3)
ax_clin = fig.add_subplot(gs[6, 1:4])
ax_clin.axis('off')
clin_handles = []
for k,c in PAL.items(): clin_handles.append(mpatches.Patch(color=c, label=k))
for k,c in PAL_STAGE_DEEP.items(): clin_handles.append(mpatches.Patch(color=c, label=k))
clin_handles.append(mpatches.Patch(color='#0e2a47', label='Male'))
clin_handles.append(mpatches.Patch(color='#a01b2b', label='Female'))
clin_handles.append(mpatches.Patch(color=plt.cm.YlOrRd(0.7), label='Age (gradient)'))
clin_handles.append(mpatches.Patch(color=plt.cm.Purples(0.7), label='TMB (gradient)'))
ax_clin.legend(handles=clin_handles, loc='center', fontsize=9.5, ncol=6,
               frameon=False)

save_panel(fig, 'Fig3B_signature_heatmap', OUT)

# ============================================================
# FIG 3E — revert to v3.1 style
# ============================================================
fig = plt.figure(figsize=(7.5, 7.5))
gs = gridspec.GridSpec(4, 4, hspace=0.05, wspace=0.05)
ax_main = fig.add_subplot(gs[1:, :3])
ax_top = fig.add_subplot(gs[0, :3], sharex=ax_main)
ax_right = fig.add_subplot(gs[1:, 3], sharey=ax_main)

x_sig = 'Cytolytic_activity'; y_sig = 'CD8_exhaustion'
pre = sigs_m[sigs_m.timepoint=='pre']

for resp in ['good','bad']:
    sub = pre[pre.response_bin==resp]
    ax_main.scatter(sub[x_sig], sub[y_sig], color=PAL[resp], s=55, alpha=0.92,
                    edgecolor='white', linewidth=0.9, label=f'{resp} (n={len(sub)})', zorder=4)
    try:
        sns.kdeplot(x=sub[x_sig], y=sub[y_sig], ax=ax_main, color=PAL[resp],
                    levels=3, alpha=0.45, linewidths=1.0)
    except: pass
    sns.kdeplot(x=sub[x_sig].dropna(), ax=ax_top, color=PAL[resp], fill=True, alpha=0.45, lw=1.6)
    sns.kdeplot(y=sub[y_sig].dropna(), ax=ax_right, color=PAL[resp], fill=True, alpha=0.45, lw=1.6)

xmid = pre[x_sig].median(); ymid = pre[y_sig].median()
ax_main.axvline(xmid, color='#5a6772', ls=':', lw=0.7, alpha=0.7)
ax_main.axhline(ymid, color='#5a6772', ls=':', lw=0.7, alpha=0.7)
xmax_, ymax_ = pre[x_sig].max(), pre[y_sig].max()
xmin_, ymin_ = pre[x_sig].min(), pre[y_sig].min()
quad_labels = [
    (xmax_*0.85, ymax_*0.85, 'Active &\nExhausted', '#0e2a47'),
    (xmin_*0.85, ymax_*0.85, 'Cold &\nExhausted', '#0e2a47'),
    (xmax_*0.85, ymin_*0.85, 'Functional\ncytotoxic', GOOD_DEEP),
    (xmin_*0.85, ymin_*0.85, 'Cold tumor', BAD_DEEP),
]
for x, y, lbl, color in quad_labels:
    ax_main.text(x, y, lbl, ha='center', va='center', fontsize=8.5, color=color,
                 fontweight='bold', alpha=0.85,
                 bbox=dict(facecolor='white', edgecolor=color, alpha=0.85, boxstyle='round,pad=0.3'))

ax_main.set_xlabel('CD8 cytolytic activity (z-score)', fontsize=11, fontweight='bold', color='#0e2a47')
ax_main.set_ylabel('CD8 exhaustion (z-score)', fontsize=11, fontweight='bold', color='#0e2a47')
ax_main.legend(loc='lower left', fontsize=10, frameon=False)

ax_top.set_xticks([]); ax_top.set_yticks([])
ax_right.set_xticks([]); ax_right.set_yticks([])
for s in ['top','right','bottom']: ax_top.spines[s].set_visible(False)
for s in ['top','right','left']: ax_right.spines[s].set_visible(False)
add_axis_spines(ax_main)
save_panel(fig, 'Fig3E_CD8_biaxial', OUT)

# ============================================================
# FIG 3F — TLS sig label close to bar
# ============================================================
TLS_GENES_CABRITA = ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3','CD79B','MS4A1','CCL18','PTGDS','CXCL8']
present_tls = [g for g in TLS_GENES_CABRITA if g in log_tpm.index]
pre = sigs_m[sigs_m.timepoint=='pre'].sort_values(['response_bin','TLS_Cabrita'], ascending=[True, False])
pre_samples = pre.sample_id.tolist()
tls_mat = log_tpm.loc[present_tls, pre_samples]
tls_z = tls_mat.sub(tls_mat.mean(axis=1), axis=0).div(tls_mat.std(axis=1), axis=0)
n_s = len(pre_samples); n_g = len(present_tls)

fig = plt.figure(figsize=(14, 7))
# Width: [annot label space: 1.4] [bar/heatmap: 14] [colorbar: 0.4]
gs = fig.add_gridspec(4, 3, height_ratios=[0.7, 0.18, 0.18, 4.5],
                      width_ratios=[1.4, 14, 0.4], hspace=0.10, wspace=0.04)

# TOP bar (col 1)
ax_top = fig.add_subplot(gs[0, 1])
tls_sig_vals = pre.TLS_Cabrita.values
colors = [PAL[r] for r in pre.response_bin]
ax_top.bar(np.arange(n_s), tls_sig_vals, color=colors, edgecolor='white', linewidth=0.5, width=0.92)
ax_top.axhline(0, color='#0e2a47', lw=0.6)
# Use ylabel directly (not transAxes external text)
ax_top.set_ylabel('TLS sig.\n(z-score)', fontsize=10.5, fontweight='bold',
                  color='#0e2a47', rotation=90, labelpad=8)
ax_top.set_xticks([]); ax_top.set_xlim(-0.5, n_s-0.5)
add_axis_spines(ax_top)
ax_top.tick_params(labelsize=8.5)

# Stats inside top-right
g = pre[pre.response_bin=='good'].TLS_Cabrita
b = pre[pre.response_bin=='bad'].TLS_Cabrita
p_tls = stats.mannwhitneyu(g, b).pvalue
ax_top.text(0.99, 0.95, f'TLS sig: good vs poor   p = {p_tls:.3f}',
            transform=ax_top.transAxes, ha='right', va='top',
            fontsize=10, color='#0e2a47', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.92, boxstyle='round,pad=0.4'))

# Annotation strips
def ann_strip3(ax, vals, label, palette_map=None, cmap=None, vmin=None, vmax=None):
    if palette_map:
        arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        v_n = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(v_n)
        if vmax is None: vmax = np.nanmax(v_n)
        norm = ((v_n-vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[0, n_s, 0, 1])
    ax.set_yticks([]); ax.set_xticks([])
    ax.set_ylabel(label, rotation=0, labelpad=8, fontsize=10, fontweight='bold',
                  color='#0e2a47', va='center', ha='right')
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, n_s)

ax_a1 = fig.add_subplot(gs[1, 1])
ann_strip3(ax_a1, pre.response_bin.values.tolist(), 'Response', palette_map=PAL)
ax_a2 = fig.add_subplot(gs[2, 1])
ann_strip3(ax_a2, pre.MHC_II.values.tolist(), 'MHC-II z',
           cmap=plt.cm.Purples, vmin=pre.MHC_II.min(), vmax=pre.MHC_II.max())

# Main heatmap (col 1, row 3)
ax_main = fig.add_subplot(gs[3, 1])
im = ax_main.imshow(tls_z.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                     aspect='auto', interpolation='nearest',
                     extent=[0, n_s, 0, n_g])
ax_main.set_yticks(np.arange(n_g)+0.5)
ax_main.set_yticklabels(present_tls[::-1], fontsize=11, fontstyle='italic', color='#0e2a47')
ax_main.invert_yaxis()
xt = np.arange(0, n_s, 2)
ax_main.set_xticks(xt+0.5)
ax_main.set_xticklabels([pre_samples[i] for i in xt], fontsize=8, rotation=90, color='#0e2a47')
ax_main.tick_params(length=2)
ax_main.set_xlim(0, n_s)
for s in ['top','right','left','bottom']: ax_main.spines[s].set_visible(False)

# Colorbar col 2 (only for heatmap row 3, not full height)
cax = fig.add_subplot(gs[3, 2])
cb = fig.colorbar(im, cax=cax)
cb.set_label('log2(TPM+1)\nz-score', fontsize=9.5, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)

save_panel(fig, 'Fig3F_TLS_Cabrita', OUT)

print('\n=== Fig 3 v3.3 (3B/3E/3F) saved ===')
