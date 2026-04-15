"""
Fig 3B v3.4 — colorbar (top-right) + signature module legend BELOW colorbar (same column).
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgb

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
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')

sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index(); sigs_m.columns = ['sample_id'] + list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

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

fig = plt.figure(figsize=(16, 9))
# Single right column for colorbar (top) + module legend (below)
# Width ratios: [pad | sig labels | module strip | heatmap | gap | right column]
gs = fig.add_gridspec(8, 6,
    height_ratios=[0.16, 0.16, 0.16, 0.16, 0.16, 0.18, 5.4, 1.0],
    width_ratios=[0.05, 1.6, 0.18, 6.0, 0.4, 1.7],
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
    ax_a = fig.add_subplot(gs[i, 3])
    ann_row_strip(ax_a, vals, **kwargs)
    ax_lbl = fig.add_subplot(gs[i, 1])
    ax_lbl.axis('off')
    ax_lbl.text(1.0, 0.5, lbl, ha='right', va='center', fontsize=10.5,
                color='#0e2a47', fontweight='bold', transform=ax_lbl.transAxes)

# Module color strip (col 2, row 6) — directly left of heatmap
ax_mod = fig.add_subplot(gs[6, 2])
row_colors = [CAT_COLORS_DEEP[SIG_CATEGORY[s]] for s in sig_cols_sorted]
arr_rc = np.array([[to_rgb(c) for c in row_colors]])
ax_mod.imshow(arr_rc.transpose(1,0,2), aspect='auto', interpolation='nearest',
              extent=[0, 1, 0, n_sigs])
ax_mod.set_xticks([]); ax_mod.set_yticks([])
for s in ['top','right','left','bottom']: ax_mod.spines[s].set_visible(False)
ax_mod.invert_yaxis()

# Signature labels in col 1 (LEFT of module strip)
ax_sig_lbl = fig.add_subplot(gs[6, 1])
ax_sig_lbl.axis('off')
ax_sig_lbl.set_xlim(0, 1); ax_sig_lbl.set_ylim(0, n_sigs)
ax_sig_lbl.invert_yaxis()
for i, s in enumerate(sig_cols_sorted):
    ax_sig_lbl.text(1.0, i+0.5, s.replace('_',' '), ha='right', va='center',
                    fontsize=10, color='#0e2a47')

# Main heatmap (col 3, row 6)
ax_h = fig.add_subplot(gs[6, 3])
im = ax_h.imshow(heat.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                 aspect='auto', interpolation='nearest',
                 extent=[0, n_samples, 0, n_sigs])
ax_h.set_xticks([]); ax_h.set_yticks([])
ax_h.invert_yaxis()
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# Right column: colorbar (TOP) + module legend (BELOW), stacked vertically
# Colorbar in rows 0-1
cax = fig.add_subplot(gs[0:2, 5])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('z-score', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.8)

# Module legend in rows 3-6 (below colorbar with one row gap)
ax_modleg = fig.add_subplot(gs[3:7, 5])
ax_modleg.axis('off')
ax_modleg.text(0.0, 0.99, 'Signature module', fontsize=10.5, color='#0e2a47',
               fontweight='bold', ha='left', va='top', transform=ax_modleg.transAxes)
mod_present = []
for s in sig_cols_sorted:
    cat = SIG_CATEGORY[s]
    if cat not in mod_present: mod_present.append(cat)
y_step = 0.10
for i, mod in enumerate(mod_present):
    y = 0.92 - i*y_step
    ax_modleg.add_patch(Rectangle((0.0, y-0.025), 0.18, 0.05, color=CAT_COLORS_DEEP[mod],
                                   transform=ax_modleg.transAxes))
    ax_modleg.text(0.22, y, mod, fontsize=9.5, va='center',
                   transform=ax_modleg.transAxes, color='#0e2a47')

# Bottom clinical legend (row 7, full width col 1-3)
ax_clin = fig.add_subplot(gs[7, 1:4])
ax_clin.axis('off')
clin_handles = []
for k,c in PAL.items(): clin_handles.append(mpatches.Patch(color=c, label=k))
for k,c in PAL_STAGE_DEEP.items(): clin_handles.append(mpatches.Patch(color=c, label=k))
clin_handles.append(mpatches.Patch(color='#0e2a47', label='Male'))
clin_handles.append(mpatches.Patch(color='#a01b2b', label='Female'))
clin_handles.append(mpatches.Patch(color=plt.cm.YlOrRd(0.7), label='Age (gradient)'))
clin_handles.append(mpatches.Patch(color=plt.cm.Purples(0.7), label='TMB (gradient)'))
ax_clin.legend(handles=clin_handles, loc='center', fontsize=9.5, ncol=6, frameon=False)

save_panel(fig, 'Fig3B_signature_heatmap', OUT)
print('=== Fig 3B v3.4 saved ===')
