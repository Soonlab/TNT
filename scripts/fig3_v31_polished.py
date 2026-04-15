"""
Figure 3 v3.1 — polished revision based on user feedback:
  - Remove all panel titles (rely on figure legend)
  - Deeper, more saturated colors
  - 3A radar: more sophisticated, no label overlap
  - 3B heatmap: colorbar moved to top-right
  - 3C volcano: arrows moved inside plot, closer to data
  - 3E biaxial: smaller point sizes
  - 3F TLS heatmap: no overlapping text
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Rectangle, Polygon, Circle
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from adjustText import adjust_text

# ---------- Saturated palette overrides ----------
GOOD_DEEP = '#0a7d6e'  # deeper teal
BAD_DEEP = '#c53e1f'   # deeper coral
BLACK_DEEP = '#0e2a47'

PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}
PAL_STAGE_DEEP = {'T2':'#7fb0c4', 'T2/T3':'#1c5d7e', 'T3':'#0e2a47', 'T4':'#a01b2b'}
CAT_COLORS_DEEP = {
    'Cytotoxic T-cell':'#057a64',
    'Antigen presentation':'#00567d',
    'IFN response':'#0099b8',
    'B-cell / TLS':'#d4a300',
    'Innate':'#d96125',
    'Regulatory':'#7a3aad',
    'Stromal/EMT':'#b03219',
    'Hypoxia':'#0e2a47',
    'Other':'#5a6772',
}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
sig_stats = pd.read_csv(ROOT/'06_rna_immune/sig_response_stats.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
tpm = pd.read_csv(ROOT/'06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)
tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')

sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index(); sigs_m.columns = ['sample_id'] + list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

UNMATCHED = [13,15,16,17,18,19,33]
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

# ============================================================
# FIG 3A — Polished TME radar (no title, no label overlap, refined aesthetics)
# ============================================================
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='polar')

radar_sigs = [
    ('CD8_proliferation','CD8 proliferation'),
    ('CD8_activation','CD8 activation'),
    ('Cytolytic_activity','Cytolytic activity'),
    ('Antigen_presentation','MHC-I antigen pres.'),
    ('MHC_II','MHC-II'),
    ('IFNg_Ayers_18','IFN-γ (Ayers)'),
    ('TLS_Cabrita','TLS (Cabrita)'),
    ('B_cell','B cell'),
    ('NK_cell','NK cell'),
    ('Treg','Treg'),
    ('TGFb_Mariathasan','TGF-β'),
    ('EMT_Mak','EMT'),
]
sig_keys = [s[0] for s in radar_sigs]
sig_labels = [s[1] for s in radar_sigs]
N = len(radar_sigs)
theta = np.linspace(0, 2*np.pi, N, endpoint=False)

pre = sigs_m[sigs_m.timepoint=='pre']
good_med = pre[pre.response_bin=='good'][sig_keys].median().values
bad_med = pre[pre.response_bin=='bad'][sig_keys].median().values

# Significance per axis
ps = []
for sk in sig_keys:
    g = pre[pre.response_bin=='good'][sk].dropna()
    b = pre[pre.response_bin=='bad'][sk].dropna()
    p = stats.mannwhitneyu(g,b).pvalue if len(g)>=3 and len(b)>=3 else 1
    ps.append(p)

theta_c = np.append(theta, theta[0])
good_c = np.append(good_med, good_med[0])
bad_c = np.append(bad_med, bad_med[0])

# Subtle radial grid (custom)
ymin = min(good_med.min(), bad_med.min()) - 0.15
ymax = max(good_med.max(), bad_med.max()) + 0.20
gridvals = np.linspace(ymin, ymax, 5)
ax.set_ylim(ymin, ymax)

# Soft dotted radial grid only at zero + max
for r in [0, ymax*0.5]:
    ax.plot(np.linspace(0, 2*np.pi, 200), [r]*200, color='#aab3bf', lw=0.7, ls=(0,(2,3)), alpha=0.65)
ax.plot(np.linspace(0, 2*np.pi, 200), [0]*200, color='#0e2a47', lw=1.1, alpha=0.7)

# Polygon fills with rich edge
ax.fill(theta_c, bad_c, color=BAD_DEEP, alpha=0.30)
ax.plot(theta_c, bad_c, color=BAD_DEEP, lw=2.6, alpha=0.95, zorder=3)
ax.fill(theta_c, good_c, color=GOOD_DEEP, alpha=0.40)
ax.plot(theta_c, good_c, color=GOOD_DEEP, lw=2.8, alpha=0.95, zorder=4)

# Larger filled markers with white edge
for tt, gv, bv in zip(theta, good_med, bad_med):
    ax.scatter(tt, gv, s=110, color=GOOD_DEEP, edgecolor='white', lw=1.4, zorder=6)
    ax.scatter(tt, bv, s=110, color=BAD_DEEP, edgecolor='white', lw=1.4, zorder=6)

# Axis labels OUTSIDE the polygon, rotated to align with spoke direction
label_r = ymax + 0.10
ax.set_xticks([])  # disable default
for tt, lbl, p in zip(theta, sig_labels, ps):
    angle_deg = np.degrees(tt)
    # Rotation: tangential to circle
    if 0 <= angle_deg <= 180:
        rot = -angle_deg
    else:
        rot = 180 - angle_deg
    # Significance star
    star = sig_symbol(p)
    star_str = f'  {star}' if (star and star != 'ns') else ''
    label_text = f'{lbl}{star_str}'
    ax.text(tt, label_r, label_text, ha='center', va='center',
            fontsize=10.5, color='#0e2a47', fontweight='bold',
            rotation=rot, rotation_mode='anchor')

# Tick numbers - small, only every other
ax.set_yticks([-0.5, 0, 0.5])
ax.set_yticklabels(['−0.5','0','+0.5'], fontsize=8.5, color='#5a6772')
ax.tick_params(pad=0)

ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
ax.spines['polar'].set_color('#0e2a47'); ax.spines['polar'].set_linewidth(1.2)
ax.grid(False)

# Custom legend (in-figure, lower right)
legend_x, legend_y = 0.88, 0.05
fig.text(legend_x, legend_y+0.05, 'Good responders', color=GOOD_DEEP, fontsize=10.5, fontweight='bold', ha='right')
fig.text(legend_x, legend_y, 'Poor responders',  color=BAD_DEEP, fontsize=10.5, fontweight='bold', ha='right')

save_panel(fig, 'Fig3A_TME_radar', OUT)

# ============================================================
# FIG 3B — Heatmap (no title; colorbar top-right aligned with module legend)
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

# Build figure with custom layout: leave top-right for colorbar + module legend
fig = plt.figure(figsize=(15, 8))
gs = fig.add_gridspec(7, 8,
    height_ratios=[0.16, 0.16, 0.16, 0.16, 0.16, 5.0, 0.6],  # 5 ann strips + main + legend bottom
    width_ratios=[0.20, 5.5, 0.20, 0.20, 0.20, 0.20, 0.50, 1.4],
    hspace=0.04, wspace=0.06)

# Compute z matrix (rows=signatures, cols=samples)
heat = mat_z[sig_cols_sorted].T
# Annotation rows go ABOVE main heatmap
def ann_row(ax, vals, label, palette_map=None, cmap=None, vmin=None, vmax=None, label_color='#0e2a47'):
    if palette_map:
        arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    else:
        v_n = pd.to_numeric(pd.Series(vals), errors='coerce')
        if vmin is None: vmin = np.nanmin(v_n)
        if vmax is None: vmax = np.nanmax(v_n)
        norm = ((v_n-vmin)/(vmax-vmin+1e-9)).clip(0,1).fillna(0)
        arr = cmap(norm.values)[:,:3][np.newaxis,...]
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[0, heat.shape[1], 0, 1])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, heat.shape[1])
    # Label on the LEFT (col 0)

# Annotation strips (rows 0-4, col 1)
ax_a1 = fig.add_subplot(gs[0,1]); ann_row(ax_a1, resp_a.values.tolist(), 'Response', PAL)
ax_a2 = fig.add_subplot(gs[1,1]); ann_row(ax_a2, ct_a.values.tolist(), 'cT stage', PAL_STAGE_DEEP)
ax_a3 = fig.add_subplot(gs[2,1]); ann_row(ax_a3, sex_a.values.tolist(), 'Sex', {'M':'#0e2a47','F':'#a01b2b'})
ax_a4 = fig.add_subplot(gs[3,1]); ann_row(ax_a4, age_a.values.tolist(), 'Age', cmap=plt.cm.YlOrRd, vmin=30, vmax=80)
ax_a5 = fig.add_subplot(gs[4,1]); ann_row(ax_a5, tmb_a.values.tolist(), 'TMB/Mb', cmap=plt.cm.Purples, vmin=0, vmax=3)
# Labels on the LEFT side (col 0)
for ax_, lbl in zip([ax_a1,ax_a2,ax_a3,ax_a4,ax_a5], ['Response','cT stage','Sex','Age','TMB/Mb']):
    ax_.text(-0.5, 0.5, lbl, ha='right', va='center', fontsize=9.5, color='#0e2a47', fontweight='bold')

# Row colors strip (col 0, row 5): module category
ax_rc = fig.add_subplot(gs[5,0])
row_colors = [CAT_COLORS_DEEP[SIG_CATEGORY[s]] for s in sig_cols_sorted]
arr_rc = np.array([[to_rgb(c) for c in row_colors]])
ax_rc.imshow(arr_rc.transpose(1,0,2), aspect='auto', interpolation='nearest',
             extent=[0,1,0,len(sig_cols_sorted)])
ax_rc.set_xticks([]); ax_rc.set_yticks([])
for s in ['top','right','left','bottom']: ax_rc.spines[s].set_visible(False)
ax_rc.set_ylim(0, len(sig_cols_sorted)); ax_rc.invert_yaxis()

# Main heatmap (row 5, col 1)
ax_h = fig.add_subplot(gs[5,1])
im = ax_h.imshow(heat.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5, aspect='auto', interpolation='nearest',
                 extent=[0, heat.shape[1], 0, heat.shape[0]])
ax_h.set_xticks([]); ax_h.set_yticks(np.arange(heat.shape[0])+0.5)
ax_h.set_yticklabels([s.replace('_',' ') for s in sig_cols_sorted[::-1]], fontsize=9)
ax_h.invert_yaxis()
ax_h.set_xlim(0, heat.shape[1])
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)
ax_h.tick_params(length=0)

# Top-right area: colorbar (z-score) + module legend, side-by-side
# Colorbar in gs[0:3, 7]
cax = fig.add_subplot(gs[0:3, 7])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('z-score', fontsize=9.5, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8)
cb.outline.set_edgecolor('#0e2a47')
cb.outline.set_linewidth(0.8)

# Module legend in gs[3:6, 7]
ax_leg = fig.add_subplot(gs[3:6, 7])
ax_leg.axis('off')
ax_leg.set_title('Signature module', fontsize=9.5, color='#0e2a47', fontweight='bold', loc='left', pad=2)
mod_present = []
for s in sig_cols_sorted:
    cat = SIG_CATEGORY[s]
    if cat not in mod_present: mod_present.append(cat)
y_step = 0.08
for i, mod in enumerate(mod_present):
    y = 0.92 - i*y_step
    ax_leg.add_patch(Rectangle((0.05, y-0.025), 0.18, 0.05, color=CAT_COLORS_DEEP[mod], transform=ax_leg.transAxes))
    ax_leg.text(0.27, y, mod, fontsize=8.5, va='center', transform=ax_leg.transAxes, color='#0e2a47')

# Bottom legend: clinical annotations (row 6, col 1)
ax_clin_leg = fig.add_subplot(gs[6,:])
ax_clin_leg.axis('off')
clin_handles = []
for k,c in PAL.items(): clin_handles.append(mpatches.Patch(color=c, label=f'{k}'))
for k,c in PAL_STAGE_DEEP.items(): clin_handles.append(mpatches.Patch(color=c, label=k))
clin_handles.append(mpatches.Patch(color='#0e2a47', label='M'))
clin_handles.append(mpatches.Patch(color='#a01b2b', label='F'))
clin_handles.append(mpatches.Patch(color=plt.cm.YlOrRd(0.7), label='Age (gradient)'))
clin_handles.append(mpatches.Patch(color=plt.cm.Purples(0.7), label='TMB (gradient)'))
ax_clin_leg.legend(handles=clin_handles, loc='center', fontsize=9, ncol=len(clin_handles), frameon=False)

save_panel(fig, 'Fig3B_signature_heatmap', OUT)

# ============================================================
# FIG 3C — Volcano (no title, arrows close to data inside plot)
# ============================================================
CC_GENES = {'MKI67','TOP2A','CCNB1','CCNB2','CDK1','CDC20','CDC25A','MCM2','MCM3','MCM4','MCM5','MCM6','MCM7',
            'BIRC5','CENPF','PLK1','AURKA','AURKB','BUB1','TYMS','UBE2C','CCNE1','PCNA','RRM2'}
DR_GENES = {'BRCA1','BRCA2','RAD51','RAD51B','RAD51C','RAD51D','PALB2','ATM','ATR','CHEK1','CHEK2','MRE11',
            'RAD50','NBN','XRCC2','XRCC3','FANCA','FANCD2','FANCI','FANCL','BLM','BRIP1','EXO1','POLD1','POLE',
            'MSH2','MSH6','MLH1','PMS2','OGG1','MUTYH'}
EMT_GENES = {'VIM','CDH2','SNAI1','SNAI2','TWIST1','TWIST2','FN1','MMP2','MMP3','MMP9','MMP14','ZEB1','ZEB2',
             'TGFB1','TGFB2','TGFBR1','COL1A1','COL3A1','COL4A1','FAP','ACTA2','S100A4','POSTN','SPARC'}
IM_GENES = {'CD8A','CD8B','GZMA','GZMB','GZMK','GZMH','PRF1','GNLY','IFNG','CD3D','CD3E','HLA-A','HLA-B','HLA-C',
            'B2M','NLRC5','TAP1','TAP2','CIITA','HLA-DRA','HLA-DRB1','MS4A1','CD79A','CD79B','BANK1','FOXP3',
            'CTLA4','PDCD1','LAG3','HAVCR2','TIGIT','CD274','CXCL9','CXCL10','CCL19','CCL21','CXCL13'}

def cat(g):
    if g in CC_GENES: return 'Cell cycle'
    if g in DR_GENES: return 'DNA repair'
    if g in IM_GENES: return 'Immune'
    if g in EMT_GENES: return 'EMT/Stromal'
    return 'Other'

deg_p = deg.copy()
deg_p['cat'] = deg_p.gene.apply(cat)
deg_p['-log10p'] = -np.log10(deg_p.pvalue.replace(0, 1e-300))

VOL_COLORS = {'Cell cycle':'#057a64','DNA repair':'#00567d','Immune':'#c11456','EMT/Stromal':'#b03219','Other':'#cccccc'}

fig, ax = plt.subplots(figsize=(8.5, 6.5))
# Background hexbin
ns = deg_p[deg_p['cat']=='Other']
ax.hexbin(ns.log2FoldChange.clip(-5,5), ns['-log10p'].clip(0,8), gridsize=55, mincnt=2,
          cmap=LinearSegmentedColormap.from_list('h',['#ffffff','#dde2e8']), zorder=1)

for c in ['Cell cycle','DNA repair','Immune','EMT/Stromal']:
    sub = deg_p[deg_p['cat']==c]
    ax.scatter(sub.log2FoldChange.clip(-5,5), sub['-log10p'].clip(0,8),
               s=28, alpha=0.92, color=VOL_COLORS[c],
               edgecolor='#0e2a47', linewidth=0.4, label=c, zorder=3)

others_sig = deg_p[(deg_p['cat']=='Other') & (deg_p.pvalue<1e-3)]
ax.scatter(others_sig.log2FoldChange.clip(-5,5), others_sig['-log10p'].clip(0,8),
           s=14, alpha=0.6, color='#5a6772', edgecolor='none', zorder=2)

# Top labels
to_label = deg_p[(deg_p['cat']!='Other') & (deg_p['-log10p']>2.5)].head(35)
texts = []
for _, r in to_label.iterrows():
    t = ax.text(np.clip(r.log2FoldChange, -5, 5), np.clip(r['-log10p'], 0, 8), r.gene,
                fontsize=8.5, color=VOL_COLORS[r['cat']], fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.3),
                expand_points=(1.3,1.4), force_points=0.4)
except: pass

ax.axvline(0, color='#0e2a47', lw=0.6, alpha=0.5)
ax.axhline(-np.log10(0.01), color='#5a6772', lw=0.7, ls='--', alpha=0.55)
ax.text(4.7, -np.log10(0.01)+0.1, 'p = 0.01', fontsize=8, color='#5a6772', ha='right')

ax.set_xlim(-5.2, 5.2); ax.set_ylim(0, 8.5)

# Direction arrows INSIDE plot near top
y_arrow = 8.0
ax.annotate('', xy=(4.7, y_arrow), xytext=(0.5, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=GOOD_DEEP, lw=2.4, mutation_scale=18))
ax.text(2.6, y_arrow+0.18, 'Up in Good responders', ha='center', va='bottom',
        fontsize=10, color=GOOD_DEEP, fontweight='bold')
ax.annotate('', xy=(-4.7, y_arrow), xytext=(-0.5, y_arrow),
            arrowprops=dict(arrowstyle='-|>', color=BAD_DEEP, lw=2.4, mutation_scale=18))
ax.text(-2.6, y_arrow+0.18, 'Up in Poor responders', ha='center', va='bottom',
        fontsize=10, color=BAD_DEEP, fontweight='bold')

ax.set_xlabel('log2 fold change (good vs poor)', fontsize=10.5, fontweight='bold', color='#0e2a47')
ax.set_ylabel('−log10(p-value)', fontsize=10.5, fontweight='bold', color='#0e2a47')
ax.legend(loc='upper left', fontsize=10, title='Pathway category', title_fontsize=10.5, frameon=False)
add_axis_spines(ax)
save_panel(fig, 'Fig3C_volcano_journal', OUT)

# ============================================================
# FIG 3D — Forest (no title, deeper colors)
# ============================================================
fig, ax = plt.subplots(figsize=(8, 7.5))
pre_st = sig_stats[sig_stats.timepoint=='pre'].sort_values('pvalue').head(20).iloc[::-1].reset_index(drop=True)
np.random.seed(42)
def bootstrap_delta(g_vals, b_vals, n=1000):
    diffs = []
    for _ in range(n):
        gi = np.random.choice(g_vals, len(g_vals), replace=True)
        bi = np.random.choice(b_vals, len(b_vals), replace=True)
        diffs.append(np.mean(gi) - np.mean(bi))
    return np.percentile(diffs, [2.5, 97.5])

cis = []
for _, r in pre_st.iterrows():
    sig_n = r.signature
    g_vals = sigs_m[(sigs_m.timepoint=='pre') & (sigs_m.response_bin=='good')][sig_n].dropna().values
    b_vals = sigs_m[(sigs_m.timepoint=='pre') & (sigs_m.response_bin=='bad')][sig_n].dropna().values
    if len(g_vals)>=3 and len(b_vals)>=3:
        ci = bootstrap_delta(g_vals, b_vals)
    else:
        ci = (0, 0)
    cis.append(ci)
pre_st['ci_low'] = [c[0] for c in cis]
pre_st['ci_high'] = [c[1] for c in cis]

y = np.arange(len(pre_st))
for i, (_, r) in enumerate(pre_st.iterrows()):
    color = GOOD_DEEP if r.delta_good_minus_bad>0 else BAD_DEEP
    ax.plot([r.ci_low, r.ci_high], [i, i], color=color, lw=2.0, alpha=0.85, solid_capstyle='round')
    ax.scatter(r.delta_good_minus_bad, i, s=170 if r.pvalue<0.05 else 100,
               color=color, edgecolor='white', linewidth=1.4, zorder=3)
    star = sig_symbol(r.pvalue)
    if star == 'ns': star = ''
    label = f'p = {r.pvalue:.3g} {star}'
    ax.text(max(r.ci_high, r.delta_good_minus_bad) + 0.07, i, label,
            va='center', ha='left', fontsize=9, color='#0e2a47')

ax.axvline(0, color='#0e2a47', lw=1.0)
ax.axvspan(-0.05, 0.05, color='#dee2e6', alpha=0.4)
ax.set_yticks(y)
ax.set_yticklabels([s.replace('_',' ') for s in pre_st.signature], fontsize=9.5, color='#0e2a47')
ax.set_xlabel('Δ z-score (good − poor),  95% bootstrap CI', fontsize=10.5, fontweight='bold', color='#0e2a47')
add_axis_spines(ax)
ax.set_xlim(min(pre_st.ci_low.min(), -0.5)*1.05, max(pre_st.ci_high.max(), 1.5)*1.6)
save_panel(fig, 'Fig3D_forest_lollipop', OUT)

# ============================================================
# FIG 3E — Bivariate (no title, smaller points)
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

ax_main.set_xlabel('CD8 cytolytic activity (z-score)', fontsize=10.5, fontweight='bold', color='#0e2a47')
ax_main.set_ylabel('CD8 exhaustion (z-score)', fontsize=10.5, fontweight='bold', color='#0e2a47')
ax_main.legend(loc='lower left', fontsize=10, frameon=False)

ax_top.set_xticks([]); ax_top.set_yticks([])
ax_right.set_xticks([]); ax_right.set_yticks([])
for s in ['top','right','bottom']: ax_top.spines[s].set_visible(False)
for s in ['top','right','left']: ax_right.spines[s].set_visible(False)
add_axis_spines(ax_main)
save_panel(fig, 'Fig3E_CD8_biaxial', OUT)

# ============================================================
# FIG 3F — TLS (no title, no overlapping text, tighter spacing)
# ============================================================
TLS_GENES_CABRITA = ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3','CD79B','MS4A1','CCL18','PTGDS','CXCL8']
present_tls = [g for g in TLS_GENES_CABRITA if g in log_tpm.index]
pre = sigs_m[sigs_m.timepoint=='pre'].sort_values(['response_bin','TLS_Cabrita'], ascending=[True, False])
pre_samples = pre.sample_id.tolist()
tls_mat = log_tpm.loc[present_tls, pre_samples]
tls_z = tls_mat.sub(tls_mat.mean(axis=1), axis=0).div(tls_mat.std(axis=1), axis=0)
n_s = len(pre_samples); n_g = len(present_tls)

fig = plt.figure(figsize=(13.5, 6.0))
gs = fig.add_gridspec(4, 2, height_ratios=[0.55, 0.16, 0.16, 4.2], width_ratios=[14, 0.4],
                      hspace=0.06, wspace=0.04)

# Top: TLS bar
ax_top = fig.add_subplot(gs[0, 0])
tls_sig_vals = pre.TLS_Cabrita.values
colors = [PAL[r] for r in pre.response_bin]
ax_top.bar(np.arange(n_s), tls_sig_vals, color=colors, edgecolor='white', linewidth=0.5, width=0.92)
ax_top.axhline(0, color='#0e2a47', lw=0.6)
ax_top.set_ylabel('TLS sig. (z)', fontsize=9.5, fontweight='bold', color='#0e2a47')
ax_top.set_xticks([]); ax_top.set_xlim(-0.5, n_s-0.5)
add_axis_spines(ax_top)
ax_top.tick_params(labelsize=8)

# Annotation strips
def ann_strip2(ax, vals, label, palette_map=None, cmap=None, vmin=None, vmax=None):
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
    ax.text(-0.5, 0.5, label, ha='right', va='center', fontsize=9.5, fontweight='bold', color='#0e2a47')
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, n_s)

ax_a1 = fig.add_subplot(gs[1, 0])
ann_strip2(ax_a1, pre.response_bin.values.tolist(), 'Response', PAL)
ax_a2 = fig.add_subplot(gs[2, 0])
ann_strip2(ax_a2, pre.MHC_II.values.tolist(), 'MHC-II z', cmap=plt.cm.Purples,
           vmin=pre.MHC_II.min(), vmax=pre.MHC_II.max())

# Main heatmap
ax_main = fig.add_subplot(gs[3, 0])
im = ax_main.imshow(tls_z.values, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                     aspect='auto', interpolation='nearest',
                     extent=[0, n_s, 0, n_g])
# Y-tick labels: place INSIDE plot using text, properly spaced (avoid overlap)
ax_main.set_yticks(np.arange(n_g)+0.5)
ax_main.set_yticklabels(present_tls[::-1], fontsize=10.5, fontstyle='italic', color='#0e2a47')
ax_main.invert_yaxis()
# X-tick labels: only every 2nd or use rotation 90 with smaller font
ax_main.set_xticks(np.arange(n_s)+0.5)
ax_main.set_xticklabels(pre_samples, fontsize=7, rotation=90, color='#0e2a47')
ax_main.tick_params(length=0)
ax_main.set_xlim(0, n_s); ax_main.set_ylim(0, n_g)
ax_main.invert_yaxis()
for s in ['top','right','left','bottom']: ax_main.spines[s].set_visible(False)

# Colorbar in right column (spans rows 0-3)
cax = fig.add_subplot(gs[:, 1])
cb = fig.colorbar(im, cax=cax)
cb.set_label('log2(TPM+1) z-score', fontsize=9.5, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.7)

# Stats annotation top-right INSIDE the figure
g = pre[pre.response_bin=='good'].TLS_Cabrita
b = pre[pre.response_bin=='bad'].TLS_Cabrita
p_tls = stats.mannwhitneyu(g, b).pvalue
ax_top.text(0.99, 0.95, f'TLS sig: good vs poor  p = {p_tls:.3f}',
            transform=ax_top.transAxes, ha='right', va='top',
            fontsize=10, color='#0e2a47', fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='#0e2a47', alpha=0.9, boxstyle='round,pad=0.4'))

save_panel(fig, 'Fig3F_TLS_Cabrita', OUT)

print('\n=== Fig 3 v3.1 polished panels saved ===')
