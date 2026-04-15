"""
Fig 4 v3.3 — final audit & polish.
  4A: add LE gene labels (top 3 per pathway), thicken response strip
  4B: tighten size legend labelspacing
  4C: add 'Source' header label above source color bar
  4E: separate colorbar from module legend (proper gridspec cells, no inset_axes)
"""
import sys; sys.path.insert(0, '/mnt/sda1/data/TNT/analysis/scripts')
from _fig_style import *
setup_style()
from pathlib import Path
import matplotlib.gridspec as gridspec
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from adjustText import adjust_text
import gseapy as gp
import warnings; warnings.filterwarnings('ignore')

GOOD_DEEP = '#0a7d6e'; BAD_DEEP = '#c53e1f'; BLACK_DEEP = '#0e2a47'
PAL = {'good':GOOD_DEEP, 'bad':BAD_DEEP}

HALLMARK_CAT = {
    'E2F_TARGETS':'Proliferation','G2M_CHECKPOINT':'Proliferation','MYC_TARGETS_V1':'Proliferation',
    'MYC_TARGETS_V2':'Proliferation','MITOTIC_SPINDLE':'Proliferation','MTORC1_SIGNALING':'Signaling',
    'DNA_REPAIR':'DNA repair','UV_RESPONSE_DN':'DNA repair','UNFOLDED_PROTEIN_RESPONSE':'Stress',
    'EPITHELIAL_MESENCHYMAL_TRANSITION':'Stromal/EMT','MYOGENESIS':'Stromal/EMT','APICAL_JUNCTION':'Stromal/EMT',
    'TGF_BETA_SIGNALING':'Stromal/EMT','HEDGEHOG_SIGNALING':'Signaling','WNT_BETA_CATENIN_SIGNALING':'Signaling',
    'NOTCH_SIGNALING':'Signaling','PI3K_AKT_MTOR_SIGNALING':'Signaling','KRAS_SIGNALING_UP':'Signaling',
    'KRAS_SIGNALING_DN':'Signaling','IL2_STAT5_SIGNALING':'Immune','IL6_JAK_STAT3_SIGNALING':'Immune',
    'INTERFERON_GAMMA_RESPONSE':'Immune','INTERFERON_ALPHA_RESPONSE':'Immune','INFLAMMATORY_RESPONSE':'Immune',
    'TNFA_SIGNALING_VIA_NFKB':'Immune','COMPLEMENT':'Immune','ALLOGRAFT_REJECTION':'Immune','COAGULATION':'Immune',
    'OXIDATIVE_PHOSPHORYLATION':'Metabolism','GLYCOLYSIS':'Metabolism','HYPOXIA':'Metabolism',
    'CHOLESTEROL_HOMEOSTASIS':'Metabolism','FATTY_ACID_METABOLISM':'Metabolism','XENOBIOTIC_METABOLISM':'Metabolism',
    'ADIPOGENESIS':'Metabolism','PEROXISOME':'Metabolism','HEME_METABOLISM':'Metabolism',
    'BILE_ACID_METABOLISM':'Metabolism','ESTROGEN_RESPONSE_EARLY':'Signaling','ESTROGEN_RESPONSE_LATE':'Signaling',
    'ANDROGEN_RESPONSE':'Signaling','P53_PATHWAY':'TP53','APOPTOSIS':'Apoptosis','PROTEIN_SECRETION':'Other',
    'UV_RESPONSE_UP':'Other','ANGIOGENESIS':'Other','PANCREAS_BETA_CELLS':'Other','SPERMATOGENESIS':'Other',
    'REACTIVE_OXYGEN_SPECIES_PATHWAY':'Stress',
}
GSEA_CAT_DEEP = {
    'Proliferation':'#057a64','DNA repair':'#00567d','Immune':'#c11456','Stromal/EMT':'#b03219',
    'Metabolism':'#d4a300','Signaling':'#7a3aad','Apoptosis':'#5a3582','Stress':'#0099b8',
    'TP53':'#ad2831','Other':'#5a6772',
}

ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels_v3'

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')
gsea_h = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
gsea_r = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
ssg = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce')
ssg_stats = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_response_stats.tsv', sep='\t')
tpm = pd.read_csv(ROOT/'06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)

deg_p = deg.dropna(subset=['log2FoldChange','gene']).sort_values('log2FoldChange', ascending=False).reset_index(drop=True)
ranks = pd.Series(deg_p.log2FoldChange.values, index=deg_p.gene.values)
hm_sets = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')

def gseapy_pathway_key(name):
    if name.startswith('HALLMARK_'):
        return name.replace('HALLMARK_','').replace('_',' ').title()
    return name

def running_es(ranks, gene_set):
    n = len(ranks)
    in_set = ranks.index.isin(gene_set).astype(int)
    n_hits = in_set.sum()
    if n_hits == 0: return None, None
    weights = np.abs(ranks.values)
    sum_w_in = (weights * in_set).sum()
    if sum_w_in == 0: return None, None
    increment = (weights * in_set) / sum_w_in
    decrement = (1 - in_set) / (n - n_hits)
    return np.cumsum(increment - decrement), in_set

# ============================================================
# 4A — leading-edge gene labels + thicker response strip
# ============================================================
top_paths = [
    ('HALLMARK_E2F_TARGETS', 'E2F targets'),
    ('HALLMARK_G2M_CHECKPOINT', 'G2M checkpoint'),
    ('HALLMARK_MYC_TARGETS_V1', 'MYC targets V1'),
    ('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT (mesenchymal)'),
]
pre_samples = rna_inv[rna_inv.timepoint=='pre'].sample_id.tolist()
pre_resp = {s: rna_inv[rna_inv.sample_id==s].response_bin.iloc[0] for s in pre_samples}
sample_order_pre = sorted(pre_samples, key=lambda s: (pre_resp[s]!='good', s))

fig = plt.figure(figsize=(15, 11))
outer_gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.20)

for idx, (pw, label) in enumerate(top_paths):
    row, col = divmod(idx, 2)
    inner = outer_gs[row, col].subgridspec(4, 1, height_ratios=[3.0, 0.32, 1.6, 0.20], hspace=0.12)
    ax_es = fig.add_subplot(inner[0])
    ax_rug = fig.add_subplot(inner[1])
    ax_hm = fig.add_subplot(inner[2])
    ax_resp = fig.add_subplot(inner[3])

    key = gseapy_pathway_key(pw)
    gset = set(hm_sets.get(key, []))
    if not gset:
        for k in hm_sets:
            if k.lower().replace(' ','').replace('-','') == pw.replace('HALLMARK_','').replace('_','').lower():
                gset = set(hm_sets[k]); break
    es, in_set = running_es(ranks, gset)
    if es is None:
        for a in [ax_es, ax_rug, ax_hm, ax_resp]: a.set_visible(False)
        continue
    row_g = gsea_h[gsea_h.pathway==pw]
    nes_val = row_g.NES.iloc[0] if len(row_g) else 0
    pval = row_g.pval.iloc[0] if len(row_g) else 1
    color = GOOD_DEEP if nes_val>0 else BAD_DEEP

    x = np.arange(len(es))
    ax_es.plot(x, es, color=color, lw=2.6)
    ax_es.fill_between(x, 0, es, color=color, alpha=0.25)
    ax_es.axhline(0, color='#0e2a47', lw=0.6, alpha=0.5)
    if nes_val>0: peak = np.argmax(es); peak_y = es[peak]
    else: peak = np.argmin(es); peak_y = es[peak]
    ax_es.scatter([peak],[peak_y], color=color, s=110, edgecolor='white', linewidth=1.8, zorder=6)
    ax_es.set_ylabel('ES', fontsize=10, color='#0e2a47', fontweight='bold')
    ax_es.set_xlim(0, len(es)); ax_es.set_xticks([])
    add_axis_spines(ax_es)
    ax_es.text(0.02, 0.96, label, transform=ax_es.transAxes, fontsize=11, fontweight='bold',
               color='#0e2a47', va='top', ha='left',
               bbox=dict(facecolor='white', edgecolor=color, alpha=0.92, boxstyle='round,pad=0.4'))
    ax_es.text(0.98, 0.04 if nes_val>0 else 0.96, f'NES = {nes_val:+.2f}\np = {pval:.2g}',
               transform=ax_es.transAxes, fontsize=10, color=color, fontweight='bold',
               ha='right', va='bottom' if nes_val>0 else 'top')

    # Leading-edge rug
    hit_pos = np.where(in_set==1)[0]
    for hp in hit_pos:
        ax_rug.axvline(hp, color=color, lw=0.65, alpha=0.85)
    ax_rug.set_xlim(0, len(es)); ax_rug.set_yticks([]); ax_rug.set_xticks([])
    ax_rug.text(-0.005, 0.5, 'Leading\nedge', transform=ax_rug.transAxes, ha='right', va='center',
                fontsize=8.5, color='#0e2a47')
    for s in ['top','right','left','bottom']: ax_rug.spines[s].set_visible(False)

    # LE gene heatmap with TOP 3 gene labels
    if nes_val>0:
        leading = [g for i,g in enumerate(ranks.index) if in_set[i]==1 and i<=peak]
    else:
        leading = [g for i,g in enumerate(ranks.index) if in_set[i]==1 and i>=peak]
    leading = [g for g in leading if g in log_tpm.index][:20]
    if leading:
        sub_mat = log_tpm.loc[leading, sample_order_pre]
        sub_z = sub_mat.sub(sub_mat.mean(axis=1), axis=0).div(sub_mat.std(axis=1), axis=0)
        ax_hm.imshow(sub_z.values, cmap='RdBu_r', vmin=-2, vmax=2, aspect='auto',
                     interpolation='nearest', extent=[0, len(sample_order_pre), 0, len(leading)])
        # Annotate top 3 leading-edge genes as text on right side
        top3 = leading[:3]
        for i, g in enumerate(top3):
            row_idx = leading.index(g)
            y_pos = (len(leading) - row_idx) - 0.5
            ax_hm.annotate(g, xy=(len(sample_order_pre), y_pos),
                           xytext=(len(sample_order_pre)+1, y_pos),
                           fontsize=8, color=color, fontweight='bold',
                           va='center', ha='left',
                           arrowprops=dict(arrowstyle='-', color=color, lw=0.7))
        ax_hm.set_yticks([]); ax_hm.set_xticks([])
        ax_hm.set_xlim(0, len(sample_order_pre))
        ax_hm.text(-0.005, 0.5, f'Top leading-edge\ngenes (n={len(leading)})',
                   transform=ax_hm.transAxes, ha='right', va='center', fontsize=8.5, color='#0e2a47')
        for s in ['top','right','left','bottom']: ax_hm.spines[s].set_visible(False)

        # Separate response strip below
        for j, s in enumerate(sample_order_pre):
            ax_resp.add_patch(Rectangle((j, 0), 1, 1, color=PAL[pre_resp[s]], clip_on=False))
        ax_resp.set_xlim(0, len(sample_order_pre)); ax_resp.set_ylim(0, 1)
        ax_resp.set_yticks([]); ax_resp.set_xticks([])
        ax_resp.text(-0.005, 0.5, 'Response', transform=ax_resp.transAxes, ha='right', va='center',
                     fontsize=8.5, color='#0e2a47')
        for s in ['top','right','left','bottom']: ax_resp.spines[s].set_visible(False)
        ax_resp.set_xlabel('Pre-treatment samples (good ← → poor)', fontsize=9, color='#0e2a47',
                           labelpad=2)

save_panel(fig, 'Fig4A_running_ES', OUT)

# ============================================================
# 4B — tighter size legend
# ============================================================
fig, ax = plt.subplots(figsize=(11, 7.5))
h = gsea_h.copy()
h['name'] = h.pathway.str.replace('HALLMARK_','')
h['category'] = h.name.map(HALLMARK_CAT).fillna('Other')
h['neglogp'] = -np.log10(h.pval.replace(0, 1e-300))

cat_handles = []
for cat, color in GSEA_CAT_DEEP.items():
    sub = h[h.category==cat]
    if len(sub)==0: continue
    sizes = sub['size'] * 1.5
    ax.scatter(sub.NES, sub.neglogp, s=sizes, c=color, alpha=0.85,
               edgecolor='#0e2a47', linewidth=0.8, zorder=3)
    cat_handles.append(mpatches.Patch(color=color, label=cat))

texts = []
sig_h = h[h.pval < 1e-3]
for _, r in sig_h.iterrows():
    name_short = r['name'].replace('_',' ').title()[:30]
    t = ax.text(r.NES, r.neglogp, name_short, fontsize=8.5, color='#0e2a47', fontweight='bold')
    texts.append(t)
try:
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='#5a6772', lw=0.4),
                expand_points=(1.5, 1.6), force_points=0.6)
except: pass

ax.axvline(0, color='#0e2a47', lw=0.7, alpha=0.5)
ax.axhline(-np.log10(0.05), color='#5a6772', lw=0.6, ls='--', alpha=0.5)
ax.text(ax.get_xlim()[1]*0.97, -np.log10(0.05)+0.5, 'p = 0.05', fontsize=8.5, color='#5a6772', ha='right')
ax.set_xlabel('Normalized enrichment score (NES)\n← UP in poor          UP in good →',
              fontsize=11, fontweight='bold', color='#0e2a47')
ax.set_ylabel('−log10(p-value)', fontsize=11, fontweight='bold', color='#0e2a47')

leg1 = ax.legend(handles=cat_handles, loc='upper left', fontsize=9, title='Category',
                 title_fontsize=10, ncol=1, bbox_to_anchor=(1.01, 1.0), frameon=False)
ax.add_artist(leg1)

# Size legend — tight spacing
size_examples = [25, 100, 200]
size_handles = [ax.scatter([], [], s=sz*1.5, color='#5a6772', edgecolor='#0e2a47', linewidth=0.8,
                            alpha=0.7, label=f'{sz} genes')
                for sz in size_examples]
ax.legend(handles=size_handles, loc='lower left', fontsize=9, title='Gene set size',
          title_fontsize=10, ncol=1, bbox_to_anchor=(1.01, 0.0), frameon=False,
          scatterpoints=1, labelspacing=1.0, handletextpad=1.0, borderpad=0.6)

add_axis_spines(ax)
save_panel(fig, 'Fig4B_Hallmark_bubble', OUT)

# ============================================================
# 4C — add 'Source' header label
# ============================================================
h_top = gsea_h.copy()
h_top['name'] = h_top.pathway.str.replace('HALLMARK_','')
h_top['source'] = 'Hallmark'
h_top['category'] = h_top.name.map(HALLMARK_CAT).fillna('Other')
r_top = gsea_r.copy()
r_top['name'] = r_top.pathway.str.replace('REACTOME_','').str.replace('_',' ').str.title()
r_top['source'] = 'Reactome'
def reactome_cat(n):
    n = n.lower()
    if 'cycle' in n or 'mitotic' in n or 'm phase' in n or 's phase' in n: return 'Proliferation'
    if 'repair' in n or 'recombination' in n or 'replication' in n or 'damage' in n: return 'DNA repair'
    if 'immune' in n or 'antigen' in n or 'mhc' in n or 'interferon' in n or 'inflammatory' in n: return 'Immune'
    if 'extracellular matrix' in n or 'ecm' in n or 'collagen' in n: return 'Stromal/EMT'
    if 'metabol' in n or 'glycolysis' in n or 'oxidative' in n: return 'Metabolism'
    if 'signal' in n: return 'Signaling'
    if 'apopt' in n: return 'Apoptosis'
    return 'Other'
r_top['category'] = r_top.name.apply(reactome_cat)
combined = pd.concat([h_top.assign(label=h_top.name.str.replace('_',' ').str.title()),
                      r_top.assign(label=r_top.name)], ignore_index=True)
combined = combined.sort_values('pval').head(25).sort_values('NES').reset_index(drop=True)

fig = plt.figure(figsize=(11.5, 8.5))
# Add a small top row for Source header
gs = fig.add_gridspec(2, 4, height_ratios=[0.05, 1.0],
                      width_ratios=[0.20, 4.5, 6.0, 0.6], wspace=0.02, hspace=0.02)

# 'Source' header label at top of col 0
ax_src_hdr = fig.add_subplot(gs[0, 0])
ax_src_hdr.axis('off')
ax_src_hdr.text(0.5, 0.0, 'Source', ha='center', va='bottom', fontsize=9.5, fontweight='bold',
                color='#0e2a47', transform=ax_src_hdr.transAxes)

# Source color bar
ax_src = fig.add_subplot(gs[1, 0])
SRC_COLOR = {'Hallmark':'#264653', 'Reactome':'#e9c46a'}
src_arr = np.array([[to_rgb(SRC_COLOR[s]) for s in combined.source]])
ax_src.imshow(src_arr.transpose(1,0,2), aspect='auto', interpolation='nearest',
              extent=[0,1,0,len(combined)])
ax_src.set_xticks([]); ax_src.set_yticks([])
ax_src.invert_yaxis()
for s in ['top','right','left','bottom']: ax_src.spines[s].set_visible(False)

# Labels
ax_lbl = fig.add_subplot(gs[1, 1])
ax_lbl.axis('off')
ax_lbl.set_xlim(0,1); ax_lbl.set_ylim(0, len(combined))
for i, r in combined.iterrows():
    ax_lbl.text(0.99, len(combined)-1-i+0.5, r.label[:55], ha='right', va='center',
                fontsize=9.5, color='#0e2a47')

# Dot column
ax_dot = fig.add_subplot(gs[1, 2])
y = np.arange(len(combined))
sizes = combined['size'] * 1.5
colors = [GSEA_CAT_DEEP[c] for c in combined.category]
ax_dot.axvline(0, color='#0e2a47', lw=0.9)
for i, r in combined.iterrows():
    color = GSEA_CAT_DEEP[r.category]
    ax_dot.plot([0, r.NES], [len(combined)-1-i, len(combined)-1-i], color=color, lw=1.2, alpha=0.45)
ax_dot.scatter(combined.NES, len(combined)-1-y, s=sizes, c=colors, alpha=0.92,
               edgecolor='#0e2a47', linewidth=0.7, zorder=3)
ax_dot.set_yticks([]); ax_dot.set_ylim(-0.5, len(combined)-0.5)
ax_dot.set_xlabel('Normalized enrichment score (NES)', fontsize=11, fontweight='bold', color='#0e2a47')
add_axis_spines(ax_dot, sides=('bottom',))
ax_dot.spines['left'].set_visible(False)

# Legend col
ax_leg = fig.add_subplot(gs[1, 3])
ax_leg.axis('off')
ax_leg.text(0.05, 0.99, 'Source', transform=ax_leg.transAxes, fontsize=10, fontweight='bold',
            color='#0e2a47', va='top', ha='left')
for i, (k, c) in enumerate(SRC_COLOR.items()):
    y_pos = 0.93 - i*0.05
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.018), 0.25, 0.035, color=c, transform=ax_leg.transAxes))
    ax_leg.text(0.34, y_pos, k, transform=ax_leg.transAxes, fontsize=9, va='center', color='#0e2a47')

ax_leg.text(0.05, 0.77, 'Category', transform=ax_leg.transAxes, fontsize=10, fontweight='bold',
            color='#0e2a47', va='top', ha='left')
cats_used = list(combined.category.unique())
for i, c in enumerate(cats_used):
    y_pos = 0.71 - i*0.055
    ax_leg.add_patch(Rectangle((0.05, y_pos-0.018), 0.25, 0.035, color=GSEA_CAT_DEEP[c],
                                transform=ax_leg.transAxes))
    ax_leg.text(0.34, y_pos, c, transform=ax_leg.transAxes, fontsize=9, va='center', color='#0e2a47')

save_panel(fig, 'Fig4C_pathway_dotplot', OUT)

# ============================================================
# 4E — proper gridspec for colorbar + module legend (no inset_axes)
# ============================================================
ssg_m = ssg.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in ssg_m.columns:
    ssg_m = ssg.reset_index(); ssg_m.columns = ['sample_id'] + list(ssg_m.columns[1:])
ssg_m = ssg_m.merge(rna_inv[['sample_id','timepoint','response_bin']], on='sample_id')
pre = ssg_m[ssg_m.timepoint=='pre']
top30 = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(30).pathway.tolist()
present = [p for p in top30 if p in pre.columns]
mat = pre.set_index('sample_id')[present]
mat_z = mat.sub(mat.mean()).div(mat.std())
order_samples = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
mat_z = mat_z.loc[order_samples]
n_s = len(order_samples); n_p = len(present)

from scipy.cluster.hierarchy import linkage, leaves_list
Z = linkage(mat_z.T.values, method='average', metric='correlation')
order_idx = leaves_list(Z)
present_ordered = [present[i] for i in order_idx]
mat_z = mat_z[present_ordered]

def pw_cat(p):
    pl = p.lower()
    if 'cycle' in pl or 'mitotic' in pl or 'phase' in pl or 'cdk' in pl: return 'Proliferation'
    if 'repair' in pl or 'recombination' in pl or 'damage' in pl or 'hdr' in pl or 'replication' in pl: return 'DNA repair'
    if 'immune' in pl or 'antigen' in pl or 'mhc' in pl or 'interferon' in pl or 'inflammatory' in pl: return 'Immune'
    if 'extracellular' in pl or 'ecm' in pl or 'collagen' in pl or 'epithelial' in pl: return 'Stromal/EMT'
    if 'metabol' in pl or 'glycolysis' in pl or 'oxidative' in pl: return 'Metabolism'
    if 'signal' in pl: return 'Signaling'
    if 'myc' in pl or 'e2f' in pl or 'g2m' in pl: return 'Proliferation'
    return 'Other'

fig = plt.figure(figsize=(15.5, 9))
# Right column split: colorbar (rows 0-1) + module legend (rows 2-onwards)
gs = fig.add_gridspec(5, 6,
    height_ratios=[0.18, 1.5, 0.15, 3.5, 0.7],
    width_ratios=[0.05, 1.6, 0.18, 6.0, 0.4, 1.5],
    hspace=0.06, wspace=0.0)

def ann_strip(ax, vals, palette_map):
    arr = np.array([[to_rgb(palette_map.get(v, '#ecf0f1')) for v in vals]])
    ax.imshow(arr, aspect='auto', interpolation='nearest', extent=[0, n_s, 0, 1])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ['top','right','left','bottom']: ax.spines[s].set_visible(False)
    ax.set_xlim(0, n_s)

# Top annotation
resp_a = pd.Series([pre.set_index('sample_id').loc[s,'response_bin'] for s in order_samples])
ax_ann = fig.add_subplot(gs[0, 3])
ann_strip(ax_ann, resp_a.tolist(), PAL)
ax_ann_lbl = fig.add_subplot(gs[0, 1])
ax_ann_lbl.axis('off')
ax_ann_lbl.text(1.0, 0.5, 'Response', ha='right', va='center', fontsize=10.5,
                color='#0e2a47', fontweight='bold', transform=ax_ann_lbl.transAxes)

# Pathway labels & module strip & heatmap (rows 1, 3 — split for visual breathing)
# Combine into one row 1+2+3 via single big row 3 for heatmap
# Actually let me keep it simpler: heatmap spans rows 1-3 so it's the main element
ax_pw_lbl = fig.add_subplot(gs[1:4, 1])
ax_pw_lbl.axis('off')
ax_pw_lbl.set_xlim(0, 1); ax_pw_lbl.set_ylim(0, n_p)
ax_pw_lbl.invert_yaxis()
for i, p in enumerate(present_ordered):
    short = p.split(' R-HSA')[0][:48]
    ax_pw_lbl.text(1.0, i+0.5, short, ha='right', va='center', fontsize=8.5, color='#0e2a47')

ax_mod = fig.add_subplot(gs[1:4, 2])
mod_colors = [GSEA_CAT_DEEP[pw_cat(p)] for p in present_ordered]
arr_mod = np.array([[to_rgb(c) for c in mod_colors]])
ax_mod.imshow(arr_mod.transpose(1,0,2), aspect='auto', interpolation='nearest', extent=[0,1,0,n_p])
ax_mod.set_xticks([]); ax_mod.set_yticks([])
for s in ['top','right','left','bottom']: ax_mod.spines[s].set_visible(False)
ax_mod.invert_yaxis()

ax_h = fig.add_subplot(gs[1:4, 3])
im = ax_h.imshow(mat_z.values.T, cmap='RdBu_r', vmin=-2.5, vmax=2.5,
                 aspect='auto', interpolation='nearest', extent=[0, n_s, 0, n_p])
ax_h.set_xticks([]); ax_h.set_yticks([])
ax_h.invert_yaxis()
for s in ['top','right','left','bottom']: ax_h.spines[s].set_visible(False)

# Right column — colorbar at row 1 (top), module legend at row 3 (below) — separated by gs[2,5] gap
cax = fig.add_subplot(gs[1, 5])
cb = fig.colorbar(im, cax=cax, orientation='vertical')
cb.set_label('z-score', fontsize=10, color='#0e2a47', fontweight='bold')
cb.ax.tick_params(labelsize=8.5)
cb.outline.set_edgecolor('#0e2a47'); cb.outline.set_linewidth(0.8)

# Module legend at row 3 col 5
ax_modleg = fig.add_subplot(gs[3, 5])
ax_modleg.axis('off')
ax_modleg.text(0.0, 1.0, 'Pathway category', fontsize=10, color='#0e2a47',
               fontweight='bold', ha='left', va='top', transform=ax_modleg.transAxes)
mod_present = []
for p in present_ordered:
    c = pw_cat(p)
    if c not in mod_present: mod_present.append(c)
y_step = 0.10
for i, mod in enumerate(mod_present):
    y = 0.92 - i*y_step
    ax_modleg.add_patch(Rectangle((0.05, y-0.025), 0.18, 0.05, color=GSEA_CAT_DEEP[mod],
                                   transform=ax_modleg.transAxes))
    ax_modleg.text(0.27, y, mod, fontsize=9.5, va='center',
                   transform=ax_modleg.transAxes, color='#0e2a47')

# Bottom legend
ax_clin = fig.add_subplot(gs[4, 1:4])
ax_clin.axis('off')
ax_clin.legend(handles=[mpatches.Patch(color=GOOD_DEEP, label='Good'),
                        mpatches.Patch(color=BAD_DEEP, label='Poor')],
               loc='center', fontsize=10, ncol=2, frameon=False)

save_panel(fig, 'Fig4E_ssGSEA_heatmap', OUT)

print('=== Fig 4 v3.3 audit polishes saved (4A, 4B, 4C, 4E) ===')
