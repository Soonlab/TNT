"""
Publication-grade figure suite for TNT manuscript.
Each panel saved as individual PDF+PNG with consistent styling.
Figure naming: Fig{N}{letter} for letter-subpanels (e.g., Fig1A.pdf).

Figure list:
  Fig1: Cohort & clinical (A-E)
  Fig2: WES landscape (A-G)
  Fig3: RNA-seq signatures & DEG (A-F)
  Fig4: GSEA (A-E)
  Fig5: Integration & ML (A-E)
  Fig6: Treatment-induced delta (A-H)
  Fig7: External validation (A-C)
  Fig8: HLA & neoantigen (A-F)
  Fig9: Clonal evolution (A-C)
"""
import pandas as pd, numpy as np, os, re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy import stats

# ---------------------------------------------------------
# Style — nature-like professional
# ---------------------------------------------------------
plt.rcParams.update({
    'font.family':'DejaVu Sans',
    'font.size':9,
    'axes.labelsize':10,
    'axes.titlesize':11,
    'axes.titleweight':'bold',
    'axes.spines.top':False,
    'axes.spines.right':False,
    'axes.linewidth':0.8,
    'xtick.major.size':3,
    'ytick.major.size':3,
    'xtick.major.width':0.8,
    'ytick.major.width':0.8,
    'xtick.direction':'out',
    'ytick.direction':'out',
    'legend.frameon':False,
    'pdf.fonttype':42,
    'ps.fonttype':42,
    'savefig.bbox':'tight',
    'savefig.dpi':300,
})

GOOD = '#2a9d8f'
BAD = '#e76f51'
NORMAL = '#8d99ae'
PALETTE = {'good':GOOD,'bad':BAD,'normal':NORMAL}
PAL_TP = {'pre':'#264653','post':'#e9c46a','normal':'#a8dadc'}

def save_panel(fig, name, out_dir):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir/f'{name}.pdf', bbox_inches='tight')
    fig.savefig(out_dir/f'{name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'  ✓ {name}')

def sig_star(p):
    if pd.isna(p): return 'ns'
    if p<0.001: return '***'
    if p<0.01: return '**'
    if p<0.05: return '*'
    if p<0.1: return '†'
    return 'ns'

def boxstrip(data, x, y, ax, order=None, palette=None, size=4, alpha=0.8,
             title=None, ylabel=None, test=True, expected_direction=None):
    """Beautiful boxplot + jitter, auto p-value label."""
    if order is None and x=='response_bin': order = ['good','bad']
    if palette is None: palette = PALETTE
    sns.boxplot(data=data, x=x, y=y, ax=ax, order=order, hue=x, palette=palette,
                width=0.55, linewidth=1.0, fliersize=0, legend=False, saturation=0.9)
    sns.stripplot(data=data, x=x, y=y, ax=ax, order=order, color='black',
                  size=size, alpha=alpha, jitter=0.15)
    if ylabel: ax.set_ylabel(ylabel)
    if title: ax.set_title(title)
    ax.set_xlabel('')
    if test and order and len(order)==2:
        v1 = pd.to_numeric(data[data[x]==order[0]][y], errors='coerce').dropna()
        v2 = pd.to_numeric(data[data[x]==order[1]][y], errors='coerce').dropna()
        if len(v1)>=3 and len(v2)>=3:
            u = stats.mannwhitneyu(v1, v2)
            p = u.pvalue
            ymax = max(v1.max(), v2.max())
            ymin = min(v1.min(), v2.min())
            span = ymax - ymin
            ax.plot([0,0,1,1], [ymax+span*0.04,ymax+span*0.07,ymax+span*0.07,ymax+span*0.04], color='black', lw=0.9)
            ax.text(0.5, ymax+span*0.08, f'p = {p:.3g}', ha='center', fontsize=8)
            ax.set_ylim(top=ymax+span*0.20)

# ---------------------------------------------------------
# Load all data
# ---------------------------------------------------------
ROOT = Path('/mnt/sda1/data/TNT/analysis')
OUT = ROOT/'figures/panels'; OUT.mkdir(parents=True, exist_ok=True)

clin = pd.read_csv(ROOT/'00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(ROOT/'00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(ROOT/'00_cohort/rna_inventory.tsv', sep='\t')

tmb = pd.read_csv(ROOT/'02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(ROOT/'02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
sbs_full = pd.read_csv(ROOT/'01_wes_signatures/sbs_activities_with_meta.tsv', sep='\t')
sbs_sum = pd.read_csv(ROOT/'01_wes_signatures/sbs_summary_key.tsv', sep='\t')

drv = pd.read_csv(ROOT/'04_wes_cnv_clonal/driver_oncoprint_matrix.tsv', sep='\t', index_col=0)
cnv = pd.read_csv(ROOT/'04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
hrd = pd.read_csv(ROOT/'04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv', sep='\t')

sigs = pd.read_csv(ROOT/'06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
sig_stats = pd.read_csv(ROOT/'06_rna_immune/sig_response_stats.tsv', sep='\t')
ssg = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0).apply(pd.to_numeric, errors='coerce')
ssg_stats = pd.read_csv(ROOT/'08_rna_pathway/ssgsea_response_stats.tsv', sep='\t')
deg = pd.read_csv(ROOT/'05_rna_deg_gsea/DEG_good_vs_bad_pre.tsv', sep='\t')
gsea_h = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')
gsea_r = pd.read_csv(ROOT/'05_rna_deg_gsea/GSEA_Reactome_pre.tsv', sep='\t')

cms = pd.read_csv(ROOT/'07_rna_cms/cms_assignments.tsv', sep='\t')
trust = pd.read_csv(ROOT/'06_rna_immune/trust4_summary.tsv', sep='\t')

hla = pd.read_csv(ROOT/'03_hla/hla_class_I_typing.tsv', sep='\t')
loh = pd.read_csv(ROOT/'03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')
neo = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
neo_delta = pd.read_csv(ROOT/'03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')

integ = pd.read_csv(ROOT/'tables/integrated_subject_master.tsv', sep='\t')
rfeat = pd.read_csv(ROOT/'tables/response_feature_stats.tsv', sep='\t')

ml = pd.read_csv(ROOT/'10_ml_predictor/ml_loocv_results.tsv', sep='\t')
rf_imp = pd.read_csv(ROOT/'10_ml_predictor/rf_feature_importance.tsv', sep='\t')

delta_22 = pd.read_csv(ROOT/'09_integration/paired_delta/delta_22sigs_response.tsv', sep='\t')
delta_ssg = pd.read_csv(ROOT/'09_integration/paired_delta/delta_ssgsea_response.tsv', sep='\t')
delta_trust = pd.read_csv(ROOT/'09_integration/paired_delta/delta_trust4_response.tsv', sep='\t')
delta_sbs = pd.read_csv(ROOT/'09_integration/paired_delta/delta_sbs_response.tsv', sep='\t')

ext_meta = pd.read_csv(ROOT/'11_external_validation/meta_analysis_manual.tsv', sep='\t')
ext_stats = pd.read_csv(ROOT/'11_external_validation/signature_stats_manual.tsv', sep='\t')

pyclone = pd.read_csv(ROOT/'04_wes_cnv_clonal/pyclone/clonal_summary.tsv', sep='\t')

print('Data loaded — generating panels...')
UNMATCHED = [13,15,16,17,18,19,33]

# =========================================================
# FIGURE 1 — Cohort & Clinical
# =========================================================
# 1A — Response pie
fig, ax = plt.subplots(figsize=(3.5,3.5))
rcounts = clin.response_bin.value_counts()
colors = [GOOD if k=='good' else BAD for k in rcounts.index]
wedges, texts, autotexts = ax.pie(rcounts, labels=[f'{k}\n(n={v})' for k,v in rcounts.items()],
                                    colors=colors, autopct='%1.0f%%', startangle=90,
                                    wedgeprops=dict(linewidth=1.5, edgecolor='white'),
                                    textprops={'fontsize':10})
for t in autotexts: t.set_color('white'); t.set_fontweight('bold')
ax.set_title('Response distribution (n=35)')
save_panel(fig, 'Fig1A_response_pie', OUT)

# 1B — cT stage × response stacked bar
fig, ax = plt.subplots(figsize=(4,3.5))
ct_tab = pd.crosstab(clin.cT, clin.response_bin, normalize='columns')*100
stages = ['T2','T2/T3','T3','T4']
ct_tab = ct_tab.reindex(stages).fillna(0)[['good','bad']]
bottom = np.zeros(2)
stage_colors = ['#264653','#2a9d8f','#e9c46a','#e76f51']
for i,t in enumerate(stages):
    ax.bar(range(2), ct_tab.loc[t].values, bottom=bottom, label=t, color=stage_colors[i],
           edgecolor='white', linewidth=1.5)
    bottom += ct_tab.loc[t].values
ax.set_xticks(range(2)); ax.set_xticklabels(['good','bad'])
ax.set_ylabel('Percentage (%)'); ax.set_ylim(0,105)
ax.legend(title='cT stage', loc='center left', bbox_to_anchor=(1.02,0.5), fontsize=9)
p_ct = stats.chi2_contingency(pd.crosstab(clin.response_bin, clin.cT)).pvalue
ax.set_title(f'Clinical T-stage by response\nχ² p = {p_ct:.3f}')
save_panel(fig, 'Fig1B_cT_stacked', OUT)

# 1C — Age × sex
fig, ax = plt.subplots(figsize=(4,3.5))
for resp in ['good','bad']:
    for sex, marker in [('M','o'),('F','^')]:
        sub = clin[(clin.response_bin==resp) & (clin.sex==sex)]
        ax.scatter([0 if resp=='good' else 1]*len(sub) + np.random.uniform(-0.12,0.12,len(sub)),
                   sub.age, marker=marker, s=70, alpha=0.8,
                   color=PALETTE[resp], edgecolor='black', linewidth=0.6, label=f'{resp}-{sex}')
ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'])
ax.set_ylabel('Age (years)')
p_age = stats.mannwhitneyu(clin[clin.response_bin=='good'].age, clin[clin.response_bin=='bad'].age).pvalue
ax.set_title(f'Age distribution by response\np = {p_age:.3f}')
ax.legend(fontsize=8, ncol=2, loc='lower right')
save_panel(fig, 'Fig1C_age_sex', OUT)

# 1D — Sample matrix
fig, ax = plt.subplots(figsize=(13, 3))
subs = sorted(clin.subject_id)
mat = np.zeros((6, len(subs)))
for j, s in enumerate(subs):
    w = wes_inv[wes_inv.subject_id==s]
    r = rna_inv[rna_inv.subject_id==s]
    mat[0,j] = 1 if (w.timepoint=='normal').any() else 0
    mat[1,j] = 1 if (w.timepoint=='pre').any() else 0
    mat[2,j] = 1 if (w.timepoint=='post').any() else 0
    mat[3,j] = 1 if (r.timepoint=='normal').any() else 0
    mat[4,j] = 1 if (r.timepoint=='pre').any() else 0
    mat[5,j] = 1 if (r.timepoint=='post').any() else 0
sns.heatmap(mat, cmap=sns.color_palette(['white','#264653'],as_cmap=False), cbar=False, ax=ax,
            yticklabels=['WES normal','WES pre','WES post','RNA normal','RNA pre','RNA post'],
            xticklabels=subs, linewidths=0.6, linecolor='#cccccc')
# response color row
ax.set_xlabel('Subject ID')
for j, s in enumerate(subs):
    resp = clin[clin.subject_id==s].response_bin.iloc[0]
    ax.add_patch(Rectangle((j, 6.15), 1, 0.15, color=PALETTE[resp], clip_on=False, transform=ax.transData))
ax.text(-0.5, 6.22, 'Response', ha='right', va='center', fontsize=9, transform=ax.transData)
ax.set_title('Sample availability matrix (35 subjects × 6 sample types)')
ax.set_xticklabels(subs, rotation=0, fontsize=7)
save_panel(fig, 'Fig1D_sample_matrix', OUT)

# 1E — Consort-style summary
fig, ax = plt.subplots(figsize=(5,5))
ax.axis('off')
boxes = [
    (0.5, 0.95, 'LARC patients enrolled\nN = 35', '#264653'),
    (0.25, 0.72, 'Good responders\nTRG 0-1\nn = 18', GOOD),
    (0.75, 0.72, 'Poor responders\nTRG 2-3\nn = 17', BAD),
    (0.18, 0.48, 'WES 77 samples\nRNA 31 samples', '#264653'),
    (0.82, 0.48, 'WES 77 samples\nRNA 25 samples', '#264653'),
    (0.5, 0.22, 'Matched pre+post (14)\nUnmatched (21)', '#8d99ae'),
    (0.5, 0.05, 'Analyses:\nWES somatic · SBS · MSI · CNV · HLA\nRNA DEG · GSEA · ssGSEA · TRUST4\nIntegration · ML · External validation', '#1d3557'),
]
for x,y,text,color in boxes:
    ax.add_patch(plt.Rectangle((x-0.22,y-0.055), 0.44, 0.11, facecolor=color, alpha=0.85, edgecolor='black', linewidth=1))
    ax.text(x, y, text, ha='center', va='center', fontsize=9, color='white', fontweight='bold')
# arrows
for (x1,y1,x2,y2) in [(0.5,0.89,0.3,0.78),(0.5,0.89,0.7,0.78),(0.27,0.66,0.21,0.54),(0.73,0.66,0.79,0.54),(0.5,0.42,0.5,0.28),(0.5,0.16,0.5,0.12)]:
    ax.annotate('', xy=(x2,y2), xytext=(x1,y1), arrowprops=dict(arrowstyle='->', lw=1.2, color='#264653'))
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.set_title('Study design (CONSORT-style)')
save_panel(fig, 'Fig1E_study_design', OUT)

print('=== Fig 1 done ===')

# =========================================================
# FIGURE 2 — WES Landscape
# =========================================================
tmb_m = tmb[tmb.subject_id.isin(set(clin.subject_id)-set(UNMATCHED))].copy()

# 2A — TMB pre
fig, ax = plt.subplots(figsize=(3.5,3.5))
pre_m = tmb_m[tmb_m.timepoint=='pre']
boxstrip(pre_m, 'response_bin', 'TMB_nonsyn_per_Mb', ax, title='TMB pre-treatment (matched)',
         ylabel='Nonsynonymous TMB (/Mb)')
ax.axhline(10, color='red', ls='--', lw=0.8, alpha=0.6, label='TMB-high (10/Mb)')
ax.legend(fontsize=7, loc='upper right')
save_panel(fig, 'Fig2A_TMB_pre', OUT)

# 2B — MSI
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(msi, 'response_bin', 'MSI_pct', ax, title='MSI percentage (all matched)',
         ylabel='MSI % (microsatellite unstable sites)')
ax.axhline(20, color='red', ls='--', lw=0.8, alpha=0.6, label='MSI-H (>20%)')
ax.legend(fontsize=7, loc='upper right')
save_panel(fig, 'Fig2B_MSI', OUT)

# 2C — Driver oncoprint (pre samples, top 15 genes)
pre_samples = wes_inv[wes_inv.timepoint=='pre'].sample_id.tolist()
wide = drv.drop(columns=['total_samples_mutated'], errors='ignore')
wide = wide[[c for c in wide.columns if c in pre_samples]]
gene_counts = (wide>0).sum(axis=1).sort_values(ascending=False)
top_genes = wide.loc[gene_counts.head(15).index]
# order samples by response then gene count
sample_order = sorted(top_genes.columns, key=lambda x: (
    wes_inv[wes_inv.sample_id==x].response_bin.iloc[0]!='good',
    -(top_genes[x]>0).sum()
))
top_genes = top_genes[sample_order]
binary = (top_genes>0).astype(int)

fig, ax = plt.subplots(figsize=(11, 4.5))
sns.heatmap(binary, cmap=['#f0f0f0','#1d3557'], cbar=False, linewidths=0.4, linecolor='white',
            ax=ax, yticklabels=True, xticklabels=False)
# response bar above
for j, s in enumerate(binary.columns):
    resp = wes_inv[wes_inv.sample_id==s].response_bin.iloc[0]
    matched = wes_inv[wes_inv.sample_id==s].subject_id.iloc[0] not in UNMATCHED
    ax.add_patch(Rectangle((j,-0.8), 1, 0.35, color=PALETTE[resp], clip_on=False))
    if not matched:
        ax.add_patch(Rectangle((j,-0.4), 1, 0.35, color='#e9c46a', clip_on=False))
ax.text(-1.5, -0.6, 'Response', ha='right', va='center', fontsize=8)
ax.text(-1.5, -0.2, 'Unmatched', ha='right', va='center', fontsize=8)
ax.set_xlim(0, binary.shape[1]); ax.set_ylim(-1, binary.shape[0])
# frequency on right
for i, gene in enumerate(binary.index):
    freq = binary.loc[gene].sum()/binary.shape[1]*100
    ax.text(binary.shape[1]+0.3, i+0.5, f'{freq:.0f}%', ha='left', va='center', fontsize=8)
ax.set_ylabel('')
ax.set_title('Driver mutation oncoprint (top 15 genes, pre-treatment samples)')
save_panel(fig, 'Fig2C_driver_oncoprint', OUT)

# 2D — SBS signature stacked bar (top contribs)
sig_cols = [c for c in sbs_full.columns if c.startswith('SBS')]
active = [c for c in sig_cols if (sbs_full[c]>0).any()]
# Only top contributors
total_contrib = sbs_full[active].sum(axis=0).sort_values(ascending=False)
keep_sigs = total_contrib.head(10).index.tolist()
sbs_plot = sbs_full[['sample_id','response_bin','timepoint'] + keep_sigs].copy()
sbs_plot = sbs_plot[sbs_plot.timepoint=='pre']
# sort by response
sbs_plot = sbs_plot.sort_values(['response_bin','sample_id'])
sbs_norm = sbs_plot[keep_sigs].div(sbs_plot[keep_sigs].sum(axis=1).replace(0,np.nan), axis=0).fillna(0)

fig, ax = plt.subplots(figsize=(11,4))
cmap = plt.cm.tab10(range(len(keep_sigs)))
bottom = np.zeros(len(sbs_plot))
for i, sig in enumerate(keep_sigs):
    ax.bar(range(len(sbs_plot)), sbs_norm[sig].values, bottom=bottom, label=sig, color=cmap[i], width=0.95, edgecolor='white', linewidth=0.3)
    bottom += sbs_norm[sig].values
# response color bar under x
for j, (_,r) in enumerate(sbs_plot.reset_index(drop=True).iterrows()):
    ax.add_patch(Rectangle((j-0.5, -0.08), 1, 0.05, color=PALETTE[r.response_bin], clip_on=False))
ax.set_xticks(range(len(sbs_plot)))
ax.set_xticklabels(sbs_plot.sample_id.tolist(), rotation=90, fontsize=6)
ax.set_ylabel('Proportion')
ax.set_ylim(0,1.02)
ax.set_xlim(-0.5, len(sbs_plot)-0.5)
ax.legend(title='Signature', loc='center left', bbox_to_anchor=(1.02,0.5), fontsize=8)
ax.set_title('SBS mutational signature composition (pre-treatment)')
save_panel(fig, 'Fig2D_SBS_stacked', OUT)

# 2E — MMR claim vs Mutect2 refit scatter
fig, ax = plt.subplots(figsize=(4,3.5))
# Map subjects: user claim = 1, 9, 12 MMR; 5, 14 SBS3
claim_mmr = {1,9,12}; claim_sbs3 = {5,14}
pre_sbs = sbs_sum[sbs_sum.timepoint=='pre'].copy()
pre_sbs['claimed_MMR'] = pre_sbs.subject_id.isin(claim_mmr)
pre_sbs['claimed_SBS3'] = pre_sbs.subject_id.isin(claim_sbs3)
# Plot MMR_prop for each pre sample, highlight claimed
all_pre = pre_sbs[~pre_sbs.subject_id.isin(UNMATCHED)]
ax.scatter(all_pre[~all_pre.claimed_MMR].subject_id, all_pre[~all_pre.claimed_MMR].MMR_prop*100,
           color='#8d99ae', s=50, alpha=0.7, label='Other subjects')
ax.scatter(all_pre[all_pre.claimed_MMR].subject_id, all_pre[all_pre.claimed_MMR].MMR_prop*100,
           color=BAD, s=120, alpha=0.9, edgecolor='black', linewidth=1.5, label='Claimed MMR (prior)', marker='*')
for _, r in all_pre.iterrows():
    if r.claimed_MMR:
        ax.annotate(f's{r.subject_id}', (r.subject_id, r.MMR_prop*100), xytext=(3,3), textcoords='offset points', fontsize=8)
ax.set_xlabel('Subject ID')
ax.set_ylabel('MMR signature proportion (%)')
ax.set_title('MMR signature reassessed\n(prior vs Mutect2 refit)')
ax.legend(fontsize=8)
save_panel(fig, 'Fig2E_MMR_reassessed', OUT)

# 2F — CIN by response
fig, ax = plt.subplots(figsize=(3.5,3.5))
cnv_m = cnv[(cnv.timepoint=='pre') & cnv.matched]
boxstrip(cnv_m, 'response_bin', 'CIN', ax, title='Chromosomal instability (pre)',
         ylabel='Fraction genome altered (CIN)')
save_panel(fig, 'Fig2F_CIN', OUT)

# 2G — HRD proxy LST
fig, ax = plt.subplots(figsize=(3.5,3.5))
hrd_m = hrd[(hrd.timepoint=='pre') & hrd.matched]
boxstrip(hrd_m, 'response_bin', 'LST', ax, title='Large-scale transitions (LST)',
         ylabel='LST count (≥10 Mb segments)')
save_panel(fig, 'Fig2G_HRD_LST', OUT)

print('=== Fig 2 done ===')

# =========================================================
# FIGURE 3 — RNA-seq: signatures + DEG
# =========================================================
# Merge signature scores with metadata
sigs_m = sigs.reset_index().rename(columns={'index':'sample_id'})
if 'sample_id' not in sigs_m.columns:
    sigs_m = sigs.reset_index()
    sigs_m.columns = ['sample_id']+list(sigs_m.columns[1:])
sigs_m = sigs_m.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')

# 3A — Key pre signatures boxplot grid
fig, axes = plt.subplots(2, 4, figsize=(14,7))
key_sigs_list = ['CD8_proliferation','CD8_activation','Antigen_presentation','NLRC5_HLA_IFNG',
                 'MHC_II','TLS_Cabrita','TGFb_Mariathasan','EMT_Mak']
for ax, sig in zip(axes.flat, key_sigs_list):
    pre = sigs_m[sigs_m.timepoint=='pre']
    boxstrip(pre, 'response_bin', sig, ax, title=sig.replace('_','\n'))
save_panel(fig, 'Fig3A_signatures_grid', OUT)

# 3B — Pre signature heatmap
fig, ax = plt.subplots(figsize=(10,7))
pre = sigs_m[sigs_m.timepoint=='pre']
sig_cols = [c for c in sigs.columns if c in sigs_m.columns]
pre_mat = pre.set_index('sample_id')[sig_cols]
pre_mat = pre_mat.sub(pre_mat.mean()).div(pre_mat.std())
order = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
pre_mat = pre_mat.loc[order]
resp_colors_list = [PALETTE[pre.set_index('sample_id').loc[s,'response_bin']] for s in order]
hm = sns.clustermap(pre_mat.T, row_cluster=True, col_cluster=False, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
               col_colors=resp_colors_list, figsize=(12,7), xticklabels=False,
               cbar_kws={'label':'z-score'}, linewidths=0.2, linecolor='white')
hm.ax_heatmap.set_xlabel('Samples (ordered by response)')
hm.fig.suptitle('Pre-treatment signature scores — heatmap', y=1.02, fontsize=12, fontweight='bold')
hm.savefig(OUT/'Fig3B_signature_heatmap.pdf', bbox_inches='tight')
hm.savefig(OUT/'Fig3B_signature_heatmap.png', dpi=300, bbox_inches='tight')
plt.close('all')
print('  ✓ Fig3B_signature_heatmap')

# 3C — DEG volcano plot
fig, ax = plt.subplots(figsize=(5,5))
deg_plot = deg.copy()
deg_plot['-log10p'] = -np.log10(deg_plot['pvalue'].replace(0, 1e-300))
deg_plot['significant'] = (deg_plot['pvalue']<0.01) & (deg_plot['log2FoldChange'].abs()>1)
ax.scatter(deg_plot[~deg_plot.significant]['log2FoldChange'], deg_plot[~deg_plot.significant]['-log10p'],
           s=4, alpha=0.4, color='#cccccc')
up = deg_plot[deg_plot.significant & (deg_plot.log2FoldChange>0)]
dn = deg_plot[deg_plot.significant & (deg_plot.log2FoldChange<0)]
ax.scatter(up.log2FoldChange, up['-log10p'], s=10, alpha=0.7, color=GOOD, label=f'Up in good (n={len(up)})')
ax.scatter(dn.log2FoldChange, dn['-log10p'], s=10, alpha=0.7, color=BAD, label=f'Up in bad (n={len(dn)})')
# Label top 10
for _, r in deg_plot.head(10).iterrows():
    ax.annotate(r['gene'], (r.log2FoldChange, r['-log10p']), fontsize=7, alpha=0.8,
                xytext=(3,3), textcoords='offset points')
ax.axvline(1, color='gray', ls='--', lw=0.5); ax.axvline(-1, color='gray', ls='--', lw=0.5)
ax.axhline(-np.log10(0.01), color='gray', ls='--', lw=0.5)
ax.set_xlabel('log2(fold change) good vs bad')
ax.set_ylabel('-log10(p-value)')
ax.set_title('Pre-treatment DEG (DESeq2, covariate-adjusted)')
ax.legend(fontsize=8, loc='upper right')
save_panel(fig, 'Fig3C_DEG_volcano', OUT)

# 3D — CD8 proliferation (main finding)
fig, ax = plt.subplots(figsize=(3.5,3.5))
pre = sigs_m[sigs_m.timepoint=='pre']
boxstrip(pre, 'response_bin', 'CD8_proliferation', ax, title='CD8 proliferation signature',
         ylabel='z-score')
save_panel(fig, 'Fig3D_CD8_proliferation', OUT)

# 3E — Top signature ranked bar (response pvals)
fig, ax = plt.subplots(figsize=(6,5))
pre_stats = sig_stats[sig_stats.timepoint=='pre'].sort_values('pvalue').head(15).iloc[::-1]
colors = [GOOD if d>0 else BAD for d in pre_stats.delta_good_minus_bad]
ax.barh(range(len(pre_stats)), -np.log10(pre_stats.pvalue), color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(pre_stats)))
ax.set_yticklabels(pre_stats.signature.str.replace('_',' '))
ax.axvline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('-log10(p-value) good vs bad (pre)')
ax.set_title('Pre-treatment signature ranking (green=↑good, orange=↑bad)')
save_panel(fig, 'Fig3E_signature_rank', OUT)

# 3F — Post-treatment MHC-II + Treg boxplots (key paired finding)
fig, axes = plt.subplots(1, 3, figsize=(10,3.5))
post = sigs_m[sigs_m.timepoint=='post']
for ax, sig in zip(axes, ['MHC_II','Treg','IFNg_Ayers_18']):
    boxstrip(post, 'response_bin', sig, ax, title=f'{sig} (post)')
save_panel(fig, 'Fig3F_post_immune', OUT)

print('=== Fig 3 done ===')

# =========================================================
# FIGURE 4 — GSEA
# =========================================================
# 4A — Hallmark GSEA barplot
fig, ax = plt.subplots(figsize=(7,7))
h = gsea_h.copy()
h['signedLogP'] = -np.log10(h.pval.replace(0,1e-300)) * np.sign(h.NES)
top = h.nsmallest(20, 'pval').sort_values('signedLogP').reset_index(drop=True)
colors = [GOOD if n>0 else BAD for n in top.NES]
ax.barh(range(len(top)), top.signedLogP, color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([t.replace('HALLMARK_','') for t in top.pathway], fontsize=8)
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('Signed -log10(p-value)   [+ ↑good, − ↑bad]')
ax.set_title('Hallmark GSEA (pre-treatment, good vs bad)')
save_panel(fig, 'Fig4A_Hallmark_GSEA', OUT)

# 4B — Reactome GSEA top
fig, ax = plt.subplots(figsize=(8,6))
r = gsea_r.copy()
r['signedLogP'] = -np.log10(r.pval.replace(0,1e-300)) * np.sign(r.NES)
top = r.nsmallest(15, 'pval').sort_values('signedLogP').reset_index(drop=True)
colors = [GOOD if n>0 else BAD for n in top.NES]
ax.barh(range(len(top)), top.signedLogP, color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([t.replace('REACTOME_','')[:60] for t in top.pathway], fontsize=7)
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('Signed -log10(p-value)')
ax.set_title('Reactome GSEA (pre-treatment)')
save_panel(fig, 'Fig4B_Reactome_GSEA', OUT)

# 4C — ssGSEA top pathway boxplots grid
fig, axes = plt.subplots(2, 3, figsize=(12,7))
ssg_m = ssg.reset_index().rename(columns={'index':'sample_id'}).merge(rna_inv[['sample_id','timepoint','response_bin']], on='sample_id')
top_pre_ssg = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(6).pathway.tolist()
for ax, pw in zip(axes.flat, top_pre_ssg):
    if pw not in ssg_m.columns: ax.set_visible(False); continue
    pre = ssg_m[ssg_m.timepoint=='pre']
    boxstrip(pre, 'response_bin', pw, ax, title=pw.split(' R-HSA')[0][:40])
save_panel(fig, 'Fig4C_ssGSEA_top_boxes', OUT)

# 4D — ssGSEA heatmap (top 30 pre pathways)
fig, ax = plt.subplots(figsize=(10,8))
pre = ssg_m[ssg_m.timepoint=='pre']
top30 = ssg_stats[ssg_stats.timepoint=='pre'].sort_values('pvalue').head(30).pathway.tolist()
present = [p for p in top30 if p in pre.columns]
mat = pre.set_index('sample_id')[present]
mat = mat.sub(mat.mean()).div(mat.std())
order = pre.sort_values(['response_bin','sample_id']).sample_id.tolist()
mat = mat.loc[order]
resp_cols = [PALETTE[pre.set_index('sample_id').loc[s,'response_bin']] for s in order]
hm = sns.clustermap(mat.T, row_cluster=True, col_cluster=False, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
               col_colors=resp_cols, figsize=(12,8), xticklabels=False,
               yticklabels=[p.split(' R-HSA')[0][:45] for p in present],
               cbar_kws={'label':'z-score'}, linewidths=0.1)
hm.fig.suptitle('Top 30 ssGSEA pathways — pre-treatment', y=1.02, fontsize=12, fontweight='bold')
hm.savefig(OUT/'Fig4D_ssGSEA_heatmap.pdf', bbox_inches='tight')
hm.savefig(OUT/'Fig4D_ssGSEA_heatmap.png', dpi=300, bbox_inches='tight')
plt.close('all')
print('  ✓ Fig4D_ssGSEA_heatmap')

# 4E — Key DNA repair / cell cycle pathway triple plot
fig, axes = plt.subplots(1, 3, figsize=(10,3.5))
for ax, pw_match in zip(axes, ['DNA Double-Strand Break Repair','HDR Thru Homologous Recombination','G2-M Checkpoint']):
    candidates = [c for c in ssg_m.columns if pw_match in c]
    if not candidates: ax.set_visible(False); continue
    pw = candidates[0]
    pre = ssg_m[ssg_m.timepoint=='pre']
    boxstrip(pre, 'response_bin', pw, ax, title=pw_match)
save_panel(fig, 'Fig4E_key_pathways', OUT)

print('=== Fig 4 done ===')

# =========================================================
# FIGURE 5 — Integration & ML
# =========================================================
# 5A — Feature correlation heatmap (already exists but regenerate clean)
fig, ax = plt.subplots(figsize=(10,9))
feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
num = integ[feats_num].apply(pd.to_numeric, errors='coerce').dropna(axis=0, how='any')
if len(num)>5:
    corr = num.corr(method='spearman')
    # Cluster
    from scipy.cluster.hierarchy import linkage, leaves_list
    Z = linkage(corr.values, method='average')
    idx = leaves_list(Z)
    corr_r = corr.iloc[idx, idx]
    mask = np.triu(np.ones_like(corr_r, dtype=bool), k=1)
    sns.heatmap(corr_r, cmap='RdBu_r', center=0, vmin=-1, vmax=1, ax=ax, mask=mask,
                xticklabels=True, yticklabels=True, cbar_kws={'label':'Spearman ρ','shrink':0.5},
                linewidths=0.1)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=6, rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6, rotation=0)
    ax.set_title('Integrated feature correlation (Spearman, 35 subjects × 37 features)')
save_panel(fig, 'Fig5A_feature_correlation', OUT)

# 5B — Response feature ranking
fig, ax = plt.subplots(figsize=(7,7))
top = rfeat.head(20).iloc[::-1]
colors = [GOOD if d>0 else BAD for d in top.delta_med]
ax.barh(range(len(top)), -np.log10(top.pvalue), color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f[:50] for f in top.feature], fontsize=8)
ax.axvline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('-log10(p-value)')
ax.set_title('Top 20 response-associated features\n(green=↑good, orange=↑bad)')
save_panel(fig, 'Fig5B_feature_ranking', OUT)

# 5C — ML LOOCV AUC
fig, ax = plt.subplots(figsize=(4,3.5))
bars = ax.bar(ml.model, ml.LOOCV_AUC, color=['#264653','#2a9d8f','#e9c46a'], edgecolor='black', linewidth=1)
for bar, v in zip(bars, ml.LOOCV_AUC):
    ax.text(bar.get_x()+bar.get_width()/2, v+0.01, f'{v:.3f}', ha='center', fontsize=9, fontweight='bold')
ax.axhline(0.5, color='red', ls='--', lw=0.8, alpha=0.6, label='Random (0.5)')
ax.set_ylim(0, 1)
ax.set_ylabel('LOOCV AUC')
ax.set_title('Machine-learning response prediction\n(37-feature integrated table)')
ax.legend(fontsize=8)
save_panel(fig, 'Fig5C_ML_AUC', OUT)

# 5D — RF feature importance
fig, ax = plt.subplots(figsize=(5,7))
top_imp = rf_imp.head(15).iloc[::-1]
ax.barh(range(len(top_imp)), top_imp.importance, color='#1d3557', edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top_imp)))
ax.set_yticklabels([f[:50] for f in top_imp.feature], fontsize=8)
ax.set_xlabel('Random Forest importance')
ax.set_title('Top 15 features (Random Forest)')
save_panel(fig, 'Fig5D_RF_importance', OUT)

# 5E — Scatter of two top features, colored by response
fig, ax = plt.subplots(figsize=(4.5,4.5))
top2 = rfeat.head(2).feature.tolist()
if len(top2)==2 and all(t in integ.columns for t in top2):
    for resp in ['good','bad']:
        sub = integ[integ.response_bin==resp]
        ax.scatter(sub[top2[0]], sub[top2[1]], color=PALETTE[resp], s=70, alpha=0.8,
                   edgecolor='black', linewidth=0.6, label=resp)
    ax.set_xlabel(top2[0][:40])
    ax.set_ylabel(top2[1][:40])
    ax.set_title('Top two response-associated features')
    ax.legend(fontsize=9)
save_panel(fig, 'Fig5E_top2_scatter', OUT)

print('=== Fig 5 done ===')

# =========================================================
# FIGURE 6 — Treatment-induced paired delta
# =========================================================
# 6A — TRUST4 BCR IGH delta boxplot
trust_m = trust.copy() if {'subject_id','timepoint','response_bin'}.issubset(trust.columns) else trust.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')
def compute_delta(df, metric):
    rows=[]
    for s in df.subject_id.unique():
        sub = df[df.subject_id==s]
        pre = sub[sub.timepoint=='pre'][metric].values
        post = sub[sub.timepoint=='post'][metric].values
        if len(pre)>0 and len(post)>0:
            rows.append({'subject_id':s,'response':sub.response_bin.iloc[0],'delta':post[0]-pre[0]})
    return pd.DataFrame(rows)

for panel_letter, metric, title in [('A','IGH_n','BCR (IGH) clonotypes Δ'),
                                     ('B','IGH_shannon','IGH Shannon diversity Δ'),
                                     ('C','TRB_shannon','TCR (TRB) Shannon Δ')]:
    fig, ax = plt.subplots(figsize=(3.5,3.5))
    dd = compute_delta(trust_m, metric)
    boxstrip(dd, 'response', 'delta', ax, title=title, ylabel='Δ (post − pre)')
    ax.axhline(0, color='gray', ls='--', lw=0.5)
    save_panel(fig, f'Fig6{panel_letter}_TRUST4_{metric}', OUT)

# 6D — 22 immune sigs delta top
fig, ax = plt.subplots(figsize=(6,5))
top = delta_22.head(12).iloc[::-1]
colors = [GOOD if d>0 else BAD for d in top.delta_good_minus_bad]
ax.barh(range(len(top)), -np.log10(top.MW_p), color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f.replace('_','\n') for f in top.feature], fontsize=7)
ax.axvline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('-log10(p) Δ(post-pre) good vs bad')
ax.set_title('Top 12 immune signatures by paired Δ')
save_panel(fig, 'Fig6D_22sigs_delta', OUT)

# 6E — ssGSEA pathway delta top
fig, ax = plt.subplots(figsize=(7,6))
top = delta_ssg.head(15).iloc[::-1]
colors = [GOOD if d>0 else BAD for d in top.delta_good_minus_bad]
ax.barh(range(len(top)), -np.log10(top.MW_p), color=colors, edgecolor='black', linewidth=0.5)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f[:50] for f in top.feature], fontsize=7)
ax.axvline(-np.log10(0.05), color='gray', ls='--', lw=0.8)
ax.set_xlabel('-log10(p) Δ(post-pre)')
ax.set_title('Top 15 ssGSEA pathways by paired Δ')
save_panel(fig, 'Fig6E_ssGSEA_delta', OUT)

# 6F — SBS delta
fig, ax = plt.subplots(figsize=(3.5,3.5))
sbs_m_delta = pd.read_csv(ROOT/'09_integration/paired_delta/delta_sbs_response.tsv', sep='\t')
# Use SBS5 as representative mutation signature
tmb_m2 = tmb.copy()
rows_m=[]
for s in tmb_m2.subject_id.unique():
    sub = tmb_m2[tmb_m2.subject_id==s]
    pre = sub[sub.timepoint=='pre']
    post = sub[sub.timepoint=='post']
    if len(pre)>0 and len(post)>0 and s not in UNMATCHED:
        rows_m.append({'subject_id':s,'response':sub.response_bin.iloc[0],
                       'delta':post.TMB_nonsyn_per_Mb.iloc[0] - pre.TMB_nonsyn_per_Mb.iloc[0]})
D = pd.DataFrame(rows_m)
boxstrip(D, 'response', 'delta', ax, title='TMB Δ (post − pre)', ylabel='Δ TMB (/Mb)')
ax.axhline(0, color='gray', ls='--', lw=0.5)
save_panel(fig, 'Fig6F_TMB_delta', OUT)

# 6G — Neoantigen binders delta
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(neo_delta, 'response', 'delta_binders', ax, title='MHC-I binders Δ\n(post − pre)', ylabel='Δ binders (<500 nM)')
ax.axhline(0, color='gray', ls='--', lw=0.5)
save_panel(fig, 'Fig6G_neoantigen_delta_binders', OUT)

# 6H — Neoantigen sites delta
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(neo_delta, 'response', 'delta_sites', ax, title='Neoantigen mutation sites Δ', ylabel='Δ sites with binder')
ax.axhline(0, color='gray', ls='--', lw=0.5)
save_panel(fig, 'Fig6H_neoantigen_delta_sites', OUT)

print('=== Fig 6 done ===')

# =========================================================
# FIGURE 7 — External validation
# =========================================================
# 7A — Per-cohort forest plot (effect sizes for each signature)
fig, axes = plt.subplots(1, 4, figsize=(14,4), sharey=True)
sig_names = ['DSB_HDR_repair','E2F_MYC_cellcycle','CD8_proliferation','EMT']
for ax, sig in zip(axes, sig_names):
    sub = ext_stats[ext_stats.signature==sig].sort_values('delta')
    y = range(len(sub))
    colors_ = [GOOD if d>0 else BAD for d in sub.delta]
    ax.scatter(sub.delta, y, s=60, color=colors_, edgecolor='black', linewidth=0.7)
    ax.axvline(0, color='black', lw=0.8)
    ax.set_yticks(y); ax.set_yticklabels(sub.gse, fontsize=8)
    ax.set_xlabel('Δ good − bad (z)')
    ax.set_title(sig.replace('_',' '), fontsize=10)
axes[0].set_ylabel('Cohort')
plt.suptitle('External validation — effect sizes per cohort (n=7 GEO LARC/CRC cohorts)', fontsize=11, fontweight='bold')
plt.tight_layout()
save_panel(fig, 'Fig7A_external_forest', OUT)

# 7B — Meta-analysis Z score
fig, ax = plt.subplots(figsize=(5,3.5))
colors = ['#264653' if z>0 else BAD for z in ext_meta.Z_stouffer]
bars = ax.barh(range(len(ext_meta)), ext_meta.Z_stouffer, color=colors, edgecolor='black', linewidth=0.7)
ax.set_yticks(range(len(ext_meta)))
ax.set_yticklabels(ext_meta.signature)
for i,(v,p) in enumerate(zip(ext_meta.Z_stouffer, ext_meta.p_meta_onesided)):
    ax.text(v + (0.05 if v>0 else -0.05), i, f'Z={v:.2f}, p={p:.2f}', va='center', fontsize=8,
            ha='left' if v>0 else 'right')
ax.axvline(0, color='black', lw=0.8); ax.axvline(1.96, color='red', ls='--', lw=0.5, alpha=0.5)
ax.axvline(-1.96, color='red', ls='--', lw=0.5, alpha=0.5)
ax.set_xlabel('Stouffer Z (expected direction)')
ax.set_title('Meta-analysis across 7 GEO cohorts')
save_panel(fig, 'Fig7B_meta_Zscore', OUT)

# 7C — GSE150082 example deep dive (strongest disagreement)
fig, ax = plt.subplots(figsize=(4,3.5))
gse_f = ROOT/'11_external_validation/GSE150082_signature_scores.tsv'
if gse_f.exists():
    g = pd.read_csv(gse_f, sep='\t', index_col=0)
    if 'response_bin' in g.columns:
        sub = g[g.response_bin.isin(['good','bad'])]
        boxstrip(sub, 'response_bin', 'DSB_HDR_repair', ax, title='GSE150082 — DSB/HDR\n(contradicts discovery)')
save_panel(fig, 'Fig7C_GSE150082_DSB', OUT)

print('=== Fig 7 done ===')

# =========================================================
# FIGURE 8 — HLA + neoantigen
# =========================================================
# 8A — HLA allele frequency (top by locus)
fig, axes = plt.subplots(1, 3, figsize=(13,4))
for ax, locus in zip(axes, ['A','B','C']):
    freq = pd.concat([hla[f'{locus}1'], hla[f'{locus}2']]).value_counts().head(8)
    ax.barh(range(len(freq)), freq.values, color='#1d3557', edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(freq)))
    ax.set_yticklabels([f.replace('HLA-','') for f in freq.index], fontsize=9)
    ax.set_xlabel('Allele count')
    ax.set_title(f'HLA-{locus} top alleles', fontsize=10)
plt.suptitle('HLA class I allele frequency (35 subjects)', fontsize=12, fontweight='bold')
plt.tight_layout()
save_panel(fig, 'Fig8A_HLA_alleles', OUT)

# 8B — HLA homozygosity by response
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(hla, 'response_bin', 'n_homozygous_loci', ax, title='HLA homozygous loci\n(class I)', ylabel='Homozygous loci (0-3)')
save_panel(fig, 'Fig8B_HLA_homozygosity', OUT)

# 8C — HLA LOH count by response
loh_summary = loh.groupby(['subject_id','sample'])['LOH_call'].sum().reset_index()
loh_summary = loh_summary.merge(wes_inv[['sample_id','timepoint','response_bin']].rename(columns={'sample_id':'sample'}), on='sample')
loh_pre = loh_summary[loh_summary.timepoint=='pre']
fig, ax = plt.subplots(figsize=(4,3.5))
boxstrip(loh_pre, 'response_bin', 'LOH_call', ax, title='HLA class I LOH events\n(pre-treatment)', ylabel='LOH loci (0-3)')
save_panel(fig, 'Fig8C_HLA_LOH', OUT)

# 8D — Neoantigen binders pre
fig, ax = plt.subplots(figsize=(3.5,3.5))
neo_m = neo[neo.matched & (neo.timepoint=='pre')]
boxstrip(neo_m, 'response_bin', 'n_sites_with_binder', ax, title='Neoantigen mutation sites (pre)', ylabel='Sites with MHC-I binder (<500 nM)')
save_panel(fig, 'Fig8D_neoantigen_pre_sites', OUT)

# 8E — Strong binders pre
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(neo_m, 'response_bin', 'n_strong_binders_50nM', ax, title='Strong MHC-I binders (<50 nM)', ylabel='Strong binder count')
save_panel(fig, 'Fig8E_strong_binders', OUT)

# 8F — PCN score
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(neo_m, 'response_bin', 'PCN_score', ax, title='Presentation-competent\nneoantigen score', ylabel='PCN score')
save_panel(fig, 'Fig8F_PCN_score', OUT)

print('=== Fig 8 done ===')

# =========================================================
# FIGURE 9 — Clonal evolution
# =========================================================
# 9A — Clusters per subject
fig, ax = plt.subplots(figsize=(5,4))
py_plot = pyclone.sort_values(['response','subject_id'])
colors = [PALETTE[r] for r in py_plot.response]
ax.bar(range(len(py_plot)), py_plot.n_clusters, color=colors, edgecolor='black', linewidth=0.5)
ax.set_xticks(range(len(py_plot)))
ax.set_xticklabels(py_plot.subject_id, fontsize=9)
ax.set_ylabel('Number of clonal clusters')
ax.set_xlabel('Subject ID')
ax.set_title('PyClone-VI clonal clusters per subject (paired pre+post)')
# response legend
handles = [plt.Rectangle((0,0),1,1, color=GOOD, label='good'),
           plt.Rectangle((0,0),1,1, color=BAD, label='bad')]
ax.legend(handles=handles, fontsize=9, loc='upper right')
save_panel(fig, 'Fig9A_pyclone_clusters', OUT)

# 9B — Dominant clone shrinkage
fig, ax = plt.subplots(figsize=(3.5,3.5))
boxstrip(pyclone, 'response', 'dominant_shrink', ax, title='Dominant clone shrinkage\n(min Δ CP pre→post)', ylabel='Cellular prevalence Δ')
ax.axhline(0, color='gray', ls='--', lw=0.5)
save_panel(fig, 'Fig9B_dominant_shrink', OUT)

# 9C — Clonal evolution scatter: n_shrinking vs n_expanding
fig, ax = plt.subplots(figsize=(4.5,4))
for resp in ['good','bad']:
    sub = pyclone[pyclone.response==resp]
    ax.scatter(sub.n_shrinking, sub.n_expanding, color=PALETTE[resp], s=100, alpha=0.8,
               edgecolor='black', linewidth=0.8, label=resp)
    for _, r in sub.iterrows():
        ax.annotate(f's{int(r.subject_id)}', (r.n_shrinking, r.n_expanding), fontsize=8, xytext=(3,3), textcoords='offset points')
ax.set_xlabel('Shrinking clusters (Δ CP < −0.2)')
ax.set_ylabel('Expanding clusters (Δ CP > +0.2)')
ax.set_title('Clonal evolution pattern per subject')
ax.legend(fontsize=9)
save_panel(fig, 'Fig9C_clonal_scatter', OUT)

print('=== Fig 9 done ===')
print('\n=== ALL FIGURE PANELS GENERATED ===')
print(f'Output directory: {OUT}')
panels = sorted(OUT.glob('Fig*.png'))
print(f'Total: {len(panels)} PNG panels + matching PDFs')
for p in panels: print('  ', p.name)
