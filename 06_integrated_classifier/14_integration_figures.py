"""
Multi-omics integration + key figures + summary tables for manuscript.
Combines: WES (TMB, SBS, MSI, driver, CNV/CIN, HLA), RNA (sigs, CMS, ssGSEA, DEG)
"""
import pandas as pd, numpy as np, os
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT='/mnt/sda1/data/TNT/analysis'
FIG=Path(f'{ROOT}/figures'); FIG.mkdir(exist_ok=True)
TAB=Path(f'{ROOT}/tables'); TAB.mkdir(exist_ok=True)

# Load all data
clin = pd.read_csv(f'{ROOT}/00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(f'{ROOT}/00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(f'{ROOT}/00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(f'{ROOT}/02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(f'{ROOT}/02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
cnv = pd.read_csv(f'{ROOT}/04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
sbs = pd.read_csv(f'{ROOT}/01_wes_signatures/sbs_summary_key.tsv', sep='\t')
hla = pd.read_csv(f'{ROOT}/03_hla/hla_class_I_typing.tsv', sep='\t')
cms = pd.read_csv(f'{ROOT}/07_rna_cms/cms_assignments.tsv', sep='\t')
sigs= pd.read_csv(f'{ROOT}/06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
ssg = pd.read_csv(f'{ROOT}/08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0)
ssg_res = pd.read_csv(f'{ROOT}/08_rna_pathway/ssgsea_response_stats.tsv', sep='\t')
deg_gsea_h = pd.read_csv(f'{ROOT}/05_rna_deg_gsea/GSEA_Hallmark_pre.tsv', sep='\t')

# Build per-subject (pre-treatment tumor) integrated table
pre_wes = tmb[tmb.timepoint=='pre'][['subject_id','TMB_nonsyn_per_Mb','n_nonsyn']]
pre_cnv = cnv[cnv.timepoint=='pre'][['subject_id','CIN','frac_amp','frac_del']]
pre_msi = msi[msi.timepoint=='pre'][['subject_id','MSI_pct']]
pre_sbs = sbs[sbs.timepoint=='pre'][['subject_id','MMR_prop','SBS5']]

integ = clin[['subject_id','response_bin','response_num','age','sex','cT','prepost_set']].copy()
integ = integ.merge(pre_wes, on='subject_id', how='left')
integ = integ.merge(pre_cnv, on='subject_id', how='left')
integ = integ.merge(pre_msi, on='subject_id', how='left')
integ = integ.merge(pre_sbs, on='subject_id', how='left')
# RNA-level pre signature scores
pre_rna_ids = rna_inv[rna_inv.timepoint=='pre'][['sample_id','subject_id']]
pre_sig = sigs.loc[pre_rna_ids.sample_id.values].reset_index().merge(pre_rna_ids, on='sample_id')
key_sigs = ['CD8_proliferation','CD8_activation','MHC_II','Antigen_presentation',
            'NLRC5_HLA_IFNG','TLS_Cabrita','IFNg_Ayers_18','TGFb_Mariathasan','EMT_Mak',
            'Hypoxia_Buffa','Stemness_mRNAsi_proxy']
pre_sig_k = pre_sig[['subject_id']+key_sigs].set_index('subject_id')
integ = integ.merge(pre_sig_k, on='subject_id', how='left')

# ssGSEA key pathways
ssg_pre = ssg.loc[pre_rna_ids.sample_id.values]
ssg_pre = ssg_pre.merge(pre_rna_ids, left_index=True, right_on='sample_id').set_index('subject_id')
key_pw = ['DNA Double-Strand Break Repair R-HSA-5693532','HDR Thru Homologous Recombination (HRR) R-HSA-5685942',
          'DNA Repair R-HSA-73894','E2F Targets','G2-M Checkpoint','Myc Targets V2',
          'Epithelial Mesenchymal Transition','TGF-beta Signaling','Hypoxia']
pw_exist = [p for p in key_pw if p in ssg_pre.columns]
ssg_pre_k = ssg_pre[pw_exist]
integ = integ.merge(ssg_pre_k, left_on='subject_id', right_index=True, how='left')

# CMS
cms_pre = cms[cms.timepoint=='pre'][['subject_id','prediction']].rename(columns={'prediction':'CMS'})
integ = integ.merge(cms_pre, on='subject_id', how='left')

integ['matched_wes'] = ~integ['subject_id'].isin([13,15,16,17,18,19,33])
integ.to_csv(f'{TAB}/integrated_subject_master.tsv', sep='\t', index=False)
print(f'Integrated: {len(integ)} subjects, {integ.shape[1]} features')

# --- Response association volcano: all features vs response ---
feats_num = [c for c in integ.columns if c not in ['subject_id','response_bin','response_num','sex','cT','prepost_set','CMS','matched_wes']]
rows=[]
for f in feats_num:
    x = pd.to_numeric(integ[f], errors='coerce')
    gg = x[integ.response_bin=='good'].dropna()
    bb = x[integ.response_bin=='bad'].dropna()
    if len(gg)<3 or len(bb)<3: continue
    u = stats.mannwhitneyu(gg, bb)
    rows.append((f, len(gg), len(bb), gg.median(), bb.median(), gg.median()-bb.median(), float(u.pvalue)))
F = pd.DataFrame(rows, columns=['feature','n_good','n_bad','med_good','med_bad','delta_med','pvalue'])
_,q,_,_=multipletests(F['pvalue'],method='fdr_bh')
F['qvalue']=q
F=F.sort_values('pvalue')
F.to_csv(f'{TAB}/response_feature_stats.tsv', sep='\t', index=False)
print('\n=== Response-associated features (p<0.1) ===')
print(F[F.pvalue<0.1].to_string(index=False))

# Volcano-style bar: -log10(p) sorted, color by direction
fig, ax = plt.subplots(figsize=(9,7))
F['neglog10p'] = -np.log10(F['pvalue'])
F['dir'] = np.sign(F['delta_med'])
top = F.head(25).iloc[::-1]
colors = ['#2a9d8f' if d>0 else '#e76f51' for d in top['dir']]
ax.barh(range(len(top)), top['neglog10p'], color=colors)
ax.set_yticks(range(len(top))); ax.set_yticklabels(top['feature'], fontsize=9)
ax.axvline(-np.log10(0.05), color='gray', linestyle='--', label='p=0.05')
ax.set_xlabel('-log10(p) good vs bad')
ax.set_title('Response association — integrated features (teal=↑good, orange=↑bad)')
ax.legend()
plt.tight_layout()
plt.savefig(f'{FIG}/Fig_response_features_barh.png', dpi=150, bbox_inches='tight')
print('Saved Fig_response_features_barh.png')

# Correlation heatmap of key features
num = integ[feats_num].apply(pd.to_numeric, errors='coerce').dropna(axis=0, how='any')
if len(num)>5:
    corr = num.corr(method='spearman')
    fig, ax = plt.subplots(figsize=(14,12))
    sns.heatmap(corr, cmap='RdBu_r', center=0, vmin=-1, vmax=1, ax=ax, xticklabels=True, yticklabels=True, cbar_kws={'shrink':.5})
    plt.xticks(rotation=90, fontsize=8); plt.yticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(f'{FIG}/Fig_feature_correlation_heatmap.png', dpi=150, bbox_inches='tight')
    print('Saved Fig_feature_correlation_heatmap.png')

# --- Summary GSEA bar (Hallmark) ---
deg_gsea_h['signedLogP'] = -np.log10(deg_gsea_h['pval']) * np.sign(deg_gsea_h['NES'])
top_path = deg_gsea_h.nsmallest(20, 'pval').copy()
top_path = top_path.iloc[::-1]
fig, ax = plt.subplots(figsize=(9,7))
colors = ['#2a9d8f' if n>0 else '#e76f51' for n in top_path['NES']]
ax.barh(range(len(top_path)), top_path['signedLogP'], color=colors)
ax.set_yticks(range(len(top_path))); ax.set_yticklabels(top_path['pathway'].str.replace('HALLMARK_',''), fontsize=9)
ax.axvline(0, color='black')
ax.set_xlabel('signed -log10(p) [+=UP in good]'); ax.set_title('GSEA Hallmark — good vs bad (pre)')
plt.tight_layout()
plt.savefig(f'{FIG}/Fig_GSEA_Hallmark.png', dpi=150, bbox_inches='tight')
print('Saved Fig_GSEA_Hallmark.png')

print('\n--- ALL key figures/tables generated ---')
