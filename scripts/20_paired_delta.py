"""
Paired pre→post delta analysis: treatment-induced change stratified by response.
Tests whether good vs bad responders show DIFFERENT magnitudes of post-pre change
across molecular features.
Features: TRUST4 TCR/BCR, RNA signatures, ssGSEA, immune deconv, TMB, SBS, CNV/CIN
"""
import pandas as pd, numpy as np, os
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns

ROOT='/mnt/sda1/data/TNT/analysis'
OUT = Path(f'{ROOT}/09_integration/paired_delta'); OUT.mkdir(parents=True, exist_ok=True)
FIG = Path(f'{ROOT}/figures'); FIG.mkdir(exist_ok=True)

rna_inv = pd.read_csv(f'{ROOT}/00_cohort/rna_inventory.tsv', sep='\t')
wes_inv = pd.read_csv(f'{ROOT}/00_cohort/wes_inventory.tsv', sep='\t')
clin = pd.read_csv(f'{ROOT}/00_cohort/clinical_master.tsv', sep='\t')

# Helper: build pre/post wide per subject for a given feature matrix
def paired_delta(mat, inv, id_col='sample_id'):
    """mat: sample_id x features. inv: sample_id,subject_id,timepoint. Returns subject_id x features delta (post-pre)."""
    m = mat.reset_index().rename(columns={'index':'sample_id', mat.index.name or 'sample_id':'sample_id'})
    if 'sample_id' not in m.columns: m = m.reset_index()
    m = m.merge(inv[['sample_id','subject_id','timepoint']], on='sample_id')
    keep = m[m.timepoint.isin(['pre','post'])]
    subjects = [s for s,g in keep.groupby('subject_id') if set(g.timepoint)=={'pre','post'}]
    numeric = keep.select_dtypes(include=[np.number]).columns.tolist()
    rows=[]
    for s in subjects:
        pre = keep[(keep.subject_id==s) & (keep.timepoint=='pre')].iloc[0]
        post = keep[(keep.subject_id==s) & (keep.timepoint=='post')].iloc[0]
        d = post[numeric] - pre[numeric]
        d['subject_id']=s
        rows.append(d)
    df = pd.DataFrame(rows)
    if 'sample_id' in df.columns: df = df.drop(columns=['sample_id'])
    if 'timepoint' in df.columns: df = df.drop(columns=['timepoint'])
    df = df.set_index('subject_id')
    return df.apply(pd.to_numeric, errors='coerce')

def response_test(delta_df, feature_list=None):
    m = delta_df.reset_index().merge(clin[['subject_id','response_bin']], on='subject_id')
    g = m[m.response_bin=='good']; b = m[m.response_bin=='bad']
    feats = feature_list if feature_list else [c for c in m.columns if c not in ['subject_id','response_bin']]
    rows=[]
    for f in feats:
        gv = g[f].dropna(); bv = b[f].dropna()
        if len(gv)<3 or len(bv)<3: continue
        u = stats.mannwhitneyu(gv, bv)
        # also test "is delta significantly != 0 in each group" (paired feature = pre→post change)
        # Wilcoxon signed-rank equivalent: delta distribution
        try: gp = stats.wilcoxon(gv)
        except: gp = None
        try: bp = stats.wilcoxon(bv)
        except: bp = None
        rows.append({'feature':f,
            'n_good':len(gv),'delta_good_median':float(gv.median()),'good_wilcox_p':None if gp is None else gp.pvalue,
            'n_bad':len(bv),'delta_bad_median':float(bv.median()),'bad_wilcox_p':None if bp is None else bp.pvalue,
            'delta_good_minus_bad':float(gv.mean()-bv.mean()),'MW_p':float(u.pvalue)})
    return pd.DataFrame(rows).sort_values('MW_p')

print('\n########## RNA: 22 immune signatures ##########')
sigs = pd.read_csv(f'{ROOT}/06_rna_immune/signature_scores.tsv', sep='\t', index_col=0)
d_sig = paired_delta(sigs, rna_inv)
print(f'paired subjects: {len(d_sig)}')
r_sig = response_test(d_sig)
r_sig.to_csv(OUT/'delta_22sigs_response.tsv', sep='\t', index=False)
print(r_sig.head(15).to_string(index=False))

print('\n########## RNA: ssGSEA 95 pathways ##########')
ssg = pd.read_csv(f'{ROOT}/08_rna_pathway/ssgsea_scores.tsv', sep='\t', index_col=0)
ssg = ssg.apply(pd.to_numeric, errors='coerce').dropna(axis=1, how='all')
d_ssg = paired_delta(ssg, rna_inv)
print(f'paired subjects: {len(d_ssg)}')
r_ssg = response_test(d_ssg)
r_ssg.to_csv(OUT/'delta_ssgsea_response.tsv', sep='\t', index=False)
print(r_ssg.head(20).to_string(index=False))

print('\n########## RNA: TRUST4 clonality ##########')
trust = pd.read_csv(f'{ROOT}/06_rna_immune/trust4_summary.tsv', sep='\t')
# Keep numeric columns minus subject_id (which we re-add via inv merge)
trust_num = trust.select_dtypes(include=[np.number]).drop(columns=['subject_id'], errors='ignore')
trust_num['sample_id'] = trust['sample_id']
trust_wide = trust_num.set_index('sample_id')
d_trust = paired_delta(trust_wide, rna_inv)
print(f'paired subjects: {len(d_trust)}')
r_trust = response_test(d_trust)
r_trust.to_csv(OUT/'delta_trust4_response.tsv', sep='\t', index=False)
print(r_trust.head(25).to_string(index=False))

print('\n########## WES: TMB + SBS ##########')
tmb = pd.read_csv(f'{ROOT}/02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
tmb_feat = ['n_total','n_nonsyn','TMB_nonsyn_per_Mb']
tmb_wide = tmb.set_index('sample_id')[tmb_feat]
d_tmb = paired_delta(tmb_wide, wes_inv)
r_tmb = response_test(d_tmb)
r_tmb.to_csv(OUT/'delta_tmb_response.tsv', sep='\t', index=False)
print(r_tmb.to_string(index=False))

sbs = pd.read_csv(f'{ROOT}/01_wes_signatures/sbs_summary_key.tsv', sep='\t')
sbs_feat = ['SBS1','SBS5','SBS6','SBS15','MMR_sum','SBS3','SBS40','SBS54']
sbs_wide = sbs.set_index('sample_id')[sbs_feat]
d_sbs = paired_delta(sbs_wide, wes_inv)
r_sbs = response_test(d_sbs)
r_sbs.to_csv(OUT/'delta_sbs_response.tsv', sep='\t', index=False)
print(r_sbs.to_string(index=False))

print('\n########## WES: CNV/CIN ##########')
cnv = pd.read_csv(f'{ROOT}/04_wes_cnv_clonal/cnv_cin_per_sample.tsv', sep='\t')
cnv_feat = ['CIN','frac_amp','frac_del','n_segments']
cnv_wide = cnv.set_index('sample_id')[cnv_feat]
d_cnv = paired_delta(cnv_wide, wes_inv)
r_cnv = response_test(d_cnv)
r_cnv.to_csv(OUT/'delta_cnv_response.tsv', sep='\t', index=False)
print(r_cnv.to_string(index=False))

# HRD proxy
hrd = pd.read_csv(f'{ROOT}/04_wes_cnv_clonal/hrd_proxy/hrd_proxy_scores.tsv', sep='\t')
hrd_feat = ['LST','LOH','TAI','HRD_sum']
hrd_wide = hrd.set_index('sample_id')[hrd_feat]
d_hrd = paired_delta(hrd_wide, wes_inv)
r_hrd = response_test(d_hrd)
r_hrd.to_csv(OUT/'delta_hrd_response.tsv', sep='\t', index=False)
print(r_hrd.to_string(index=False))

# FIGURE: TRUST4 delta boxplots for key metrics
fig, axes = plt.subplots(2,4, figsize=(16,7))
m = d_trust.reset_index().merge(clin[['subject_id','response_bin']], on='subject_id')
for ax, metric in zip(axes.flat, ['TRB_n','TRB_shannon','TRB_gini','TRB_top1',
                                    'IGH_n','IGH_shannon','IGK_n','IGL_n']):
    if metric not in m.columns: ax.set_visible(False); continue
    sns.boxplot(data=m, x='response_bin', y=metric, ax=ax, order=['good','bad'],
                hue='response_bin', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
    sns.stripplot(data=m, x='response_bin', y=metric, ax=ax, order=['good','bad'], color='black', size=5)
    ax.axhline(0, color='gray', ls='--')
    row = r_trust[r_trust.feature==metric]
    p = row.MW_p.values[0] if len(row) else np.nan
    ax.set_title(f'{metric} Δ (post−pre)\np={p:.3f}' if not np.isnan(p) else metric, fontsize=10)
    ax.set_xlabel('')
plt.tight_layout()
plt.savefig(FIG/'Fig_paired_delta_TRUST4.png', dpi=150, bbox_inches='tight')
print('\nSaved Fig_paired_delta_TRUST4.png')

# FIGURE: top ssGSEA deltas
fig, ax = plt.subplots(figsize=(9,7))
top = r_ssg.head(20).iloc[::-1]
colors = ['#2a9d8f' if d>0 else '#e76f51' for d in top.delta_good_minus_bad]
ax.barh(range(len(top)), -np.log10(top.MW_p), color=colors)
ax.set_yticks(range(len(top))); ax.set_yticklabels([s[:60] for s in top.feature], fontsize=8)
ax.axvline(-np.log10(0.05), color='gray', linestyle='--')
ax.set_xlabel('-log10(p) post-pre Δ good vs bad')
ax.set_title('Top treatment-induced delta (teal=higher Δ in good, orange=higher Δ in bad)')
plt.tight_layout()
plt.savefig(FIG/'Fig_paired_delta_ssgsea.png', dpi=150, bbox_inches='tight')
print('Saved Fig_paired_delta_ssgsea.png')

print('\nAll paired-delta outputs in', OUT)
