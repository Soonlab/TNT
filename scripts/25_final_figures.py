"""
Generate final figure panels: Fig1 cohort, Fig6 neoantigen + paired delta, updated Fig5 integration.
"""
import pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
from pathlib import Path
from matplotlib.gridspec import GridSpec
from scipy import stats

ROOT = '/mnt/sda1/data/TNT/analysis'
FIG = Path(f'{ROOT}/figures/panels'); FIG.mkdir(parents=True, exist_ok=True)

# ---------- Load data ----------
clin = pd.read_csv(f'{ROOT}/00_cohort/clinical_master.tsv', sep='\t')
wes_inv = pd.read_csv(f'{ROOT}/00_cohort/wes_inventory.tsv', sep='\t')
rna_inv = pd.read_csv(f'{ROOT}/00_cohort/rna_inventory.tsv', sep='\t')
tmb = pd.read_csv(f'{ROOT}/02_wes_tmb_msi/tmb_per_sample.tsv', sep='\t')
msi = pd.read_csv(f'{ROOT}/02_wes_tmb_msi/msi/msi_summary_paired.tsv', sep='\t')
neo = pd.read_csv(f'{ROOT}/03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv', sep='\t')
trust = pd.read_csv(f'{ROOT}/06_rna_immune/trust4_summary.tsv', sep='\t')
delta_trust = pd.read_csv(f'{ROOT}/09_integration/paired_delta/delta_trust4_response.tsv', sep='\t')
neo_delta = pd.read_csv(f'{ROOT}/03_wes_hla_neoantigen/neoantigen_paired_delta.tsv', sep='\t')

# ========================================================
# Figure 1 — Cohort overview
# ========================================================
fig = plt.figure(figsize=(14, 8))
gs = GridSpec(3, 3, figure=fig, hspace=0.6, wspace=0.4)

# 1A: pie response
ax = fig.add_subplot(gs[0,0])
rcounts = clin.response_bin.value_counts()
ax.pie(rcounts, labels=[f'{k} (n={v})' for k,v in rcounts.items()],
       colors=['#2a9d8f' if k=='good' else '#e76f51' for k in rcounts.index], autopct='%1.0f%%')
ax.set_title('A. Response distribution\n(TRG0-1 vs TRG2-3)', fontsize=11, fontweight='bold')

# 1B: stage stack
ax = fig.add_subplot(gs[0,1])
ct_crosstab = pd.crosstab(clin.cT, clin.response_bin, normalize='columns')*100
ct_crosstab = ct_crosstab[['good','bad']] if 'good' in ct_crosstab else ct_crosstab
bottom = np.zeros(2)
colors=['#264653','#2a9d8f','#e9c46a','#e76f51']
for i,t in enumerate(['T2','T2/T3','T3','T4']):
    if t in ct_crosstab.index:
        ax.bar(range(2), ct_crosstab.loc[t].values, bottom=bottom, label=t, color=colors[i])
        bottom += ct_crosstab.loc[t].values
ax.set_xticks(range(2)); ax.set_xticklabels(['good','bad'])
ax.set_ylabel('Percentage'); ax.legend(title='cT stage', fontsize=8, loc='upper right')
ax.set_title('B. Clinical T stage by response\n(p=0.086)', fontsize=11, fontweight='bold')

# 1C: age/sex
ax = fig.add_subplot(gs[0,2])
for i, resp in enumerate(['good','bad']):
    for j, sex in enumerate(['M','F']):
        sub = clin[(clin.response_bin==resp) & (clin.sex==sex)]
        ax.scatter([i+(j-0.5)*0.3]*len(sub), sub.age, alpha=0.7,
                   color='#2a9d8f' if resp=='good' else '#e76f51',
                   marker='o' if sex=='M' else '^',
                   s=50)
ax.set_xticks([0,1]); ax.set_xticklabels(['good','bad'])
ax.set_ylabel('Age (years)'); ax.set_title('C. Age × sex (○M △F)', fontsize=11, fontweight='bold')

# 1D: sample matrix
ax = fig.add_subplot(gs[1,:])
# Build per-subject matrix: pre/post/normal × WES/RNA
subs = sorted(clin.subject_id)
mat = np.zeros((len(subs), 6))
for i, s in enumerate(subs):
    w = wes_inv[wes_inv.subject_id==s]
    r = rna_inv[rna_inv.subject_id==s]
    mat[i,0] = len(w[w.timepoint=='normal'])  # WES normal
    mat[i,1] = len(w[w.timepoint=='pre'])
    mat[i,2] = len(w[w.timepoint=='post'])
    mat[i,3] = len(r[r.timepoint=='normal'])  # RNA normal
    mat[i,4] = len(r[r.timepoint=='pre'])
    mat[i,5] = len(r[r.timepoint=='post'])
sns.heatmap(mat.T, cmap='Greys', cbar=False, ax=ax,
            yticklabels=['WES N','WES pre','WES post','RNA N','RNA pre','RNA post'],
            xticklabels=[str(s) for s in subs],
            linewidths=0.3, linecolor='white')
ax.set_xticklabels(subs, rotation=0, fontsize=7)
ax.set_title('D. Sample matrix per subject (black=present)', fontsize=11, fontweight='bold')
# add response color annotation below
resp_colors = [('#2a9d8f' if clin[clin.subject_id==s].response_bin.iloc[0]=='good' else '#e76f51') for s in subs]
for i, c in enumerate(resp_colors):
    ax.add_patch(plt.Rectangle((i, 6.2), 1, 0.2, color=c, clip_on=False))

# 1E: TMB + MSI summary
ax = fig.add_subplot(gs[2,0])
pre_m = tmb[(tmb.timepoint=='pre') & (~tmb.subject_id.isin([13,15,16,17,18,19,33]))]
sns.boxplot(data=pre_m, x='response_bin', y='TMB_nonsyn_per_Mb', ax=ax, order=['good','bad'],
            hue='response_bin', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
sns.stripplot(data=pre_m, x='response_bin', y='TMB_nonsyn_per_Mb', ax=ax, order=['good','bad'],
              color='black', size=4)
ax.axhline(10, color='gray', ls='--', label='TMB-high (10/Mb)')
ax.set_title('E. TMB pre (matched)\np=0.186', fontsize=10, fontweight='bold')

# 1F: MSI
ax = fig.add_subplot(gs[2,1])
sns.boxplot(data=msi[msi.timepoint=='pre'], x='response_bin', y='MSI_pct', ax=ax, order=['good','bad'],
            hue='response_bin', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
sns.stripplot(data=msi[msi.timepoint=='pre'], x='response_bin', y='MSI_pct', ax=ax, order=['good','bad'],
              color='black', size=4)
ax.axhline(20, color='red', ls='--', label='MSI-H (>20%)')
ax.set_title('F. MSI status\nall <0.2% (MSS)', fontsize=10, fontweight='bold')

# 1G: Neoantigen pre (this study's new deliverable)
ax = fig.add_subplot(gs[2,2])
neo_pre = neo[(neo.timepoint=='pre') & (neo.matched)]
sns.boxplot(data=neo_pre, x='response_bin', y='n_sites_with_binder', ax=ax, order=['good','bad'],
            hue='response_bin', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
sns.stripplot(data=neo_pre, x='response_bin', y='n_sites_with_binder', ax=ax, order=['good','bad'],
              color='black', size=4)
ax.set_title('G. Neoantigen sites (pre)\ngood 73 vs bad 66, p=0.082', fontsize=10, fontweight='bold')
ax.set_ylabel('Sites with MHC-I binder (<500nM)')

plt.suptitle('Fig 1. Cohort overview and basic molecular landscape', fontsize=13, fontweight='bold', y=1.00)
plt.savefig(FIG/'Fig1_cohort.png', dpi=150, bbox_inches='tight')
plt.savefig(FIG/'Fig1_cohort.pdf', bbox_inches='tight')
print('Saved Fig1_cohort')

# ========================================================
# Figure 6 — Paired pre→post immune + neoantigen
# ========================================================
fig, axes = plt.subplots(2, 4, figsize=(16, 8))

trust_m = trust.copy() if {'subject_id','timepoint','response_bin'}.issubset(trust.columns) else trust.merge(rna_inv[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')
# Filter paired
subs_rna = set(trust_m.groupby('subject_id').filter(lambda g: set(g.timepoint)>={'pre','post'}).subject_id)

# Top row: TRUST4 paired deltas for TRB + IGH
for ax, metric, label in zip(axes[0], ['IGH_n','IGH_shannon','TRB_shannon','TRB_n'],
                              ['BCR (IGH) clonotypes','IGH Shannon diversity','TCR (TRB) Shannon','TRB clonotypes']):
    sub = trust_m[trust_m.subject_id.isin(subs_rna)]
    deltas=[]
    for s in sorted(subs_rna):
        ss = sub[sub.subject_id==s]
        pre = ss[ss.timepoint=='pre'][metric].values
        post = ss[ss.timepoint=='post'][metric].values
        resp = ss.response_bin.iloc[0]
        if len(pre)>0 and len(post)>0:
            deltas.append({'subject':s,'response':resp,'delta':post[0]-pre[0]})
    dd = pd.DataFrame(deltas)
    sns.boxplot(data=dd, x='response', y='delta', ax=ax, order=['good','bad'],
                hue='response', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
    sns.stripplot(data=dd, x='response', y='delta', ax=ax, order=['good','bad'], color='black', size=5)
    ax.axhline(0, color='gray', ls='--')
    row = delta_trust[delta_trust.feature==metric]
    p = row.MW_p.values[0] if len(row) else np.nan
    ax.set_title(f'{label}\nΔ (post-pre)  p={p:.3f}', fontsize=10)
    ax.set_xlabel('')

# Bottom row: neoantigen delta
for ax, metric, label in zip(axes[1][:4], ['delta_binders','delta_sites','delta_strong','delta_PCN'],
                              ['Δ MHC-I binders','Δ mutation sites (binder)','Δ strong binders','Δ PCN score']):
    sns.boxplot(data=neo_delta, x='response', y=metric, ax=ax, order=['good','bad'],
                hue='response', palette={'good':'#2a9d8f','bad':'#e76f51'}, legend=False)
    sns.stripplot(data=neo_delta, x='response', y=metric, ax=ax, order=['good','bad'], color='black', size=5)
    ax.axhline(0, color='gray', ls='--')
    # p from MW
    g = neo_delta[neo_delta.response=='good'][metric].dropna()
    b = neo_delta[neo_delta.response=='bad'][metric].dropna()
    if len(g)>=3 and len(b)>=3:
        p = stats.mannwhitneyu(g,b).pvalue
    else: p = np.nan
    ax.set_title(f'{label}\nΔ (post-pre)  p={p:.3f}', fontsize=10)
    ax.set_xlabel('')

plt.suptitle('Fig 6. Treatment-induced immune activation and neoantigen clearance in good responders', fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(FIG/'Fig6_paired_delta.png', dpi=150, bbox_inches='tight')
plt.savefig(FIG/'Fig6_paired_delta.pdf', bbox_inches='tight')
print('Saved Fig6_paired_delta')

print('\n=== Final figures in', FIG, '===')
for f in sorted(FIG.glob('*.png')): print(' ', f.name)
