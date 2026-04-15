"""
Parse pvacseq all_epitopes tsv per sample, count neoantigens by binding strength,
integrate HLA LOH + response, paired pre→post delta analysis.
"""
import pandas as pd, numpy as np, os
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests

P = Path('/mnt/sda1/data/TNT/analysis/03_wes_hla_neoantigen/pvacseq')
OUT = Path('/mnt/sda1/data/TNT/analysis/03_wes_hla_neoantigen')
INV = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')
LOH = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_hla/loh_lite/hla_loh_lite_results.tsv', sep='\t')

rows=[]
for d in sorted(P.iterdir()):
    if not d.is_dir(): continue
    ep_dir = d/'MHC_Class_I'
    all_ep = list(ep_dir.glob('*_DNA.MHC_I.all_epitopes.tsv'))
    if not all_ep: continue
    sid = d.name
    df = pd.read_csv(all_ep[0], sep='\t', low_memory=False)
    if len(df)==0:
        rows.append({'sample_id':sid,'n_mutation_sites':0,'n_candidate_peptides':0,
                     'n_binders_500nM':0,'n_strong_binders_50nM':0,'n_neoantigens_per_site':0}); continue
    # Best-per-mutation-site: pick minimum MHCflurry MT IC50 per (Chromosome, Start, Reference, Variant)
    df['Best IC50'] = pd.to_numeric(df['MHCflurry MT IC50 Score'].replace('NA',np.nan), errors='coerce')
    by_site = df.groupby(['Chromosome','Start','Reference','Variant'])['Best IC50'].min().reset_index()
    n_sites = len(by_site)
    n_pep_total = len(df)
    # Binding cutoffs (MHCflurry MT IC50)
    binders = df[df['Best IC50']<500]
    strong = df[df['Best IC50']<50]
    # unique mutation sites with at least one binder
    n_site_binders = binders[['Chromosome','Start','Reference','Variant']].drop_duplicates().shape[0]
    n_site_strong = strong[['Chromosome','Start','Reference','Variant']].drop_duplicates().shape[0]
    rows.append({'sample_id':sid,
                 'n_mutation_sites':n_sites,
                 'n_candidate_peptides':n_pep_total,
                 'n_binders_500nM':len(binders),
                 'n_strong_binders_50nM':len(strong),
                 'n_sites_with_binder':n_site_binders,
                 'n_sites_with_strong':n_site_strong,
                 'neoantigens_per_site':len(binders)/n_sites if n_sites else 0})

S = pd.DataFrame(rows)
S = S.merge(INV[['sample_id','subject_id','timepoint','response_bin']], on='sample_id')
# HLA LOH: binary per sample (any LOH call)
loh_samp = LOH.groupby('sample')['LOH_call'].any().reset_index()
loh_samp.columns=['sample_id','HLA_LOH']
S = S.merge(loh_samp, on='sample_id', how='left')
S['HLA_LOH'] = S['HLA_LOH'].fillna(False)
S['matched'] = ~S['subject_id'].isin([13,15,16,17,18,19,33])

# Presentation-competent neoantigen: n_sites_with_binder × (1-0.33*LOH)
S['PCN_score'] = S['n_sites_with_binder'] * (1 - S['HLA_LOH'].astype(int)*0.33)

S.to_csv(OUT/'neoantigen_summary_by_sample.tsv', sep='\t', index=False)
print(f'{len(S)} samples parsed')
print(S[['sample_id','response_bin','timepoint','n_mutation_sites','n_sites_with_binder','n_sites_with_strong','HLA_LOH','PCN_score']].head(10).to_string(index=False))

# ------ Response association ------
print('\n=== Pre-treatment (matched) good vs bad ===')
pre = S[(S.timepoint=='pre') & S.matched]
for col in ['n_mutation_sites','n_candidate_peptides','n_binders_500nM','n_strong_binders_50nM',
            'n_sites_with_binder','n_sites_with_strong','neoantigens_per_site','PCN_score']:
    g = pre[pre.response_bin=='good'][col].dropna()
    b = pre[pre.response_bin=='bad'][col].dropna()
    if len(g)>=3 and len(b)>=3:
        u = stats.mannwhitneyu(g,b)
        print(f'  {col}: good med={g.median():.1f} (n={len(g)}) vs bad med={b.median():.1f} (n={len(b)})  p={u.pvalue:.3f}')

# ------ Paired pre-post delta ------
print('\n=== Pre→Post Δ good vs bad ===')
pair = S[S.timepoint.isin(['pre','post']) & S.matched]
subs = [s for s,g in pair.groupby('subject_id') if set(g.timepoint)>={'pre','post'}]
rows_d=[]
for s in subs:
    sub = pair[pair.subject_id==s]
    pre = sub[sub.timepoint=='pre'].iloc[0]
    post = sub[sub.timepoint=='post'].iloc[0]
    rows_d.append({'subject_id':s,'response':sub.response_bin.iloc[0],
                   'delta_binders':post.n_binders_500nM - pre.n_binders_500nM,
                   'delta_sites':post.n_sites_with_binder - pre.n_sites_with_binder,
                   'delta_strong':post.n_strong_binders_50nM - pre.n_strong_binders_50nM,
                   'delta_PCN':post.PCN_score - pre.PCN_score})
D = pd.DataFrame(rows_d)
D.to_csv(OUT/'neoantigen_paired_delta.tsv', sep='\t', index=False)
print(D.to_string(index=False))
print()
for col in ['delta_binders','delta_sites','delta_strong','delta_PCN']:
    g = D[D.response=='good'][col].dropna()
    b = D[D.response=='bad'][col].dropna()
    if len(g)>=3 and len(b)>=3:
        u = stats.mannwhitneyu(g,b)
        print(f'  {col}: good med={g.median():.0f} (n={len(g)}) vs bad med={b.median():.0f} (n={len(b)})  p={u.pvalue:.3f}')

# ------ Final summary with all sample info ------
print('\n=== Sample distribution (neoantigen count ranges) ===')
print(S.groupby(['timepoint','response_bin'])['n_sites_with_binder'].describe()[['count','mean','50%','min','max']])
