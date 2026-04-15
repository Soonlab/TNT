"""
Build unified somatic variant master table, compute TMB, driver mutation oncoprint.
"""
import pandas as pd, numpy as np, os, re
from pathlib import Path

TSV = Path('/mnt/sda1/data/TNT/analysis/02_wes_mutect2/variant_tables')
OUT = Path('/mnt/sda1/data/TNT/analysis/02_wes_tmb_msi'); OUT.mkdir(parents=True, exist_ok=True)
DRV = Path('/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal'); DRV.mkdir(parents=True, exist_ok=True)
INV = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')

# Standardize columns
cols = ['CHROM','POS','REF','ALT','EFFECT','IMPACT','GENE','HGVS_P','FEATUREID','FILTER','AD','DP','AF']
rename = {'ANN[0].EFFECT':'EFFECT','ANN[0].IMPACT':'IMPACT','ANN[0].GENE':'GENE',
          'ANN[0].HGVS_P':'HGVS_P','ANN[0].FEATUREID':'FEATUREID',
          'GEN[0].AD':'AD','GEN[0].DP':'DP','GEN[0].AF':'AF'}

NONSYN = {'missense_variant','stop_gained','stop_lost','start_lost','frameshift_variant',
          'inframe_insertion','inframe_deletion','splice_acceptor_variant','splice_donor_variant',
          'protein_altering_variant','initiator_codon_variant'}
ALL_CODING = NONSYN | {'synonymous_variant','stop_retained_variant','splice_region_variant'}

# SureSelect V5 target ~50Mb
EXOME_MB = 50.0

UNMATCHED_SUBJ = {13,15,16,17,18,19,33}

all_vars = []
for f in sorted(TSV.glob('*.tsv')):
    sid = f.stem
    df = pd.read_csv(f, sep='\t')
    df = df.rename(columns=rename)
    df['sample_id'] = sid
    all_vars.append(df)

M = pd.concat(all_vars, ignore_index=True)
# parse VAF / DP
def first_float(x):
    try: return float(str(x).split(',')[0])
    except: return np.nan
M['AF_f'] = M['AF'].apply(first_float)
M['DP_f'] = M['DP'].apply(first_float)

# Effect classification (first effect)
M['EFFECT_primary'] = M['EFFECT'].astype(str).str.split('&').str[0]
M['is_nonsyn'] = M['EFFECT_primary'].isin(NONSYN)
M['is_coding'] = M['EFFECT_primary'].isin(ALL_CODING)
M['is_high_impact'] = M['IMPACT'].isin(['HIGH','MODERATE'])

# Merge with meta
M = M.merge(INV[['sample_id','subject_id','timepoint','response_bin','response_num']], on='sample_id', how='left')
M['matched'] = ~M['subject_id'].isin(UNMATCHED_SUBJ)
M.to_csv(OUT/'variant_master.tsv.gz', sep='\t', index=False, compression='gzip')
print('Variant master:', M.shape, '→', OUT/'variant_master.tsv.gz')

# TMB per sample
tmb = M.groupby('sample_id').agg(
    n_total=('CHROM','size'),
    n_coding=('is_coding','sum'),
    n_nonsyn=('is_nonsyn','sum'),
    n_high_impact=('is_high_impact','sum'),
).reset_index()
tmb['TMB_nonsyn_per_Mb'] = tmb['n_nonsyn']/EXOME_MB
tmb['TMB_coding_per_Mb'] = tmb['n_coding']/EXOME_MB
tmb['TMB_all_per_Mb'] = tmb['n_total']/EXOME_MB
tmb = tmb.merge(INV[['sample_id','subject_id','timepoint','response_bin','response_num']], on='sample_id')
tmb['matched'] = ~tmb['subject_id'].isin(UNMATCHED_SUBJ)
# Hypermutator threshold commonly 10/Mb (TCGA), MSI-high often >20/Mb
tmb['hypermutator_10'] = tmb['TMB_nonsyn_per_Mb']>=10
tmb['msi_like_20'] = tmb['TMB_nonsyn_per_Mb']>=20
tmb.to_csv(OUT/'tmb_per_sample.tsv', sep='\t', index=False)

# Response association (matched + pre only, to remove unmatched bias)
from scipy import stats
print('\n=== TMB good vs bad by timepoint (matched samples only) ===')
for tp in ['pre','post']:
    sub = tmb[(tmb.timepoint==tp) & tmb.matched]
    g = sub[sub.response_bin=='good']['TMB_nonsyn_per_Mb']
    b = sub[sub.response_bin=='bad']['TMB_nonsyn_per_Mb']
    if len(g)>=2 and len(b)>=2:
        u = stats.mannwhitneyu(g,b)
        print(f'  {tp}: good n={len(g)} median={g.median():.2f} IQR=({g.quantile(.25):.1f}-{g.quantile(.75):.1f}) | bad n={len(b)} median={b.median():.2f} IQR=({b.quantile(.25):.1f}-{b.quantile(.75):.1f}) | p={u.pvalue:.3f}')

# ----- Driver mutation panel -----
CRC_DRIVERS = ['TP53','APC','KRAS','BRAF','SMAD4','PIK3CA','NRAS','FBXW7','TCF7L2','RNF43','AMER1',
               'SMAD2','CTNNB1','FAM123B','ACVR2A','SOX9','ARID1A','BCOR','ATM','POLE','POLD1',
               'MLH1','MSH2','MSH6','PMS2','CDKN2A','STK11','PTEN','ERBB2','ERBB3','MTOR',
               'KMT2D','CREBBP','EP300','NOTCH1','CDC27','BMPR2','GNAS','HRAS','MYC','CCND1']

drv = M[M['GENE'].isin(CRC_DRIVERS) & M['is_nonsyn']].copy()
drv_table = drv.groupby(['sample_id','GENE']).agg(n=('CHROM','size'), effects=('EFFECT_primary', lambda x: ';'.join(sorted(set(x))))).reset_index()
drv_table.to_csv(DRV/'driver_mutations.tsv', sep='\t', index=False)

# Wide oncoprint table
wide = drv_table.pivot_table(index='GENE', columns='sample_id', values='n', fill_value=0, aggfunc='sum')
# Sort genes by total carriers
wide['total_samples_mutated'] = (wide>0).sum(axis=1)
wide = wide.sort_values('total_samples_mutated', ascending=False)
wide.to_csv(DRV/'driver_oncoprint_matrix.tsv', sep='\t')
print('\n=== Top CRC drivers (sample count) ===')
print(wide[['total_samples_mutated']].head(15))

# Response association per driver (pre matched only)
pre_matched = set(tmb[(tmb.timepoint=='pre') & tmb.matched].sample_id)
pre_drv = drv_table[drv_table.sample_id.isin(pre_matched)]
pre_meta = tmb[(tmb.timepoint=='pre') & tmb.matched][['sample_id','response_bin']]

rows=[]
for gene in wide.index[:25]:
    carriers = set(pre_drv[pre_drv.GENE==gene].sample_id)
    good_total = (pre_meta.response_bin=='good').sum()
    bad_total = (pre_meta.response_bin=='bad').sum()
    good_mut = pre_meta[(pre_meta.sample_id.isin(carriers)) & (pre_meta.response_bin=='good')].shape[0]
    bad_mut  = pre_meta[(pre_meta.sample_id.isin(carriers)) & (pre_meta.response_bin=='bad')].shape[0]
    if good_mut+bad_mut<1: continue
    table = [[good_mut, good_total-good_mut],[bad_mut, bad_total-bad_mut]]
    odds, p = stats.fisher_exact(table)
    rows.append((gene, good_mut, good_total, bad_mut, bad_total, odds, p))
driv_stat = pd.DataFrame(rows, columns=['GENE','good_mut','good_n','bad_mut','bad_n','OR','p']).sort_values('p')
driv_stat.to_csv(DRV/'driver_response_fisher.tsv', sep='\t', index=False)
print('\n=== Driver mutation response association (pre, matched) ===')
print(driv_stat.head(15).to_string(index=False))
print('\nAll outputs in 02_wes_tmb_msi/ and 04_wes_cnv_clonal/')
