"""
ssGSEA per-sample pathway scores via gseapy.ssgsea
Hallmark + curated Reactome immune/cell cycle/ECM sets
"""
import pandas as pd, numpy as np, os
import gseapy as gp
from scipy import stats
from statsmodels.stats.multitest import multipletests

OUT = '/mnt/sda1/data/TNT/analysis/08_rna_pathway'
os.makedirs(OUT, exist_ok=True)

tpm = pd.read_csv('/mnt/sda1/data/TNT/analysis/06_rna_immune/tpm_symbol.tsv', sep='\t', index_col=0)
log_tpm = np.log2(tpm+1)
print('log_tpm:', log_tpm.shape)

# Hallmark gene sets
hm_file = gp.get_library(name='MSigDB_Hallmark_2020', organism='Human')
print(f'Hallmark pathways: {len(hm_file)}')

# Key Reactome
reactome_all = gp.get_library(name='Reactome_2022', organism='Human')
keep_re = {k for k in reactome_all if any(s in k for s in
    ['Cell Cycle','DNA Repair','Homology','Double-Strand Break','Interferon',
     'Antigen Processing','MHC','Extracellular Matrix','S Phase','M Phase',
     'Homologous Recombination','Mismatch Repair','Base Excision','Nucleotide Excision'])}
print(f'Reactome subset: {len(keep_re)}')
reactome_sub = {k: reactome_all[k] for k in keep_re}

all_gs = {**hm_file, **reactome_sub}
print(f'Total pathways: {len(all_gs)}')

# ssGSEA
ss = gp.ssgsea(data=log_tpm, gene_sets=all_gs, sample_norm_method='rank',
               min_size=10, max_size=500, outdir=None, no_plot=True, processes=8)
score = ss.res2d.pivot(index='Term', columns='Name', values='NES')
score = score.T  # samples x pathways
score.index.name='sample_id'
score.to_csv(f'{OUT}/ssgsea_scores.tsv', sep='\t')
print('score matrix:', score.shape)

# Response association
meta = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv', sep='\t')
m = score.reset_index().merge(meta[['sample_id','timepoint','response_bin']], on='sample_id')
rows=[]
for tp in ['pre','post']:
    sub = m[m.timepoint==tp]
    if len(sub)<8: continue
    for p in score.columns:
        g = sub[sub.response_bin=='good'][p].dropna()
        b = sub[sub.response_bin=='bad'][p].dropna()
        if len(g)<3 or len(b)<3: continue
        u = stats.mannwhitneyu(g,b)
        rows.append((tp, p, len(g), len(b), g.mean(), b.mean(), g.mean()-b.mean(), u.pvalue))
res = pd.DataFrame(rows, columns=['timepoint','pathway','n_good','n_bad','mean_good','mean_bad','delta','pvalue'])
res['qvalue']=np.nan
for tp, idx in res.groupby('timepoint').groups.items():
    _,q,_,_=multipletests(res.loc[idx,'pvalue'],method='fdr_bh')
    res.loc[idx,'qvalue']=q
res=res.sort_values(['timepoint','pvalue'])
res.to_csv(f'{OUT}/ssgsea_response_stats.tsv', sep='\t', index=False)

print('\n=== pre top 15 ===')
print(res[res.timepoint=='pre'].head(15).to_string(index=False))
print('\n=== post top 10 ===')
print(res[res.timepoint=='post'].head(10).to_string(index=False))
