"""
Build TPM matrix (symbol-level, collapsed), compute immune/pathway signatures,
and test response association for TNT RNA-seq cohort.
Outputs: 06_rna_immune/{tpm_symbol.tsv, logtpm_symbol.tsv, signature_scores.tsv, sig_response_stats.tsv}
"""
import pandas as pd, numpy as np, os
from scipy import stats

RNA_X='/mnt/sda1/data/TNT/TNT_RNAseq/result_RNAseq/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.xlsx'
META='/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv'
OUT='/mnt/sda1/data/TNT/analysis/06_rna_immune'
os.makedirs(OUT, exist_ok=True)

print('Loading expression xlsx...')
expr = pd.read_excel(RNA_X, sheet_name=0)
tpm_cols = [c for c in expr.columns if c.endswith('_TPM')]
cnt_cols = [c for c in expr.columns if c.endswith('_Read_Count')]
print(f' {len(tpm_cols)} TPM cols, {len(cnt_cols)} count cols')

# Symbol collapse (max TPM per symbol; drop empty symbols)
sym = expr[['Gene_Symbol'] + tpm_cols].dropna(subset=['Gene_Symbol'])
sym.columns = ['symbol'] + [c.replace('_TPM','') for c in tpm_cols]
sym = sym.groupby('symbol').max()
print(f' TPM matrix: {sym.shape[0]} symbols x {sym.shape[1]} samples')
sym.to_csv(f'{OUT}/tpm_symbol.tsv', sep='\t')

log_tpm = np.log2(sym + 1)
log_tpm.to_csv(f'{OUT}/logtpm_symbol.tsv', sep='\t')

# -------- Signature definitions --------
SIGS = {
 # CD8 function (Tirosh/Bindea + canonical)
 'CD8_activation': ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1','IFNG','NKG7','GNLY','CD3D','CD3E'],
 'CD8_proliferation': ['MKI67','TOP2A','STMN1','TYMS','TUBB','UBE2C','BIRC5','CCNB1','CCNB2','CDK1','MCM2','MCM5','PCNA','CENPF','KIF20A'],
 'CD8_exhaustion': ['PDCD1','LAG3','HAVCR2','TIGIT','TOX','CTLA4','ENTPD1','EOMES','BATF'],
 'Cytolytic_activity': ['GZMA','PRF1'],
 # Antigen presentation
 'Antigen_presentation': ['B2M','HLA-A','HLA-B','HLA-C','TAP1','TAP2','TAPBP','CALR','CANX','ERAP1','ERAP2','PSMB8','PSMB9','PSME1','PSME2','NLRC5'],
 'MHC_II': ['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','HLA-DMA','HLA-DMB','CIITA'],
 'NLRC5_HLA_IFNG': ['NLRC5','HLA-A','HLA-B','HLA-C','B2M','IFNG','IRF1','STAT1','STAT2'],
 'IFNg_Ayers_18': ['IFNG','STAT1','CCR5','CXCL9','CXCL10','CXCL11','IDO1','PRF1','GZMA','GZMB','MHC-class-II'.replace('MHC-class-II','HLA-DRA'),'CD8A','HLA-E','NKG7','HLA-DQA1','CD274','LAG3','TIGIT'],
 # TLS (Cabrita 12-gene)
 'TLS_Cabrita': ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3','CD79B','MS4A1','CCL18','CXCL8','PTGDS'],
 # Checkpoints
 'Checkpoint_inhibitory': ['CD274','PDCD1LG2','PDCD1','CTLA4','LAG3','HAVCR2','TIGIT','BTLA','VISTA','IDO1'],
 # Resistance / immune-exclusion
 'TGFb_Mariathasan': ['ACTA2','ACTG2','ADAM12','ADAM19','CNN1','COL1A1','COL1A2','COL3A1','COL4A1','COL5A1','COL5A2','COL5A3','COL6A1','COL6A2','COL6A3','COL7A1','COL8A1','COL8A2','COL10A1','COL11A1','COL12A1','COL14A1','COL15A1','COL16A1','FAP','LRRC15','MMP2','MMP11','MMP14','POSTN','SPARC','TGFB1','TGFB2','TGFB3','TGFBR1','TGFBR2','THBS2','VCAN'],
 'EMT_Mak': ['VIM','CDH2','FOXC2','SNAI1','SNAI2','TWIST1','FN1','ITGB6','MMP2','MMP3','MMP9','SOX10','GSC','ZEB1','ZEB2','TWIST2'],
 'Epithelial': ['CDH1','DSP','OCLN','CLDN3','CLDN4','CLDN7','KRT8','KRT18','KRT19','EPCAM'],
 'Hypoxia_Buffa': ['ACOT7','ADM','ALDOA','ANKRD37','ANLN','BNIP3','C20orf20','CA9','CDKN3','COL4A5','DCBLD1','DDIT4','DTYMK','ENO1','FAM162A','GAPDH','HIG2','HK2','KCTD11','KIF20A','LDHA','LRRC42','MAFF','MCTS1','MIF','MRPL13','MRPL15','MRPS17','NDRG1','NP','P4HA1','P4HA2','PFKP','PGAM1','PGK1','PSMA7','PSRC1','PTP4A3','SEC61G','SHCBP1','SLC16A1','SLC25A32','SLC2A1','TPI1','TUBA1B','TUBA1C','UTP11L','VEGFA','YKT6'],
 'Stemness_mRNAsi_proxy': ['SOX2','NANOG','POU5F1','KLF4','MYC','LIN28A','LIN28B','PROM1','ALDH1A1','CD44'],
 # CAF / stromal
 'CAF_iCAF': ['HAS1','HAS2','CXCL12','CXCL1','CXCL2','IL6','IL8','LIF','PDGFRA'],
 'CAF_myCAF': ['ACTA2','TAGLN','MYL9','MYLK','TPM1','TPM2','POSTN','CTGF'],
 # Expanded immune cells / mini-deconvolution approximation
 'Mac_M1': ['NOS2','IL12A','IL12B','TNF','IL1A','IL1B','CXCL9','CXCL10','CXCL11','IRF5'],
 'Mac_M2': ['CD163','MRC1','MSR1','CD68','IL10','TGFB1','ARG1','CCL22','CCL17'],
 'Treg': ['FOXP3','IL2RA','CTLA4','IKZF2','CCR8','ENTPD1'],
 'NK_cell': ['NKG7','KLRD1','KLRF1','NCAM1','FCGR3A','PRF1','GNLY','KIR2DL3'],
 'B_cell': ['CD19','CD20','MS4A1','CD79A','CD79B','BANK1','FCRL2'],
}

# z-score per-signature (mean z across member genes present)
def score_sig(log_tpm, genes):
    present = [g for g in genes if g in log_tpm.index]
    if len(present) < 2:
        return None, present
    sub = log_tpm.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0,np.nan), axis=0)
    return z.mean(axis=0), present

scores={}
sig_info=[]
for sname, genes in SIGS.items():
    s, present = score_sig(log_tpm, genes)
    if s is None:
        print(f'SKIP {sname} (not enough genes)'); continue
    scores[sname] = s
    sig_info.append((sname, len(genes), len(present), ','.join([g for g in genes if g not in present])))

score_df = pd.DataFrame(scores)
score_df.index.name='sample_id'
score_df.to_csv(f'{OUT}/signature_scores.tsv', sep='\t')
pd.DataFrame(sig_info, columns=['signature','n_genes_defined','n_genes_found','missing']).to_csv(
    f'{OUT}/signature_gene_coverage.tsv', sep='\t', index=False)

# ---- Response association ----
meta = pd.read_csv(META, sep='\t')
merged = score_df.reset_index().merge(meta[['sample_id','subject_id','timepoint','response_bin','response_num','prepost_set','cT','sex','age']], on='sample_id')

results=[]
for tp in ['pre','post','normal','all']:
    sub = merged if tp=='all' else merged[merged.timepoint==tp]
    if len(sub) < 5: continue
    g = sub[sub.response_bin=='good']; b = sub[sub.response_bin=='bad']
    for sname in score_df.columns:
        if len(g)<2 or len(b)<2: continue
        u = stats.mannwhitneyu(g[sname], b[sname], alternative='two-sided')
        results.append((tp, sname, len(g), len(b), g[sname].mean(), b[sname].mean(),
                        g[sname].mean()-b[sname].mean(), u.statistic, u.pvalue))
res = pd.DataFrame(results, columns=['timepoint','signature','n_good','n_bad','mean_good','mean_bad','delta_good_minus_bad','U','pvalue'])
# BH FDR within each timepoint
from statsmodels.stats.multitest import multipletests
res['qvalue']=np.nan
for tp, idx in res.groupby('timepoint').groups.items():
    _,q,_,_=multipletests(res.loc[idx,'pvalue'],method='fdr_bh')
    res.loc[idx,'qvalue']=q
res=res.sort_values(['timepoint','pvalue'])
res.to_csv(f'{OUT}/sig_response_stats.tsv', sep='\t', index=False)

print('\n=== Pre-treatment: top signatures (good vs bad) ===')
print(res[res.timepoint=='pre'].head(15).to_string(index=False))
print('\nFiles written in', OUT)
