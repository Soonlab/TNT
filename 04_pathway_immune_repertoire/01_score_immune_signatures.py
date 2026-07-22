"""
Score purified immune signatures (CD8_cytotoxic, Tcell_infiltration, Bcell_infiltration,
TLS_Cabrita) on TNT RNA-seq via ssGSEA, then collapse to per-subject pre-treatment values
that can be merged into integrated_subject_master.tsv.

Why: external validation v3 (script 32) showed CD8_cytotoxic axis is the only signal that
robustly reproduces (Z=+2.74, P=0.006, 9 cohorts N=721) while the legacy CD8_proliferation
in the master table is actually contaminated with cell-cycle genes (MKI67/TOP2A/MCM/CCN/CDK).
We add the *purified* immune signatures so the discovery LASSO can also discover them.

Output:
  scores_ssgsea_immune.tsv             — sample-level scores (rows = TNT_RNA_*, cols = sigs)
  signature_gene_coverage.tsv          — how many genes from each sig were measured in TPM
  pre_subject_immune_scores.tsv        — per-subject pre-treatment scores ready to merge
"""
import os, sys
import numpy as np
import pandas as pd
import gseapy as gp

OUT = '/data/data/TNT/analysis/260418_add'
TPM = '/data/data/TNT/analysis/06_rna_immune/tpm_symbol.tsv'
META = '/data/data/TNT/analysis/00_cohort/rna_inventory.tsv'
os.makedirs(OUT, exist_ok=True)

# Purified signatures (matching scripts/32_external_validation_v3_CD8axis.py)
SIGS = {
    'CD8_cytotoxic': ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1','IFNG','NKG7','GNLY',
                     'CXCL9','CXCL10','CXCL11','TBX21','EOMES','KLRK1','KLRD1'],
    'Tcell_infiltration': ['CD3D','CD3E','CD3G','CD2','CD4','CD8A','CD8B','LCK','ZAP70','ITK'],
    'Bcell_infiltration': ['CD19','CD20','MS4A1','CD79A','CD79B','CD22','TCL1A','FCRL5','BLK','FCER2'],
    # Cabrita 2020 12-chemokine TLS signature — robust across solid tumours
    'TLS_Cabrita':       ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3'],
}

print('Loading TPM matrix...')
tpm = pd.read_csv(TPM, sep='\t', index_col=0)
print(f'  TPM shape = {tpm.shape}  ({tpm.shape[0]} genes x {tpm.shape[1]} samples)')
log_tpm = np.log2(tpm + 1)

# Coverage report
cov_rows = []
for sig, genes in SIGS.items():
    found = [g for g in genes if g in log_tpm.index]
    miss  = [g for g in genes if g not in log_tpm.index]
    cov_rows.append({'signature': sig, 'n_total': len(genes), 'n_found': len(found),
                     'pct': round(100*len(found)/len(genes),1),
                     'missing_genes': ','.join(miss) if miss else ''})
cov = pd.DataFrame(cov_rows)
cov.to_csv(f'{OUT}/signature_gene_coverage.tsv', sep='\t', index=False)
print('\n=== signature gene coverage ===')
print(cov.to_string(index=False))

print('\nRunning ssGSEA (rank-norm, gseapy)...')
ss = gp.ssgsea(data=log_tpm, gene_sets=SIGS, sample_norm_method='rank',
               min_size=3, max_size=500, outdir=None, no_plot=True, processes=4)
score = ss.res2d.pivot(index='Term', columns='Name', values='NES').T
score.index.name = 'sample_id'
score = score.astype(float)
print(f'  score matrix: {score.shape}')
score.to_csv(f'{OUT}/scores_ssgsea_immune.tsv', sep='\t')

# Collapse to per-subject pre-treatment
meta = pd.read_csv(META, sep='\t')
print(f'\nRNA inventory: {len(meta)} rows; timepoints = {meta.timepoint.value_counts().to_dict()}')
pre = meta[meta.timepoint == 'pre'][['sample_id','subject_id']].copy()
print(f'  pre samples in inventory: {len(pre)}')

merged = pre.merge(score.reset_index(), on='sample_id', how='inner')
print(f'  pre samples with ssGSEA scores: {len(merged)}')
sub = merged.set_index('subject_id')[list(SIGS.keys())]
# z-score within cohort so the new features are on comparable scale to existing ssGSEA cols
sub_z = (sub - sub.mean(axis=0)) / sub.std(axis=0)
sub_z.columns = [c for c in sub.columns]
sub_z.reset_index().to_csv(f'{OUT}/pre_subject_immune_scores.tsv', sep='\t', index=False)
print(f'\nWrote per-subject pre-treatment immune scores: {OUT}/pre_subject_immune_scores.tsv')
print(sub_z.head().to_string())
print(f'\nDone. {len(sub_z)} subjects scored on {len(SIGS)} purified immune signatures.')
