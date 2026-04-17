"""
Score CD8_cytotoxic / Tcell_infiltration / Bcell_infiltration / TLS_Cabrita on the FULL
RNA-seq matrix (pre + post + normal), then build per-subject pre-post Δ for the paired
subset.

Why: existing 09_integration/paired_delta/paired_feature_long.tsv has 8 features
(missense, SBS5, neo_binders, neo_sites, MHC_II, Treg, CD8_exhaustion, IGH_n, TRB_shannon)
but does NOT have the externally-validated CD8_cytotoxic axis. We need the new immune sig
deltas too to test convergence (baseline tumor-intrinsic vs Δ CD8_cyt).

Output:
  ssgsea_immune_all_samples.tsv     — sample-level scores (rows = TNT_RNA_*)
  paired_immune_delta_per_subject.tsv  — per paired subject, Δ for 4 new sigs
"""
import os, numpy as np, pandas as pd, gseapy as gp

OUT = '/data/data/TNT/analysis/260418_add'
TPM = '/data/data/TNT/analysis/06_rna_immune/tpm_symbol.tsv'
META = '/data/data/TNT/analysis/00_cohort/rna_inventory.tsv'

SIGS = {
    'CD8_cytotoxic':       ['CD8A','CD8B','GZMA','GZMB','GZMH','GZMK','PRF1','IFNG','NKG7','GNLY',
                            'CXCL9','CXCL10','CXCL11','TBX21','EOMES','KLRK1','KLRD1'],
    'Tcell_infiltration':  ['CD3D','CD3E','CD3G','CD2','CD4','CD8A','CD8B','LCK','ZAP70','ITK'],
    'Bcell_infiltration':  ['CD19','CD20','MS4A1','CD79A','CD79B','CD22','TCL1A','FCRL5','BLK','FCER2'],
    'TLS_Cabrita':         ['CCL19','CCL21','CXCL13','CCR7','CXCR5','SELL','LAMP3'],
}

print('Loading TPM...')
tpm = pd.read_csv(TPM, sep='\t', index_col=0)
log_tpm = np.log2(tpm + 1)
print(f'  TPM = {tpm.shape}')

print('Running ssGSEA on all 56 samples...')
ss = gp.ssgsea(data=log_tpm, gene_sets=SIGS, sample_norm_method='rank',
               min_size=3, max_size=500, outdir=None, no_plot=True, processes=4)
score = ss.res2d.pivot(index='Term', columns='Name', values='NES').T.astype(float)
score.index.name = 'sample_id'
score.to_csv(f'{OUT}/ssgsea_immune_all_samples.tsv', sep='\t')
print(f'  scored: {score.shape}')

# Z-score across all samples for comparability with legacy paired_feature_long table
score_z = (score - score.mean(axis=0)) / score.std(axis=0)

# Build per-subject paired Δ
meta = pd.read_csv(META, sep='\t')
m = score_z.reset_index().merge(meta[['sample_id','subject_id','timepoint','response_bin']],
                                on='sample_id')
print(f'\nTimepoint distribution after merge: {m.timepoint.value_counts().to_dict()}')

# Subjects with both pre AND post
pre = m[m.timepoint=='pre'][['subject_id','response_bin'] + list(SIGS.keys())].copy()
post = m[m.timepoint=='post'][['subject_id'] + list(SIGS.keys())].copy()
pre.columns = ['subject_id','response_bin'] + [f'{c}_pre' for c in SIGS]
post.columns = ['subject_id'] + [f'{c}_post' for c in SIGS]
paired = pre.merge(post, on='subject_id', how='inner')
print(f'paired (pre + post both present): {len(paired)} subjects')
print(paired[['subject_id','response_bin']].sort_values('subject_id').to_string(index=False))

for s in SIGS:
    paired[f'{s}_delta'] = paired[f'{s}_post'] - paired[f'{s}_pre']

paired.to_csv(f'{OUT}/paired_immune_delta_per_subject.tsv', sep='\t', index=False)
print(f'\nWrote {OUT}/paired_immune_delta_per_subject.tsv')
print('\n=== summary by response ===')
for sig in SIGS:
    g = paired[paired.response_bin=='good'][f'{sig}_delta'].values
    b = paired[paired.response_bin=='bad'][f'{sig}_delta'].values
    print(f'  {sig:22s}  good Δmed={np.median(g):+.3f} (n={len(g)})   bad Δmed={np.median(b):+.3f} (n={len(b)})')
