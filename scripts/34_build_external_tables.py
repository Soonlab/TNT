"""Build Table 3 (main) and Table S7 (supplementary) for external validation v3."""
import pandas as pd
from pathlib import Path

OUT  = Path('/mnt/sda1/data/TNT/analysis/11_external_validation')
TBL  = Path('/mnt/sda1/data/TNT/analysis/tables'); TBL.mkdir(parents=True, exist_ok=True)

stats = pd.read_csv(OUT/'v3_signature_response_stats.tsv', sep='\t')
meta  = pd.read_csv(OUT/'v3_meta_overall.tsv', sep='\t')
summ  = pd.read_csv(OUT/'v3_cohort_summary.tsv', sep='\t')

SIG_ORDER = ['CD8_cytotoxic','Tcell_infiltration','Bcell_infiltration',
             'Tumor_cellcycle','DSB_HDR_repair','E2F_MYC_cellcycle','EMT']

# ---- Table 3 (main): meta summary with concordance ----
rows = []
for _,m in meta.iterrows():
    sig = m.signature
    sub = stats[stats.signature==sig]
    n_pos = int((sub.delta > 0).sum())
    n_neg = int((sub.delta < 0).sum())
    rows.append({
        'Signature': sig,
        'Genes in signature (example)': {
            'CD8_cytotoxic':'CD8A/B, GZMA/B/H/K, PRF1, IFNG, NKG7, GNLY, CXCL9/10/11, TBX21, EOMES',
            'Tcell_infiltration':'CD3D/E/G, CD2, CD4, CD8A/B, LCK, ZAP70, ITK',
            'Bcell_infiltration':'CD19, MS4A1, CD79A/B, CD22, TCL1A, FCRL5, BLK, FCER2',
            'Tumor_cellcycle':'MKI67, TOP2A, MCM2/5, CCNB1/B2, CDK1, PCNA, CENPF, BUB1, PLK1, AURKA/B',
            'DSB_HDR_repair':'BRCA1/2, RAD51, PALB2, ATM/ATR, CHEK1/2, FANCA/D2/I/L, MRE11, RAD50',
            'E2F_MYC_cellcycle':'E2F1/2/3, MYC, MCM3/4/6/7, CCNE1/E2, CDK2/4/6, CDC20/25A',
            'EMT':'VIM, CDH2, SNAI1/2, TWIST1/2, ZEB1/2, FN1, MMP2/3/9, COL1A1/1A2/3A1, FAP, ACTA2'
        }.get(sig,''),
        'n_cohorts': m.n_cohorts,
        'Concordant (Δ>0)': f'{n_pos} / {m.n_cohorts}',
        'Stouffer Z': f'{m.Z:+.2f}',
        'p_meta': f'{m.p_meta:.3f}',
        'Significance': '** (P<0.01)' if m.p_meta<0.01 else ('* (P<0.05)' if m.p_meta<0.05 else ('trend (P<0.15)' if m.p_meta<0.15 else 'ns'))
    })
tab3 = pd.DataFrame(rows).set_index('Signature').loc[SIG_ORDER].reset_index()
tab3.to_csv(TBL/'Table3_external_meta_summary.tsv', sep='\t', index=False)

# ---- Table S7 (supp): per-cohort per-signature full detail ----
# Wide layout: rows=cohort, cols = signature-level (delta, p) + cohort metadata
pivot_delta = stats.pivot(index='gse', columns='signature', values='delta').round(3)
pivot_p     = stats.pivot(index='gse', columns='signature', values='pvalue').round(4)
pivot_delta.columns = [f'{c}__delta' for c in pivot_delta.columns]
pivot_p.columns     = [f'{c}__p'     for c in pivot_p.columns]
cohort_info = summ.set_index('gse')[['n_samples','n_good','n_bad','resp_col','scale','regimen','n_probes']]
full = cohort_info.join(pivot_delta).join(pivot_p)
# reorder columns: info, then per-signature (delta,p) pairs
sig_cols = []
for s in SIG_ORDER:
    if f'{s}__delta' in full.columns:
        sig_cols += [f'{s}__delta', f'{s}__p']
full = full[['n_samples','n_good','n_bad','regimen','resp_col','scale'] + sig_cols]
# Add response-scale provenance column
scale_notes = {
    'GSE150082':'Explicit Good/Poor (author-assigned)',
    'GSE35452' :'Responder / Non-Responder (author)',
    'GSE119409':'Sensitive / Resistant (Mandard TRG1-2 vs 3-5 per source)',
    'GSE45404' :'Responder / Non-Responder (class field)',
    'GSE94104' :'Rödel TRG 3 = good, 1-2 = bad',
    'GSE56699' :'RCRG 3-class; complete+partial = good, poor = bad',
    'GSE46862' :'4-class TO/MO/MI/NT; TO+MO = good, MI+NT = bad',
    'GSE133057':'AJCC TRG 0-1 = good, 2-3 = bad',
    'GSE87211' :'Cancer recurrence surrogate; 0=no recurrence=good, 1=recurrence=bad',
}
full['response_scale_note'] = full.index.map(scale_notes)
full = full.loc[[g for g in list(scale_notes) if g in full.index]]
full.to_csv(TBL/'TableS7_external_percohort_signatures.tsv', sep='\t')

# ---- Print both ----
print('=== Table 3 (main): external meta-analysis summary ===')
print(tab3.to_string(index=False))
print()
print('=== Table S7 (supp): per-cohort, per-signature detail ===')
# Print compact view
display = full[['n_samples','n_good','n_bad','regimen','response_scale_note',
                'CD8_cytotoxic__delta','CD8_cytotoxic__p',
                'Tcell_infiltration__delta','Bcell_infiltration__delta',
                'Tumor_cellcycle__delta','DSB_HDR_repair__delta','EMT__delta']]
print(display.to_string())

print('\nWritten:')
print(' ', TBL/'Table3_external_meta_summary.tsv')
print(' ', TBL/'TableS7_external_percohort_signatures.tsv')
