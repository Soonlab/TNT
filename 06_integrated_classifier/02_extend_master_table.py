"""
Merge purified immune signatures into integrated_subject_master.tsv.

Output:
  integrated_subject_master_v2.tsv    — original 37 cols + 4 new immune sigs (41 cols)
  master_v2_diff.txt                  — diff summary (which subjects gained which features)
"""
import pandas as pd
import numpy as np

OUT  = '/data/data/TNT/analysis/260418_add'
SRC  = '/data/data/TNT/analysis/tables/integrated_subject_master.tsv'
ADD  = f'{OUT}/pre_subject_immune_scores.tsv'

m  = pd.read_csv(SRC, sep='\t')
ad = pd.read_csv(ADD, sep='\t')
ad['subject_id'] = ad['subject_id'].astype(str)
m['subject_id']  = m['subject_id'].astype(str)

# Drop new TLS_Cabrita: master already has a TLS_Cabrita from the legacy 22-signature
# pipeline (06_rna_immune/signature_scores.tsv). Avoid name collision; keep the existing
# one so we are not silently changing what TLS_Cabrita means in the table.
if 'TLS_Cabrita' in ad.columns and 'TLS_Cabrita' in m.columns:
    print('  TLS_Cabrita already in master -> dropping new one to avoid collision')
    ad = ad.drop(columns=['TLS_Cabrita'])

print(f'master: {m.shape}, new immune scores: {ad.shape}')
print(f'subjects in master: {sorted(m.subject_id.tolist(), key=lambda x: int(x))}')
print(f'subjects with new RNA scores: {sorted(ad.subject_id.tolist(), key=lambda x: int(x))}')
miss = set(m.subject_id) - set(ad.subject_id)
print(f'subjects missing pre-RNA (will be NaN for new features): {sorted(miss, key=lambda x: int(x))}')

merged = m.merge(ad, on='subject_id', how='left')
print(f'merged: {merged.shape}')

# Sanity: which existing columns vs new columns
new_cols = [c for c in ad.columns if c != 'subject_id']
print(f'new columns added: {new_cols}')

# Save
out_path = f'{OUT}/integrated_subject_master_v2.tsv'
merged.to_csv(out_path, sep='\t', index=False)
print(f'\nWrote {out_path}')

# Diff report
with open(f'{OUT}/master_v2_diff.txt','w') as f:
    f.write(f'integrated_subject_master.tsv     -> v2\n')
    f.write(f'rows:    {m.shape[0]} -> {merged.shape[0]}\n')
    f.write(f'columns: {m.shape[1]} -> {merged.shape[1]} (+{merged.shape[1]-m.shape[1]})\n\n')
    f.write(f'new columns (purified immune signatures, ssGSEA z-scored within cohort):\n')
    for c in new_cols:
        nmiss = merged[c].isna().sum()
        f.write(f'  {c}  (n_missing = {nmiss})\n')
    f.write(f'\nsubjects missing pre-treatment RNA (NaN for new features): {sorted(miss, key=lambda x: int(x))}\n')
    f.write(f'  (median imputation in nested-CV pipeline will fill these)\n')

print(open(f'{OUT}/master_v2_diff.txt').read())
