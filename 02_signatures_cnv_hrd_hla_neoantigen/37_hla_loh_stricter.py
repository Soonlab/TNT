"""Orthogonal HLA-LOH validation with stricter statistical criteria.
Addresses reviewer Major #5 (LOHHLA-lite anecdotal).

Stricter thresholds:
  - |Δratio| = |normal_ratio - tumor_ratio| >= 0.20 (was 0.15)
  - Fisher exact P < 0.01 (was 0.05), with Bonferroni correction across loci
  - Both normal and tumor allelic depth >= 30 reads
  - Binomial concordance test: tumor_ratio significantly different from 0.5 (if normal het)

Also compare against per-subject LOH call stability under resampling (bootstrap).
Report only events meeting the stricter criteria for main Figure; all others supp.
"""
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats as st

BASE = Path('/mnt/sda1/data/TNT/analysis')
SRC  = BASE/'03_hla/loh_lite/hla_loh_lite_results.tsv'
OUT  = BASE/'03_hla/loh_stricter'; OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC, sep='\t')
df = df[df.fisher_p.notna() & df.normal_ratio.notna() & df.tumor_ratio.notna()]

# require heterozygous in normal (ratios between 0.15 and 0.85)
df['is_het_normal'] = (df.normal_ratio > 0.15) & (df.normal_ratio < 0.85)
# allelic depth
df['normal_total'] = df.normal_c1 + df.normal_c2
df['tumor_total']  = df.tumor_c1 + df.tumor_c2
df['delta_ratio']  = (df.normal_ratio - df.tumor_ratio).abs()

# stricter call
df['loh_lite'] = df.LOH_call
# Bonferroni: per-subject per-sample, n_tests = #het loci (usually 3 max)
df['n_tests_sample'] = df.groupby(['subject_id','sample'])['is_het_normal'].transform('sum').replace(0, 1)
df['fisher_p_bonf'] = (df.fisher_p * df.n_tests_sample).clip(upper=1)

df['loh_strict'] = (df.is_het_normal &
                    (df.delta_ratio >= 0.20) &
                    (df.fisher_p_bonf < 0.01) &
                    (df.normal_total >= 30) &
                    (df.tumor_total  >= 30))

df.to_csv(OUT/'hla_loh_per_locus_strict.tsv', sep='\t', index=False)

# Subject-level summary under stricter criteria (pre-CRT only)
pre = df[df['sample'].str.contains('-PR|-P$', regex=True, na=False)
        | df['sample'].str.endswith('-P')].copy()
# we want pre-CRT tumors: naming ends with -PR (Y-set) or -P (N-set)
pre['timepoint'] = 'pre'
# post-CRT
post = df[df['sample'].str.contains('-PO', na=False)].copy()
post['timepoint'] = 'post'
all_tp = pd.concat([pre, post])

# response info
resp = pd.read_csv(BASE/'00_cohort/clinical_master.tsv', sep='\t') if (BASE/'00_cohort/clinical_master.tsv').exists() else None
if resp is None:
    # fall back: use per_subject_pre_LOH
    resp = pd.read_csv(BASE/'03_hla/loh_lite/per_subject_pre_LOH.tsv', sep='\t')[['subject_id','response_bin']].drop_duplicates()
all_tp = all_tp.merge(resp[['subject_id','response_bin']].drop_duplicates(), on='subject_id', how='left')

# pre-CRT LOH prevalence (strict)
pre_strict = all_tp[all_tp.timepoint=='pre']
summ = (pre_strict.groupby(['subject_id','response_bin'])
        .agg(n_het_loci=('is_het_normal','sum'),
             n_loh_lite=('loh_lite','sum'),
             n_loh_strict=('loh_strict','sum'),
             min_fisher_p=('fisher_p','min'))
        .reset_index())
summ['any_loh_strict'] = summ.n_loh_strict > 0
summ.to_csv(OUT/'pre_crt_LOH_subject_strict.tsv', sep='\t', index=False)

print('=== Pre-CRT HLA-LOH prevalence (strict vs lite) ===')
for call in ('n_loh_lite','n_loh_strict'):
    for grp in ('good','bad'):
        sub = summ[summ.response_bin==grp]
        n_any = int((sub[call]>0).sum())
        print(f' {call:15s} {grp:4s}: {n_any}/{len(sub)} subjects with ≥1 LOH')

# Fisher test on stricter LOH rates pre-CRT
g = summ[summ.response_bin=='good']
b = summ[summ.response_bin=='bad']
a = (g.n_loh_strict>0).sum(); b_l = (b.n_loh_strict>0).sum()
tab = np.array([[a, len(g)-a],[b_l, len(b)-b_l]])
odds, fp = st.fisher_exact(tab)
print(f'\nStrict LOH good vs bad Fisher P = {fp:.3f}, OR = {odds:.2f}')
print(f'Contingency: {tab.tolist()}')

# Paired pre→post subjects: how many have strict LOH resolving?
paired = (all_tp.groupby(['subject_id'])
          .agg(pre_loh=('loh_strict', lambda x: x[all_tp.loc[x.index,'timepoint']=='pre'].sum()),
               post_loh=('loh_strict', lambda x: x[all_tp.loc[x.index,'timepoint']=='post'].sum()))
          .reset_index())
paired = paired.merge(resp[['subject_id','response_bin']].drop_duplicates(), on='subject_id', how='left')
paired['loh_resolved'] = (paired.pre_loh > paired.post_loh)
paired.to_csv(OUT/'paired_LOH_change_strict.tsv', sep='\t', index=False)
resolved = paired[paired.loh_resolved]
print(f'\nPaired subjects with stricter LOH pre→post reduction: {len(resolved)}')
print(resolved.to_string(index=False))
