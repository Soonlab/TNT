"""Prepare data for Fig 8 v07-refresh panels (v2 cohort).

Adds to 260423_added/source_data/:
  - hla_allele_topN.tsv            : top-10 allele freq per locus (A/B/C) with hom/het split
  - hla_loh_heatmap_long.tsv       : per-subject × per-locus LOH call for heatmap (Fig 8C)
  - hla_loh_per_allele_fisher.tsv  : per-allele Fisher good vs bad contingency (Fig 8C bar)
  - neo_v2_lollipop.tsv            : per-subject pre binders + strong + HLA zygosity/LOH annotation (Fig 8F)
  - neo_v2_paired_per_subject.tsv  : per-subject paired pre/post sites + strong (Fig 8E traces)
  - hla_homozygosity_mw.tsv        : MW test on n_homozygous_loci (Fig 8B frame — matches old Fig 8B)
"""
import pandas as pd, numpy as np
from pathlib import Path
from collections import Counter
from scipy import stats as st

BASE = Path('/mnt/sda1/data/TNT/analysis')
OUT  = BASE / '260423_added' / 'source_data'

hla   = pd.read_csv(BASE/'03_hla/hla_class_I_typing.tsv', sep='\t')
loh   = pd.read_csv(BASE/'03_hla/loh_stricter/hla_loh_per_locus_strict.tsv', sep='\t')
clin  = pd.read_csv(BASE/'00_cohort/clinical_master.tsv', sep='\t')
per_sample = pd.read_csv(OUT/'neo_v2_per_sample.tsv', sep='\t')

# ============ Fig 8A: HLA allele top-N with hom/het split ============
rows = []
for locus, cols in [('HLA-A', ['A1','A2']), ('HLA-B', ['B1','B2']), ('HLA-C', ['C1','C2'])]:
    counter_het = Counter()
    counter_hom = Counter()
    for _, r in hla.iterrows():
        a1 = r[cols[0]]; a2 = r[cols[1]]
        if pd.isna(a1): a1 = ''
        if pd.isna(a2): a2 = ''
        is_hom = (a1 == a2) or (a2 == '')
        alleles = [a for a in [a1, a2] if a]
        for a in alleles:
            if is_hom:
                counter_hom[a] += 1
            else:
                counter_het[a] += 1
    # combine
    all_alleles = set(counter_het) | set(counter_hom)
    for a in all_alleles:
        rows.append({'locus': locus, 'allele': a,
                     'het_count': counter_het[a], 'hom_count': counter_hom[a],
                     'total': counter_het[a] + counter_hom[a]})
allele_df = pd.DataFrame(rows).sort_values(['locus','total'], ascending=[True, False])
allele_df.to_csv(OUT/'hla_allele_topN.tsv', sep='\t', index=False)

# ============ Fig 8B: MW test on n_homozygous_loci ============
hom_subj = hla[['subject_id','n_homozygous_loci','response_bin']].copy()
g = hom_subj[hom_subj.response_bin == 'good']['n_homozygous_loci'].values
b = hom_subj[hom_subj.response_bin == 'bad']['n_homozygous_loci'].values
mw_p = st.mannwhitneyu(g, b, alternative='two-sided').pvalue
pd.DataFrame([{
    'test':'Mann-Whitney U two-sided on n_homozygous_loci',
    'n_good': len(g), 'n_bad': len(b),
    'good_mean': float(np.mean(g)), 'bad_mean': float(np.mean(b)),
    'good_median': float(np.median(g)), 'bad_median': float(np.median(b)),
    'MW_p': mw_p
}]).to_csv(OUT/'hla_homozygosity_mw.tsv', sep='\t', index=False)
hom_subj.to_csv(OUT/'hla_homozygosity_per_subject.tsv', sep='\t', index=False)

# ============ Fig 8C: LOH heatmap + per-allele Fisher ============
# Pre-CRT het loci only → binary LOH call (use lite for main heatmap to match original visual)
loh['timepoint'] = loh['sample'].apply(
    lambda s: 'pre' if (s.endswith('-PR') or s.endswith('-P')) else
              ('post' if s.endswith('-PO') else 'other'))
pre = loh[(loh.timepoint == 'pre') & loh.is_het_normal].copy()
pre = pre.merge(clin[['subject_id','response_bin']], on='subject_id', how='left')
heat = pre[['subject_id','response_bin','locus','loh_lite','loh_strict',
            'normal_ratio','tumor_ratio']].copy()
heat.to_csv(OUT/'hla_loh_heatmap_long.tsv', sep='\t', index=False)

# Per-allele Fisher
g_tests = pre[pre.response_bin=='good']
b_tests = pre[pre.response_bin=='bad']
rows = []
for call_col, label in [('loh_lite','LOHHLA-lite'), ('loh_strict','strict')]:
    g_pos = int(g_tests[call_col].sum()); g_n = len(g_tests)
    b_pos = int(b_tests[call_col].sum()); b_n = len(b_tests)
    odds, p = st.fisher_exact([[g_pos, g_n-g_pos], [b_pos, b_n-b_pos]])
    rows.append({'call': label, 'good_pos': g_pos, 'good_total': g_n,
                 'bad_pos': b_pos,  'bad_total': b_n,
                 'good_pct': g_pos/g_n*100, 'bad_pct': b_pos/b_n*100,
                 'OR': odds, 'fisher_p': p})
pd.DataFrame(rows).to_csv(OUT/'hla_loh_per_allele_fisher.tsv', sep='\t', index=False)

# ============ Fig 8E: per-subject paired pre/post for 2 metrics ============
paired_long = []
for sid in [1,2,3,4,5,6,7,8,9,10,11,12,14]:
    sub = per_sample[per_sample.subject_id == sid]
    pre_r = sub[sub.timepoint=='pre']; post_r = sub[sub.timepoint=='post']
    if len(pre_r)!=1 or len(post_r)!=1: continue
    pr, po = pre_r.iloc[0], post_r.iloc[0]
    paired_long.append({
        'subject_id': sid, 'response': pr['response'],
        'pre_sites': pr.n_binder_sites, 'post_sites': po.n_binder_sites,
        'pre_strong': pr.n_strong_binders, 'post_strong': po.n_strong_binders,
    })
pd.DataFrame(paired_long).to_csv(OUT/'neo_v2_paired_per_subject.tsv', sep='\t', index=False)

# ============ Fig 8F: per-subject lollipop w/ annotations ============
# Pre-CRT tumor per subject (one row per subject). Use per_sample pre rows.
pre_per = per_sample[per_sample.timepoint == 'pre'].copy()
# Merge with HLA homozygosity per locus
pre_per = pre_per.merge(
    hla[['subject_id','homozygous_A','homozygous_B','homozygous_C']],
    on='subject_id', how='left')
# LOH flags per subject (strict/lite per locus → bool)
loh_per_subj = loh[(loh.timepoint=='pre') & loh.is_het_normal].groupby(['subject_id','locus']).agg(
    loh_strict=('loh_strict','any'),
    loh_lite=('loh_lite','any')).reset_index()
# pivot to wide per locus
loh_strict_wide = loh_per_subj.pivot(index='subject_id', columns='locus', values='loh_strict').fillna(False)
loh_lite_wide   = loh_per_subj.pivot(index='subject_id', columns='locus', values='loh_lite').fillna(False)
loh_strict_wide.columns = [f'loh_strict_{c.replace("HLA-","")}' for c in loh_strict_wide.columns]
loh_lite_wide.columns   = [f'loh_lite_{c.replace("HLA-","")}'   for c in loh_lite_wide.columns]
loh_strict_wide = loh_strict_wide.reset_index()
loh_lite_wide   = loh_lite_wide.reset_index()
pre_per = pre_per.merge(loh_strict_wide, on='subject_id', how='left')
pre_per = pre_per.merge(loh_lite_wide,   on='subject_id', how='left')
# Any LOH flag
for c in ['loh_strict_A','loh_strict_B','loh_strict_C','loh_lite_A','loh_lite_B','loh_lite_C']:
    if c not in pre_per.columns: pre_per[c] = False
pre_per['any_loh_strict'] = pre_per[['loh_strict_A','loh_strict_B','loh_strict_C']].any(axis=1)
pre_per['any_loh_lite']   = pre_per[['loh_lite_A','loh_lite_B','loh_lite_C']].any(axis=1)
pre_per = pre_per.sort_values(['response','n_binder_sites']).reset_index(drop=True)
pre_per.to_csv(OUT/'neo_v2_lollipop.tsv', sep='\t', index=False)

print('v07-refresh data prep done.  Files:')
for f in ['hla_allele_topN.tsv','hla_homozygosity_mw.tsv','hla_homozygosity_per_subject.tsv',
          'hla_loh_heatmap_long.tsv','hla_loh_per_allele_fisher.tsv',
          'neo_v2_paired_per_subject.tsv','neo_v2_lollipop.tsv']:
    p = OUT/f; print(' ', f, p.stat().st_size if p.exists() else 'MISSING', 'bytes')
