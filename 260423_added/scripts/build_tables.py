"""Build 5 supplementary tables (xlsx) for §3.6/§3.7."""
import pandas as pd
from pathlib import Path

SRC = Path('/mnt/sda1/data/TNT/analysis/260423_added/source_data')
OUT = Path('/mnt/sda1/data/TNT/analysis/260423_added/tables')
OUT.mkdir(parents=True, exist_ok=True)

# ============ §3.6 tables ============
# 1. neoantigen per-sample (v2)
per_sample = pd.read_csv(SRC/'neo_v2_per_sample.tsv', sep='\t')
per_sample = per_sample.rename(columns={
    'sample_id':'Sample ID', 'subject_id':'Subject',
    'timepoint':'Timepoint', 'response':'TNT response',
    'n_binders':'Total binder peptides (IC50<500nM)',
    'n_binder_sites':'Mutation sites with ≥1 binder',
    'n_strong_binders':'Strong binders (IC50<50nM)',
    'PCN':'Peptide copy number (PCN)',
    'n_variants':'Variants passed to pVACseq',
})
with pd.ExcelWriter(OUT/'Table_36_neoantigen_per_sample_v2.xlsx',
                     engine='openpyxl') as w:
    per_sample.to_excel(w, sheet_name='per_sample', index=False)

# 2. preCRT summary
pre_sum = pd.read_csv(SRC/'neo_v2_preCRT_summary.tsv', sep='\t')
pre_sum = pre_sum[['label','n_good','n_bad','good_median','bad_median',
                    'good_mean','bad_mean','MW_p_twosided']].rename(columns={
    'label':'Metric','n_good':'n good','n_bad':'n bad',
    'good_median':'Good median','bad_median':'Bad median',
    'good_mean':'Good mean','bad_mean':'Bad mean',
    'MW_p_twosided':'MW P (two-sided)'})
with pd.ExcelWriter(OUT/'Table_36_preCRT_summary_v2.xlsx',
                     engine='openpyxl') as w:
    pre_sum.to_excel(w, sheet_name='preCRT_summary', index=False)

# 3. paired Δ summary
pair_sum = pd.read_csv(SRC/'neo_v2_paired_delta_summary.tsv', sep='\t')
pair_sum = pair_sum[['label','n_good','n_bad','good_median','bad_median',
                      'good_mean','bad_mean','MW_p_twosided',
                      'wilcoxon_good_p','wilcoxon_bad_p']].rename(columns={
    'label':'Metric','n_good':'n good','n_bad':'n bad',
    'good_median':'Good Δ median','bad_median':'Bad Δ median',
    'good_mean':'Good Δ mean','bad_mean':'Bad Δ mean',
    'MW_p_twosided':'MW P (two-sided)',
    'wilcoxon_good_p':'Within-good Wilcoxon P',
    'wilcoxon_bad_p':'Within-bad Wilcoxon P'})
bca = pd.read_csv(SRC/'neo_v2_bca_ci.tsv', sep='\t')
bca = bca.rename(columns={'label':'Metric','group':'Group','n':'n',
                          'median':'Δ median','ci_lo':'BCa 95% CI lower',
                          'ci_hi':'BCa 95% CI upper'})[
    ['Metric','Group','n','Δ median','BCa 95% CI lower','BCa 95% CI upper']]
pair_raw = pd.read_csv(SRC/'neo_v2_paired_delta.tsv', sep='\t')
pair_raw = pair_raw.rename(columns={
    'subject_id':'Subject','response':'TNT response',
    'pre_binders':'Pre binders','post_binders':'Post binders',
    'pre_sites':'Pre sites','post_sites':'Post sites',
    'pre_strong':'Pre strong','post_strong':'Post strong',
    'pre_PCN':'Pre PCN','post_PCN':'Post PCN',
    'delta_binders':'Δ binders','delta_sites':'Δ sites',
    'delta_strong':'Δ strong','delta_PCN':'Δ PCN'})
with pd.ExcelWriter(OUT/'Table_36_paired_delta_summary.xlsx',
                     engine='openpyxl') as w:
    pair_sum.to_excel(w, sheet_name='summary_stats', index=False)
    bca.to_excel     (w, sheet_name='BCa_bootstrap_CI', index=False)
    pair_raw.to_excel(w, sheet_name='per_subject_raw', index=False)

# ============ §3.7 tables ============
# 4. HLA-LOH per-locus (strict + lite side-by-side)
loc = pd.read_csv(SRC/'hla_loh_per_locus.tsv', sep='\t')
loc_pivot = loc.pivot_table(
    index='locus', columns='call',
    values=['good_pos','good_n_het','bad_pos','bad_n_het','good_freq','bad_freq']
).reset_index()
loc_pivot.columns = ['_'.join([str(x) for x in c if x]) if isinstance(c, tuple) else c
                     for c in loc_pivot.columns]
# also include subject-level summary + tests
subj = pd.read_csv(SRC/'hla_loh_subject_summary.tsv', sep='\t')
subj = subj.rename(columns={
    'subject_id':'Subject','response_bin':'TNT response',
    'n_het_loci':'# het loci','n_loh_lite':'# lite LOH',
    'n_loh_strict':'# strict LOH','any_loh_lite':'Any lite LOH',
    'any_loh_strict':'Any strict LOH'})
tests = pd.read_csv(SRC/'hla_loh_tests.tsv', sep='\t').rename(columns={
    'call':'Call','good_pos':'Good + ','good_neg':'Good − ',
    'bad_pos':'Bad + ','bad_neg':'Bad − ',
    'good_freq':'Good freq','bad_freq':'Bad freq',
    'OR':'Odds ratio','fisher_p':'Fisher P (two-sided)'})
clear = pd.read_csv(SRC/'hla_loh_clearance.tsv', sep='\t')
clear = clear.rename(columns={
    'subject_id':'Subject','locus':'Locus',
    'pre_ratio':'Pre tumor ratio','post_ratio':'Post tumor ratio',
    'pre_imbalance':'Pre imbalance','post_imbalance':'Post imbalance',
    'pre_strict':'Pre strict LOH','post_strict':'Post strict LOH',
    'pre_lite':'Pre lite LOH','post_lite':'Post lite LOH',
    'response_bin':'TNT response'})
with pd.ExcelWriter(OUT/'Table_37_HLA_LOH_per_locus.xlsx',
                     engine='openpyxl') as w:
    loc.rename(columns={'locus':'Locus','call':'Call',
                        'good_pos':'Good LOH+','good_n_het':'Good het total',
                        'bad_pos':'Bad LOH+','bad_n_het':'Bad het total',
                        'good_freq':'Good freq','bad_freq':'Bad freq'}
               ).to_excel(w, sheet_name='per_locus_long', index=False)
    subj.to_excel (w, sheet_name='per_subject_summary', index=False)
    tests.to_excel(w, sheet_name='Fisher_tests_any_LOH', index=False)
    clear.to_excel(w, sheet_name='paired_clearance', index=False)

# 5. homozygosity
hom = pd.read_csv(SRC/'hla_homozygosity.tsv', sep='\t')
hom = hom.rename(columns={'subject_id':'Subject','response_bin':'TNT response',
                          'homozygous_A':'Homozygous HLA-A',
                          'homozygous_B':'Homozygous HLA-B',
                          'homozygous_C':'Homozygous HLA-C',
                          'n_homozygous_loci':'# homozygous loci',
                          'any_homozygous':'Any homozygous'})
hom_t = pd.read_csv(SRC/'hla_homozygosity_tests.tsv', sep='\t').rename(columns={
    'metric':'Metric','good_pos':'Good + ','good_neg':'Good − ',
    'bad_pos':'Bad + ','bad_neg':'Bad − ',
    'good_freq':'Good freq','bad_freq':'Bad freq',
    'OR':'Odds ratio','fisher_p':'Fisher P (two-sided)'})
with pd.ExcelWriter(OUT/'Table_37_homozygosity.xlsx',
                     engine='openpyxl') as w:
    hom.to_excel  (w, sheet_name='per_subject', index=False)
    hom_t.to_excel(w, sheet_name='Fisher_tests', index=False)

print('Tables written:')
for f in sorted(OUT.glob('*.xlsx')):
    print(' ', f.name, f.stat().st_size, 'bytes')
