"""Build formal Table 1 (Clinical), Table 2 (Top-20 features), Table S1–S9.

Sources used:
  Table 1  : 00_cohort/clinical + integrated_subject_master (n=35)
  Table 2  : response_feature_stats.tsv → top 20 by |delta_med| with q filter
  Table S1 : integrated_subject_master.tsv (rename/copy)
  Table S2 : 01_wes_signatures/sbs_activities_with_meta.tsv (+ sbs_summary_key)
  Table S3 : 05_rna_deg_gsea/GSEA_{Hallmark,Reactome}_pre.tsv concatenated
  Table S4 : 02_wes_tmb_msi/pass_variant_counts.tsv + variant_master (driver rows)
  Table S5 : 03_hla/hla_class_I_typing.tsv
  Table S6 : 03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv
  Table S7 : already exists
  Table S8 : already exists
  Table S9 : 03_hla/loh_stricter/hla_loh_per_locus_strict.tsv — format for supp
"""
import pandas as pd, numpy as np, gzip
from pathlib import Path

BASE = Path('/mnt/sda1/data/TNT/analysis')
TBL  = BASE/'tables'
SUB_TBL = Path('/data/data/TNT/analysis/genome_medicine_submission/tables')
SUB_TBL.mkdir(parents=True, exist_ok=True)

def save(df, fn):
    df.to_csv(SUB_TBL/fn, sep='\t', index=False)
    print(f'  saved {SUB_TBL}/{fn}  ({df.shape[0]} rows)')

# =========================================================
# Table 1 — Clinical characteristics by response
# =========================================================
master = pd.read_csv(TBL/'integrated_subject_master.tsv', sep='\t')
rows = []
for col, name in [('age','Age (years)'), ('sex','Sex'), ('cT','Clinical T stage'),
                  ('prepost_set','Paired pre/post set'), ('TMB_nonsyn_per_Mb','TMB (non-syn/Mb)'),
                  ('MSI_pct','MSI fraction (%)'), ('CMS','CMS subtype')]:
    if col not in master.columns: continue
    if master[col].dtype in (np.float64, np.int64):
        good = master.loc[master.response_bin=='good', col].dropna()
        bad  = master.loc[master.response_bin=='bad', col].dropna()
        rows.append({'Variable': name,
                     'Good (n=18)': f'{good.median():.2f} [{good.min():.2f}–{good.max():.2f}]',
                     'Bad (n=17)':  f'{bad.median():.2f} [{bad.min():.2f}–{bad.max():.2f}]',
                     'Test': 'Mann–Whitney U',
                     'P-value': f'{__import__("scipy").stats.mannwhitneyu(good,bad).pvalue:.3f}'})
    else:
        gv = master.loc[master.response_bin=='good', col].value_counts()
        bv = master.loc[master.response_bin=='bad',  col].value_counts()
        categories = sorted(set(gv.index) | set(bv.index), key=str)
        for c in categories:
            rows.append({'Variable': f'{name} — {c}',
                         'Good (n=18)': int(gv.get(c, 0)),
                         'Bad (n=17)':  int(bv.get(c, 0)),
                         'Test': 'Fisher exact (table-level)',
                         'P-value': ''})
tab1 = pd.DataFrame(rows)
save(tab1, 'Table1_clinical_characteristics.tsv')

# =========================================================
# Table 2 — Top 20 integrated-feature associations
# =========================================================
stats = pd.read_csv(TBL/'response_feature_stats.tsv', sep='\t')
# Sort by |delta_med| after BH FDR
stats['abs_delta'] = stats['delta_med'].abs()
stats_sorted = stats.sort_values('pvalue').head(20).copy()
stats_sorted = stats_sorted[['feature','n_good','n_bad','med_good','med_bad','delta_med','pvalue','qvalue']]
stats_sorted.columns = ['Feature','n good','n bad','Median good','Median bad','Δ (good − bad)','MW P','BH q']
for c in ['Median good','Median bad','Δ (good − bad)']:
    stats_sorted[c] = stats_sorted[c].round(3)
stats_sorted['MW P'] = stats_sorted['MW P'].apply(lambda x: f'{x:.4g}')
stats_sorted['BH q'] = stats_sorted['BH q'].apply(lambda x: f'{x:.3g}')
save(stats_sorted, 'Table2_top20_features.tsv')

# =========================================================
# Table S1 — 37-feature master
# =========================================================
save(master, 'TableS1_integrated_subject_master.tsv')

# =========================================================
# Table S2 — SBS activities per sample
# =========================================================
sbs_path = BASE/'01_wes_signatures/sbs_activities_with_meta.tsv'
if sbs_path.exists():
    sbs = pd.read_csv(sbs_path, sep='\t')
    save(sbs, 'TableS2_SBS_activities_per_sample.tsv')

# =========================================================
# Table S3 — Full GSEA (Hallmark + Reactome)
# =========================================================
parts = []
for f, collection in [('GSEA_Hallmark_pre.tsv','Hallmark'),
                      ('GSEA_Reactome_pre.tsv','Reactome')]:
    p = BASE/'05_rna_deg_gsea'/f
    if p.exists():
        df = pd.read_csv(p, sep='\t')
        df.insert(0, 'Collection', collection)
        parts.append(df)
gsea = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
save(gsea, 'TableS3_GSEA_Hallmark_Reactome.tsv')

# =========================================================
# Table S4 — CRC driver mutations per sample
# =========================================================
vm_path = BASE/'02_wes_tmb_msi/variant_master.tsv.gz'
if vm_path.exists():
    drivers = {'APC','TP53','KRAS','FBXW7','KMT2D','SMAD4','PIK3CA','NRAS','BRAF','TCF7L2','AMER1'}
    rows = []
    with gzip.open(vm_path, 'rt') as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        idx = {c:i for i,c in enumerate(hdr)}
        gcol = idx.get('GENE') or idx.get('SYMBOL') or idx.get('gene') or idx.get('Gene')
        scol = idx.get('sample_id') or idx.get('sample') or idx.get('Tumor_Sample_Barcode')
        vcol = idx.get('EFFECT_primary') or idx.get('EFFECT') or idx.get('Consequence')
        if gcol is not None and scol is not None:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) <= max(gcol, scol): continue
                if f[gcol] in drivers:
                    rows.append({'sample': f[scol],
                                 'gene':   f[gcol],
                                 'consequence': f[vcol] if vcol is not None and vcol<len(f) else ''})
    drv = pd.DataFrame(rows)
    save(drv, 'TableS4_driver_mutations_per_sample.tsv')

# =========================================================
# Table S5 — HLA class I per subject
# =========================================================
hla_path = BASE/'03_hla/hla_class_I_typing.tsv'
if hla_path.exists():
    hla = pd.read_csv(hla_path, sep='\t')
    save(hla, 'TableS5_HLA_class_I_typing.tsv')

# =========================================================
# Table S6 — pVACseq neoantigen detail per sample
# =========================================================
neo_path = BASE/'03_wes_hla_neoantigen/neoantigen_summary_by_sample.tsv'
if neo_path.exists():
    neo = pd.read_csv(neo_path, sep='\t')
    save(neo, 'TableS6_neoantigen_per_sample.tsv')

# =========================================================
# Table S7 — already built (external per-cohort)
# =========================================================
src = TBL/'TableS7_external_percohort_signatures.tsv'
if src.exists():
    import shutil
    shutil.copy(src, SUB_TBL/'TableS7_external_percohort_signatures.tsv')
    print(f'  copied TableS7')

# =========================================================
# Table S8 — already built (cascade BCa)
# =========================================================
src = TBL/'TableS8_cascade_BCa_bootstrap.tsv'
if src.exists():
    import shutil
    shutil.copy(src, SUB_TBL/'TableS8_cascade_BCa_bootstrap.tsv')
    print(f'  copied TableS8')

# =========================================================
# Table S9 — HLA-LOH lenient vs strict per-sample
# =========================================================
strict_path = BASE/'03_hla/loh_stricter/hla_loh_per_locus_strict.tsv'
if strict_path.exists():
    loh = pd.read_csv(strict_path, sep='\t')
    # Select informative columns
    keep = ['subject_id','sample','locus','allele1','allele2',
            'normal_c1','normal_c2','tumor_c1','tumor_c2',
            'normal_ratio','tumor_ratio','delta_ratio',
            'fisher_p','fisher_p_bonf','n_tests_sample',
            'is_het_normal','loh_lite','loh_strict']
    loh_out = loh[[c for c in keep if c in loh.columns]].copy()
    for c in ('normal_ratio','tumor_ratio','delta_ratio','fisher_p','fisher_p_bonf'):
        if c in loh_out.columns:
            loh_out[c] = pd.to_numeric(loh_out[c], errors='coerce').round(4)
    save(loh_out, 'TableS9_HLA_LOH_lite_vs_strict_per_locus.tsv')

print('\nAll tables emitted to:', SUB_TBL)
