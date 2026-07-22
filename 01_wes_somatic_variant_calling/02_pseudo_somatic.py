"""
Pseudo-somatic variant extraction from Macrogen germline annotated xlsx.
Strategy:
  - Tumor (PR/PO/P) variants MINUS matched Normal (N) variants at same CHR:POS:REF:ALT
  - Filter: FILTER==PASS, gnomAD AF < 0.001 (or NA = novel), DP >= 10, AD_alt >= 3
  - For unmatched tumors: skip T-N subtraction; rely on gnomAD+1000G+ESP rarity
Outputs per subject:
  02_wes_somatic/variants/{sample}.somatic.tsv
  02_wes_tmb_msi/tmb_summary.tsv
"""
import pandas as pd, numpy as np, os, re, sys
from pathlib import Path

WES_DIR = Path('/mnt/sda1/data/TNT/TNT_WES')
INV = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')
OUT = Path('/mnt/sda1/data/TNT/analysis/02_wes_tmb_msi'); OUT.mkdir(parents=True, exist_ok=True)
VARS = OUT.parent / '02_wes_somatic_pseudo'; VARS.mkdir(parents=True, exist_ok=True)

# SureSelect V5 exome ~50Mb target
EXOME_MB = 50.0

# Fields to read from annotated xlsx (large file with 195 cols — load only what we need)
KEEP_COLS = ['CHROM','POS','REF','ALT','DP','AD','QUAL','FILTER','Zygosity',
             'Effect','Putative_Impact','Gene_Name','HGVS.c','HGVS.p',
             'Transcript_BioType','AA_pos',
             'dbSNP156_ID','p3_1000G_AF','ESP6500_MAF_ALL',
             'CLINVAR_CLNSIG','gnomAD_exomes_AF','gnomAD_exomes_EAS_AF',
             'SIFT_pred','PROVEAN_pred','FATHMM_pred','MetaSVM_pred','MetaLR_pred']

# Somatic coding effects (nonsynonymous)
CODING_EFFECTS_NONSYN = {
  'missense_variant', 'stop_gained', 'stop_lost', 'start_lost',
  'frameshift_variant', 'inframe_insertion', 'inframe_deletion',
  'splice_acceptor_variant', 'splice_donor_variant',
  'protein_altering_variant', 'initiator_codon_variant',
}
# Looser set incl silent for signature analysis (SBS signatures use all SNVs)
ALL_CODING = CODING_EFFECTS_NONSYN | {'synonymous_variant', 'stop_retained_variant'}

def parse_ad(ad):
    try:
        parts = str(ad).split(',')
        ref_c = int(parts[0]); alt_c = int(parts[1])
        return ref_c, alt_c
    except: return None, None

def parse_float(x):
    try:
        if pd.isna(x) or str(x).strip() in ('','.','NA'): return np.nan
        return float(x)
    except: return np.nan

def load_variants(sample_row):
    p = Path(sample_row['dir'])
    sid = sample_row['sample_id']  # e.g., '1-PR'
    xlsx = p / f'{sid}_DNA_SNP_Indel_ANNO.xlsx'
    if not xlsx.exists():
        print(f'  missing xlsx: {xlsx}', file=sys.stderr); return None
    df = pd.read_excel(xlsx, usecols=lambda c: c in KEEP_COLS)
    df['CHROM'] = df['CHROM'].astype(str)
    df['POS'] = df['POS'].astype(int)
    df['REF'] = df['REF'].astype(str)
    df['ALT'] = df['ALT'].astype(str)
    df['key'] = df['CHROM']+':'+df['POS'].astype(str)+':'+df['REF']+':'+df['ALT']
    ad = df['AD'].apply(parse_ad)
    df['AD_ref'] = [x[0] for x in ad]
    df['AD_alt'] = [x[1] for x in ad]
    df['VAF'] = df.apply(lambda r: (r['AD_alt']/(r['AD_ref']+r['AD_alt'])) if r['AD_ref'] is not None and (r['AD_ref']+r['AD_alt'])>0 else np.nan, axis=1)
    df['gnomAD_AF'] = df['gnomAD_exomes_AF'].apply(parse_float)
    df['g1000_AF'] = df['p3_1000G_AF'].apply(parse_float)
    df['esp_AF'] = df['ESP6500_MAF_ALL'].apply(parse_float)
    df['sample'] = sid
    return df

def apply_somatic_filters(tumor_df, normal_keys=None, af_max=0.001, min_dp=10, min_alt=3, min_vaf=0.05, max_normal_vaf=0.02):
    df = tumor_df.copy()
    pass_mask = (df['FILTER'].astype(str).str.contains('PASS', na=False))
    # population AF filter: all three must be <af_max or NA (novel)
    af_mask = ((df['gnomAD_AF'].isna() | (df['gnomAD_AF']<af_max)) &
               (df['g1000_AF'].isna()  | (df['g1000_AF']<af_max)) &
               (df['esp_AF'].isna()    | (df['esp_AF']<af_max)))
    depth_mask = (df['DP']>=min_dp) & (df['AD_alt']>=min_alt) & (df['VAF']>=min_vaf)
    somatic = df[pass_mask & af_mask & depth_mask].copy()
    # T-N subtraction
    if normal_keys is not None and len(normal_keys)>0:
        somatic = somatic[~somatic['key'].isin(normal_keys)]
    return somatic

def main():
    # Build subject -> normal keys
    normals = {}
    for _, r in INV[INV.timepoint=='normal'].iterrows():
        try:
            ndf = load_variants(r)
            if ndf is None: continue
            # normal keys to subtract (loose: any variant in N, regardless of AF)
            nkeys = set(ndf['key'])
            normals[int(r['subject_id'])] = nkeys
            print(f'  normal subj {r["subject_id"]}: {len(nkeys)} variants')
        except Exception as e:
            print(f'  err normal {r["sample_id"]}: {e}', file=sys.stderr)

    # Process tumor (PR/PO/P)
    summary = []
    tumors = INV[INV.timepoint.isin(['pre','post'])]
    for _, r in tumors.iterrows():
        sid = r['sample_id']; subj = int(r['subject_id'])
        try:
            tdf = load_variants(r)
            if tdf is None:
                summary.append((sid, subj, r['timepoint'], r['response_bin'], 'no_xlsx', 0,0,0,0)); continue
            n_total = len(tdf)
            nkeys = normals.get(subj)
            matched = nkeys is not None
            som = apply_somatic_filters(tdf, normal_keys=nkeys)
            # Add effect classification
            som['is_coding'] = som['Effect'].astype(str).apply(lambda e: any(c in e for c in ALL_CODING))
            som['is_nonsyn'] = som['Effect'].astype(str).apply(lambda e: any(c in e for c in CODING_EFFECTS_NONSYN))
            n_som = len(som); n_cod = som['is_coding'].sum(); n_ns = som['is_nonsyn'].sum()
            tmb_ns = n_ns / EXOME_MB
            tmb_all = n_som / EXOME_MB
            som.to_csv(VARS/f'{sid}.somatic.tsv', sep='\t', index=False)
            summary.append((sid, subj, r['timepoint'], r['response_bin'], 'matched' if matched else 'unmatched',
                           n_total, n_som, n_cod, n_ns, tmb_ns, tmb_all))
            print(f'  tumor {sid} subj{subj} {"matched" if matched else "TUMOR-ONLY"}: total={n_total} somatic={n_som} coding={n_cod} nonsyn={n_ns} TMB_ns={tmb_ns:.2f}')
        except Exception as e:
            print(f'  err tumor {sid}: {e}', file=sys.stderr)
            summary.append((sid, subj, r['timepoint'], r['response_bin'], f'err:{e}',0,0,0,0,0,0))
    cols = ['sample_id','subject_id','timepoint','response_bin','mode','n_total_variants','n_pseudo_somatic','n_coding','n_nonsyn','TMB_nonsyn_per_Mb','TMB_all_somatic_per_Mb']
    pd.DataFrame(summary, columns=cols).to_csv(OUT/'tmb_summary.tsv', sep='\t', index=False)
    print('\nSaved:', OUT/'tmb_summary.tsv')

if __name__=='__main__':
    main()
