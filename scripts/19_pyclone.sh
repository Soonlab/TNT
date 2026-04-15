#!/bin/bash
# PyClone-VI clonal evolution analysis for 13 paired pre+post subjects
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate pyclone
# Use wes_somatic bcftools (pyclone env has broken libgsl)
export PATH=/home/soon/miniconda3/envs/wes_somatic/bin:$PATH

OUT=/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/pyclone
mkdir -p $OUT

# Input preparation: from Mutect2 PASS + CNVkit .call.cns segments
# PyClone-VI input format: mutation_id sample_id ref_counts alt_counts major_cn minor_cn normal_cn tumour_content
python3 << 'PYEOF'
import pandas as pd, numpy as np, os, gzip
from pathlib import Path

M2 = Path('/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered')
CNV = Path('/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/cnvkit')
OUT = Path('/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/pyclone')

# Paired subjects (Y group with both PR+PO): from metadata, 14 subjects (1-14)
# subj 13 lacks N so PR/PO were tumor-only — include but tag separately
PAIRED = list(range(1, 15))

import subprocess as sp
def parse_vcf_pass(path, sid):
    # Use bcftools to extract CHR,POS,REF,ALT,AD,DP
    res = sp.run(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT\t[%AD\t%DP\t%AF]\n', str(path)],
                 capture_output=True, text=True)
    rows=[]
    for line in res.stdout.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        if len(parts)<7: continue
        chrom,pos,ref,alt,ad,dp,af = parts
        if len(ref)>1 or len(alt)>1: continue  # SNV only
        try:
            ad1 = int(ad.split(',')[0]); ad2 = int(ad.split(',')[1])
            dp_i = int(dp); af_f = float(af.split(',')[0])
            if dp_i < 10 or ad2 < 3: continue
        except: continue
        rows.append({'mutation_id':f'{chrom}:{pos}:{ref}:{alt}','chrom':chrom,'pos':int(pos),
                     'ref_counts':ad1,'alt_counts':ad2,'sample_id':sid,'vaf':af_f})
    return pd.DataFrame(rows)

def load_cn(sid):
    f = CNV / f'{sid}_DNA.call.cns'
    if not f.exists(): return None
    c = pd.read_csv(f, sep='\t')
    return c

def assign_cn(muts, cn):
    # assign CN to each mut by position
    out=[]
    for _, m in muts.iterrows():
        seg = cn[(cn['chromosome']==m['chrom']) & (cn['start']<=m['pos']) & (cn['end']>=m['pos'])]
        if len(seg)==0:
            major, minor = 1, 1
        else:
            tot = int(seg.iloc[0]['cn'])
            major = max(1, tot-1); minor = 0 if tot<=1 else 1
        out.append({'major_cn':major,'minor_cn':minor})
    return pd.DataFrame(out)

all_rows=[]
for subj in PAIRED:
    for tp_sid in [f'{subj}-PR', f'{subj}-PO']:
        vcf = M2 / f'{tp_sid}.pass.vcf.gz'
        if not vcf.exists(): continue
        muts = parse_vcf_pass(vcf, tp_sid)
        if len(muts)==0: continue
        cn = load_cn(tp_sid)
        if cn is None: continue
        cn_add = assign_cn(muts, cn)
        muts = pd.concat([muts.reset_index(drop=True), cn_add], axis=1)
        muts['normal_cn']=2; muts['tumour_content']=0.5
        all_rows.append(muts)

all_df = pd.concat(all_rows, ignore_index=True)
# Only keep muts that are in both PR AND PO of same subject
def subj_from_sid(s): return int(s.split('-')[0])
all_df['subject_id'] = all_df['sample_id'].apply(subj_from_sid)
keep_muts=[]
for subj, sub in all_df.groupby('subject_id'):
    samps = set(sub['sample_id'])
    if len(samps) < 2: continue
    # union of mutations across pre/post
    pyin = sub[['mutation_id','sample_id','ref_counts','alt_counts','major_cn','minor_cn','normal_cn','tumour_content']]
    pyin.to_csv(OUT/f'pyclone_input_subj{subj}.tsv', sep='\t', index=False)
PYEOF

# Run PyClone-VI per subject
for f in $OUT/pyclone_input_subj*.tsv; do
  subj=$(basename $f .tsv | sed 's/pyclone_input_subj//')
  H5=$OUT/fit_subj${subj}.h5
  RES=$OUT/results_subj${subj}.tsv
  [ -s "$RES" ] && continue
  n=$(($(wc -l < $f)-1))
  [ $n -lt 30 ] && { echo "subj$subj too few ($n)"; continue; }
  echo "[$(date +%T)] PyClone-VI subj $subj n=$n"
  pyclone-vi fit -i $f -o $H5 -c 10 -d beta-binomial -r 10 2>&1 | tail -2
  pyclone-vi write-results-file -i $H5 -o $RES 2>&1 | tail -2
done
echo "=== done ==="
ls $OUT/results_subj*.tsv 2>/dev/null | wc -l
