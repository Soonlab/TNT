#!/bin/bash
set -e
OUT=/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/pyclone
mkdir -p $OUT

RNA_PY=/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin/python3
PYCLONE=/home/soon/miniconda3/envs/pyclone/bin/pyclone-vi
BCF=/home/soon/miniconda3/envs/wes_somatic/bin/bcftools

$RNA_PY - << PYEOF
import pandas as pd, numpy as np, os, subprocess as sp
from pathlib import Path
M2 = Path('/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered')
CNV = Path('/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/cnvkit')
OUT = Path('/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/pyclone')
BCF = '$BCF'

PAIRED = list(range(1, 15))
def parse(path, sid):
    # Use -s with sample name extraction; auto-pick first sample (tumor) — use query -u to unpack
    # Safer: use view with sample list. First get sample names
    sh = sp.run([BCF,'query','-l',str(path)], capture_output=True, text=True)
    samples = sh.stdout.strip().split('\n')
    # Pick tumor sample: first one NOT ending in -N
    tumor_sm = next((s for s in samples if not s.endswith('-N') and not s.endswith('-N_DNA')), samples[0])
    r = sp.run([BCF,'view','-s',tumor_sm, str(path)], capture_output=True, text=True)
    # Then query from stdin
    r2 = sp.run([BCF,'query','-f','%CHROM\t%POS\t%REF\t%ALT\t[%AD\t%DP\t%AF]\n'],
                input=r.stdout, capture_output=True, text=True)
    rows=[]
    for line in r2.stdout.strip().split('\n'):
        if not line: continue
        p = line.split('\t')
        if len(p)<7: continue
        chrom,pos,ref,alt,ad,dp,af = p[:7]
        if len(ref)>1 or len(alt)>1: continue
        try:
            ad1,ad2 = int(ad.split(',')[0]), int(ad.split(',')[1])
            dp_i = int(dp)
            if dp_i<10 or ad2<3: continue
        except: continue
        rows.append({'mutation_id':f'{chrom}:{pos}:{ref}:{alt}','chrom':str(chrom),'pos':int(pos),
                     'ref_counts':ad1,'alt_counts':ad2,'sample_id':sid})
    return pd.DataFrame(rows)
def load_cn(sid):
    f = CNV / f'{sid}_DNA.call.cns'
    if not f.exists(): return None
    return pd.read_csv(f, sep='\t')
def assign_cn(muts, cn):
    out=[]
    for _, m in muts.iterrows():
        seg = cn[(cn['chromosome'].astype(str)==str(m['chrom'])) & (cn['start']<=m['pos']) & (cn['end']>=m['pos'])]
        if len(seg)==0: major, minor = 1, 1
        else:
            tot = int(seg.iloc[0]['cn'])
            major = max(1, tot-1); minor = 0 if tot<=1 else 1
        out.append({'major_cn':major,'minor_cn':minor})
    return pd.DataFrame(out)
for subj in PAIRED:
    combined=[]
    for tp_sid in [f'{subj}-PR', f'{subj}-PO']:
        vcf = M2 / f'{tp_sid}.pass.vcf.gz'
        if not vcf.exists(): continue
        muts = parse(vcf, tp_sid)
        if len(muts)==0: continue
        cn = load_cn(tp_sid)
        if cn is None: continue
        cn_add = assign_cn(muts, cn)
        muts = pd.concat([muts.reset_index(drop=True), cn_add], axis=1)
        muts['normal_cn']=2; muts['tumour_content']=0.5
        combined.append(muts)
    if len(combined)<2: continue
    all_df = pd.concat(combined, ignore_index=True)
    pyin = all_df[['mutation_id','sample_id','ref_counts','alt_counts','major_cn','minor_cn','normal_cn','tumour_content']]
    pyin.to_csv(OUT/f'pyclone_input_subj{subj}.tsv', sep='\t', index=False)
    print(f'subj {subj}: {len(pyin)} rows, {pyin.mutation_id.nunique()} unique muts')
PYEOF

for f in $OUT/pyclone_input_subj*.tsv; do
    subj=$(basename $f .tsv | sed 's/pyclone_input_subj//')
    H5=$OUT/fit_subj${subj}.h5
    RES=$OUT/results_subj${subj}.tsv
    [ -s "$RES" ] && continue
    n=$(($(wc -l < $f)-1))
    [ $n -lt 30 ] && { echo "subj$subj too few ($n)"; continue; }
    echo "[$(date +%T)] fit subj$subj n=$n"
    $PYCLONE fit -i $f -o $H5 -c 10 -d beta-binomial -r 5 2>&1 | tail -2
    $PYCLONE write-results-file -i $H5 -o $RES 2>&1 | tail -2
done
echo "=== done ==="
ls $OUT/results_subj*.tsv 2>/dev/null | wc -l
