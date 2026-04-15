#!/bin/bash
# HLA LOH lite: per-subject HLA ref + realign MHC reads + allele count imbalance
# Skips novoalign; uses bwa-mem against IMGT HLA fasta
set -e
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate hla_typing

WES=/mnt/sda1/data/TNT/TNT_WES
HLA_REF=/home/soon/miniconda3/envs/hla_typing/bin/data/hla_reference_dna.fasta
OUT=/mnt/sda1/data/TNT/analysis/03_hla/loh_lite
mkdir -p $OUT/subject_refs $OUT/realign $OUT/counts
TMP=$OUT/tmp; mkdir -p $TMP

HLA_TSV=/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv

# Build alt→ID map
/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin/python3 << 'PYEOF'
import pandas as pd, os, re
h = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv', sep='\t')
# Parse IMGT ID mapping from hla reference fasta
mapping = {}
with open('/home/soon/miniconda3/envs/hla_typing/bin/data/hla_reference_dna.fasta') as f:
    for line in f:
        if line.startswith('>'):
            m = re.match(r'>(\S+)\s+(HLA-\S+)', line)
            if m:
                imgt, allele = m.group(1), m.group(2)
                if '_' not in imgt:
                    mapping[allele] = imgt
# For each subject, build list of 6 IMGT IDs (A1, A2, B1, B2, C1, C2)
for _, r in h.iterrows():
    subj = int(r['subject_id'])
    ids = []
    for col in ['A1','A2','B1','B2','C1','C2']:
        val = r[col]  # e.g. HLA-A*24:02:01:01
        imgt = mapping.get(val, None)
        # Try prefix match if exact fails
        if imgt is None:
            prefix = val
            while ':' in prefix and imgt is None:
                prefix = prefix.rsplit(':',1)[0]
                imgt = mapping.get(prefix, None)
        ids.append((col, val, imgt))
    # Save per-subject list
    with open(f'/mnt/sda1/data/TNT/analysis/03_hla/loh_lite/subject_refs/subj{subj}_alleles.tsv','w') as f:
        f.write('locus\tallele\timgt_id\n')
        for c,v,i in ids:
            f.write(f'{c}\t{v}\t{i}\n')
PYEOF

echo "Built per-subject allele maps"

# Build per-subject HLA fasta
build_ref() {
  local subj=$1
  local AL=$OUT/subject_refs/subj${subj}_alleles.tsv
  local REF=$OUT/subject_refs/subj${subj}.fa
  [ -s "$REF" ] && return 0
  [ ! -s "$AL" ] && return 1
  # Extract unique IMGT IDs from tsv
  ids=$(tail -n +2 $AL | awk -F'\t' '$3!="None" && $3!="" {print $3}' | sort -u)
  if [ -z "$ids" ]; then echo "no IDs for subj$subj"; return 1; fi
  samtools faidx $HLA_REF $ids > $REF 2>/dev/null
  # Create simplified headers (just IMGT ID)
  awk '/^>/{split($1,a,"[>]"); print ">"a[2]; next}{print}' $REF > ${REF}.tmp && mv ${REF}.tmp $REF
  samtools faidx $REF
  bwa index $REF 2>/dev/null
}
export -f build_ref
export OUT HLA_REF

# Realign MHC reads of tumor + normal to subject HLA fasta + count
process_sample() {
  local subj=$1
  local sid=$2   # e.g. 1-N or 1-PR
  local REF=$OUT/subject_refs/subj${subj}.fa
  local BAM=$WES/${sid}_DNA/${sid}_DNA.recal.bam
  local BAMOUT=$OUT/realign/${sid}.bam
  local CNT=$OUT/counts/${sid}.counts.tsv
  [ -s "$CNT" ] && return 0
  [ ! -s "$REF" ] && { echo "no ref for $sid"; return 1; }
  # Extract reads from MHC region (6:28-33Mb), collate, bwa-mem align
  samtools view -b $BAM 6:28510120-33480577 | \
    samtools collate -O -@ 2 - | \
    samtools fastq -1 $TMP/${sid}_1.fq -2 $TMP/${sid}_2.fq -0 /dev/null -s /dev/null -N - 2>/dev/null
  bwa mem -t 4 $REF $TMP/${sid}_1.fq $TMP/${sid}_2.fq 2>/dev/null | \
    samtools sort -@ 2 -o $BAMOUT -
  samtools index $BAMOUT
  # Count mapped reads per allele (primary, MQ>=30)
  samtools view -F 2308 -q 30 $BAMOUT | awk '{print $3}' | sort | uniq -c | \
    awk 'BEGIN{OFS="\t"}{print $2,$1}' > $CNT
  rm -f $TMP/${sid}_1.fq $TMP/${sid}_2.fq
}
export -f process_sample
export WES TMP

# Iterate matched-normal subjects
MATCHED_SUBJ=$(tail -n +2 /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | awk -F'\t' '$3=="normal" {print $2}' | sort -u)

for subj in $MATCHED_SUBJ; do
  build_ref $subj
done

# Process N + tumors in parallel
for subj in $MATCHED_SUBJ; do
  for sid in ${subj}-N ${subj}-PR ${subj}-PO ${subj}-P; do
    BAM=$WES/${sid}_DNA/${sid}_DNA.recal.bam
    [ -s "$BAM" ] || continue
    echo "$subj $sid"
  done
done | parallel -j 12 --colsep ' ' --line-buffer process_sample {1} {2}

echo "=== Done realign + count ==="
ls $OUT/counts/ | wc -l

# Summarize LOH per subject × tumor timepoint
/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin/python3 << 'PYEOF'
import pandas as pd, numpy as np, os
from pathlib import Path
from scipy import stats

OUT = Path('/mnt/sda1/data/TNT/analysis/03_hla/loh_lite')
inv = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')
hla = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv', sep='\t')

results=[]
for subj, sub_hla in hla.groupby('subject_id'):
    if not (OUT/'subject_refs'/f'subj{subj}_alleles.tsv').exists(): continue
    alleles = pd.read_csv(OUT/'subject_refs'/f'subj{subj}_alleles.tsv', sep='\t')
    # Normal counts
    ncnt_f = OUT/'counts'/f'{subj}-N.counts.tsv'
    if not ncnt_f.exists(): continue
    ncnt = dict([l.strip().split('\t')[:2] for l in open(ncnt_f) if l.strip()])
    # Each tumor timepoint
    for tp_sid in [f'{subj}-PR', f'{subj}-PO', f'{subj}-P']:
        tcnt_f = OUT/'counts'/f'{tp_sid}.counts.tsv'
        if not tcnt_f.exists(): continue
        tcnt = dict([l.strip().split('\t')[:2] for l in open(tcnt_f) if l.strip()])
        # Per-locus: compare allele1 vs allele2 ratio in tumor vs normal
        for locus in ['A','B','C']:
            a1 = alleles[alleles.locus==f'{locus}1']['imgt_id'].values[0]
            a2 = alleles[alleles.locus==f'{locus}2']['imgt_id'].values[0]
            if a1==a2: continue  # homozygous
            n1 = int(ncnt.get(str(a1),0)); n2 = int(ncnt.get(str(a2),0))
            t1 = int(tcnt.get(str(a1),0)); t2 = int(tcnt.get(str(a2),0))
            if n1+n2 < 20 or t1+t2 < 20: continue
            # Tumor allele ratio: t1/(t1+t2). Normal ratio expected ~0.5.
            t_ratio = t1/(t1+t2)
            n_ratio = n1/(n1+n2)
            # Test: tumor significantly different from expected 0.5 AND from normal
            # Use Fisher exact: tumor vs normal
            tab = [[t1, t2],[n1, n2]]
            try:
                _, p_tn = stats.fisher_exact(tab)
            except: p_tn = np.nan
            # LOH call: |tumor - 0.5| > 0.15 AND p<0.05
            loh = (abs(t_ratio - 0.5) > 0.15) and (p_tn < 0.05)
            results.append({'subject_id':subj,'sample':tp_sid,'locus':f'HLA-{locus}',
                           'allele1':a1,'allele2':a2,
                           'normal_c1':n1,'normal_c2':n2,'tumor_c1':t1,'tumor_c2':t2,
                           'normal_ratio':n_ratio,'tumor_ratio':t_ratio,
                           'fisher_p':p_tn,'LOH_call':loh})

df = pd.DataFrame(results)
df.to_csv(OUT/'hla_loh_lite_results.tsv', sep='\t', index=False)
print(f'{len(df)} locus-tumor tests')
print(f'LOH events called: {df.LOH_call.sum()}')

# Per-subject LOH status
loh_per = df.groupby(['subject_id','sample'])['LOH_call'].sum().reset_index()
loh_per = loh_per.merge(inv[['sample_id','response_bin']].rename(columns={'sample_id':'sample'}), on='sample', how='left')
print('\n=== LOH by sample ===')
print(loh_per.sort_values('LOH_call', ascending=False).head(20).to_string(index=False))

# Response association
pre_loh = df[df['sample'].str.endswith('-PR') | df['sample'].str.endswith('-P')]
pre_any = pre_loh.groupby('subject_id')['LOH_call'].max().reset_index()
clin = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/clinical_master.tsv', sep='\t')
pre_any = pre_any.merge(clin[['subject_id','response_bin']], on='subject_id')
tab = pd.crosstab(pre_any.response_bin, pre_any.LOH_call.astype(bool))
print('\n=== HLA LOH pre × response ===')
print(tab)
if tab.shape == (2,2):
    _, p = stats.fisher_exact(tab.values)
    print(f'Fisher p={p:.3f}')
PYEOF
