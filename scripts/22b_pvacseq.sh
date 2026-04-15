#!/bin/bash
set -e
source /home/soon/miniconda3/etc/profile.d/conda.sh

VEP_CACHE=/mnt/sda1/data/TNT/refs/vep_cache
FASTA=/mnt/sda1/data/TNT/refs/hg38_nochr/genome.fa
IN=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered
OUT=/mnt/sda1/data/TNT/analysis/03_hla/neoantigen
mkdir -p $OUT/vep $OUT/predict

WES_INV=/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv

# Step 1: VEP annotation (via conda activate within subshell)
conda activate vep
matched_tumors=$(awk -F'\t' 'NR>1 && $3!="normal" && $2!~/^(13|15|16|17|18|19|33)$/ {print $1}' $WES_INV)

annotate() {
    local s=$1
    local V=$IN/${s}.pass.vcf.gz
    local OUTV=$OUT/vep/${s}.vep.vcf
    [ -s "${OUTV}.gz" ] && return 0
    vep --input_file $V --format vcf --output_file $OUTV --vcf \
        --cache --dir_cache $VEP_CACHE --cache_version 115 \
        --assembly GRCh38 --species homo_sapiens --offline \
        --fasta $FASTA --symbol --terms SO --tsl --hgvs --transcript_version --fork 2 \
        2> $OUT/vep/${s}.log
    /home/soon/miniconda3/envs/wes_somatic/bin/bgzip -f $OUTV
    /home/soon/miniconda3/envs/wes_somatic/bin/tabix -p vcf ${OUTV}.gz
    echo "annotated $s ($(zcat ${OUTV}.gz | grep -v ^# | wc -l) variants)"
}
export -f annotate
export IN OUT VEP_CACHE FASTA

# Run in parallel, 4 at a time
echo "$matched_tumors" | tr ' ' '\n' | /home/soon/miniconda3/envs/wes_somatic/bin/parallel -j 8 --line-buffer annotate {}

echo "=== VEP done: $(ls $OUT/vep/*.vep.vcf.gz 2>/dev/null | wc -l) VCFs ==="

# Step 2: pvacseq
conda deactivate
conda activate rnaseq_arabidopsis

/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin/python3 << 'PYEOF'
import pandas as pd, subprocess as sp, os
from pathlib import Path

OUT = Path('/mnt/sda1/data/TNT/analysis/03_hla/neoantigen')
hla = pd.read_csv('/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv', sep='\t')
inv = pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')

tumors = inv[(inv.timepoint!='normal') & (~inv.subject_id.isin([13,15,16,17,18,19,33]))]
for _, r in tumors.iterrows():
    sid = r['sample_id']; subj = r['subject_id']
    vep_vcf = OUT/'vep'/f'{sid}.vep.vcf.gz'
    if not vep_vcf.exists(): continue
    outd = OUT/'predict'/sid
    if (outd/'MHC_Class_I'/f'{sid}.all_epitopes.aggregated.tsv').exists(): continue
    outd.mkdir(parents=True, exist_ok=True)
    hla_sub = hla[hla.subject_id==subj]
    if len(hla_sub)==0: continue
    hs = hla_sub.iloc[0]
    alleles = []
    for c in ['A1','A2','B1','B2','C1','C2']:
        a = hs[c]
        if pd.isna(a): continue
        parts = str(a).split(':')
        if len(parts)>=2:
            alleles.append(':'.join(parts[:2]))
    alleles = list(dict.fromkeys(alleles))
    if not alleles: continue
    allele_str = ','.join(alleles)
    proc = sp.run(['/home/soon/miniconda3/envs/wes_somatic/bin/bcftools','query','-l',str(vep_vcf)],
                  capture_output=True, text=True)
    samples = proc.stdout.strip().split('\n')
    tumor_sm = next((s for s in samples if not s.endswith('-N') and not s.endswith('-N_DNA')), samples[0])
    cmd = ['/home/soon/miniconda3/envs/pvactools/bin/pvacseq','run',
           str(vep_vcf), tumor_sm, allele_str,
           'MHCflurry',
           str(outd),
           '-e1','8,9,10,11',
           '-t','2',
           '--pass-only',
           '--iedb-install-directory','/home/soon/miniconda3/envs/pvactools']
    print(f'[{sid}] HLA={allele_str}', flush=True)
    r = sp.run(cmd, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        print(f'  FAIL: {r.stderr[-400:]}', flush=True)
    else:
        print(f'  OK', flush=True)
PYEOF
echo "=== pvacseq done ==="
