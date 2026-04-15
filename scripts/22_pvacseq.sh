#!/bin/bash
# pVACseq neoantigen prediction: VEP annotate matched pass VCFs + HLA per subject + MHCflurry binding
set -e

RNA=/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin
VEP=/home/soon/miniconda3/envs/vep/bin/vep
PVACSEQ=$RNA/pvacseq
VEP_CACHE=/mnt/sda1/data/TNT/refs/vep_cache
FASTA=/mnt/sda1/data/TNT/refs/hg38_nochr/genome.fa

IN=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered
OUT=/mnt/sda1/data/TNT/analysis/03_hla/neoantigen
mkdir -p $OUT/vep $OUT/predict
WES_INV=/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv
HLA_TSV=/mnt/sda1/data/TNT/analysis/03_hla/hla_class_I_typing.tsv

# Step 1: VEP annotation (need no-chr matching FASTA)
annotate() {
    local s=$1  # e.g. 1-PR
    local V=$IN/${s}.pass.vcf.gz
    local OUTV=$OUT/vep/${s}.vep.vcf
    [ -s "${OUTV}.gz" ] && return 0
    [ ! -s "$V" ] && return 1
    $VEP --input_file $V --format vcf --output_file $OUTV --vcf \
         --cache --dir_cache $VEP_CACHE --cache_version 115 \
         --assembly GRCh38 --species homo_sapiens --offline \
         --fasta $FASTA \
         --plugin Frameshift --plugin Wildtype \
         --symbol --terms SO --tsl --hgvs --transcript_version \
         --fork 2 2>&1 | tail -2
    /home/soon/miniconda3/envs/wes_somatic/bin/bgzip -f $OUTV
    /home/soon/miniconda3/envs/wes_somatic/bin/tabix -p vcf ${OUTV}.gz
    echo "vep $s"
}
export -f annotate
export IN OUT VEP VEP_CACHE FASTA

# Only matched tumors (unmatched unreliable for neoantigen)
matched_tumors=$(awk -F'\t' 'NR>1 && $3!="normal" && $2!~/^(13|15|16|17|18|19|33)$/ {print $1}' $WES_INV)

for s in $matched_tumors; do annotate $s; done

echo "=== VEP annotation done ==="
ls $OUT/vep/*.vep.vcf.gz 2>/dev/null | wc -l

# Step 2: pvacseq run per tumor with subject HLA
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
    if (outd/'MHC_Class_I'/f'{sid}.all_epitopes.tsv').exists(): continue
    outd.mkdir(parents=True, exist_ok=True)
    # Format HLA: comma-separated
    hla_sub = hla[hla.subject_id==subj]
    if len(hla_sub)==0: continue
    hs = hla_sub.iloc[0]
    alleles = []
    for c in ['A1','A2','B1','B2','C1','C2']:
        a = hs[c]
        if pd.isna(a): continue
        # Strip 3-field to 2-field: HLA-A*24:02:01:01 → HLA-A*24:02
        parts = str(a).split(':')
        if len(parts)>=2:
            alleles.append(':'.join(parts[:2]))
    alleles = list(dict.fromkeys(alleles))  # unique preserve order
    if not alleles: continue
    allele_str = ','.join(alleles)
    # Sample name in VCF (tumor)
    proc = sp.run(['/home/soon/miniconda3/envs/wes_somatic/bin/bcftools','query','-l',str(vep_vcf)], capture_output=True, text=True)
    samples = proc.stdout.strip().split('\n')
    tumor_sm = next((s for s in samples if not s.endswith('-N') and not s.endswith('-N_DNA')), samples[0])
    cmd = ['/home/soon/miniconda3/envs/rnaseq_arabidopsis/bin/pvacseq','run',
           str(vep_vcf), tumor_sm, allele_str,
           'MHCflurry',  # predictor
           str(outd),
           '-e1','8,9,10,11',
           '--iedb-install-directory','/home/soon/miniconda3/envs/rnaseq_arabidopsis',
           '-t','2',
           '--pass-only']
    print(f'[{sid}] HLA={allele_str} → running pvacseq')
    r = sp.run(cmd, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        print(f'  FAILED: {r.stderr[-500:]}')
    else:
        print(f'  OK')
PYEOF
echo "=== pvacseq done ==="
