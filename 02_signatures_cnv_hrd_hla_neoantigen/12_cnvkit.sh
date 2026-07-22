#!/bin/bash
# CNV analysis using CNVkit batch mode
# Uses 28 normals to build flat reference, then runs per-tumor
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_arabidopsis  # has cnvkit

REF=/mnt/sda1/data/TNT/refs/hg38_nochr/genome.fa
BED=/mnt/sda1/data/TNT/refs/hg38_nochr/exome_targets.bed
WES=/mnt/sda1/data/TNT/TNT_WES
OUT=/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/cnvkit
mkdir -p $OUT
INV=/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv

bam() { echo "$WES/${1}_DNA/${1}_DNA.recal.bam"; }

# Build BAM lists
tumors=$(awk -F'\t' 'NR>1 && $3!="normal" {print $1}' $INV)
normals=$(awk -F'\t' 'NR>1 && $3=="normal" {print $1}' $INV)
tumor_bams=""
for s in $tumors; do tumor_bams="$tumor_bams $(bam $s)"; done
normal_bams=""
for s in $normals; do normal_bams="$normal_bams $(bam $s)"; done

cd $OUT

# batch: builds reference + runs CN calling on tumors
cnvkit.py batch $tumor_bams --normal $normal_bams \
  --targets $BED \
  --fasta $REF \
  --access /home/soon/miniconda3/envs/rnaseq_arabidopsis/share/cnvkit-*/data/access-5k-mappable.hg38.bed 2>/dev/null || \
cnvkit.py batch $tumor_bams --normal $normal_bams \
  --targets $BED \
  --fasta $REF \
  --output-reference $OUT/reference.cnn \
  --output-dir $OUT \
  --processes 16 2>&1 | tail -10
echo "=== CNVkit done ==="
ls $OUT/*.cns 2>/dev/null | wc -l
