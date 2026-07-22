#!/bin/bash
# Annotate 49 PASS VCFs with snpEff + extract to table
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate wes_somatic

IN=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered
OUT=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/annotated
TSV=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/variant_tables
LOG=/mnt/sda1/data/TNT/analysis/logs
mkdir -p $OUT $TSV $LOG
DATA=/mnt/sda1/data/TNT/refs/snpeff_data

export _JAVA_OPTIONS="-Xmx6g"

annotate() {
  local v=$1
  local s=$(basename $v .pass.vcf.gz)
  local A=$OUT/${s}.snpeff.vcf
  [ -s "${A}.gz" ] && return 0
  snpEff -dataDir $DATA -noStats GRCh38.99 $v > $A 2> $LOG/snpeff_${s}.log
  bgzip -f $A
  tabix -p vcf ${A}.gz
  # Extract table: CHROM POS REF ALT TYPE GENE EFFECT IMPACT HGVS_p VAF DP AF FILTER
  SnpSift extractFields ${A}.gz CHROM POS REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].HGVS_P" "ANN[0].FEATUREID" FILTER "GEN[0].AD" "GEN[0].DP" "GEN[0].AF" > $TSV/${s}.tsv 2> $LOG/snpsift_${s}.log
  echo "annotated $s"
}
export -f annotate
export IN OUT TSV LOG DATA

ls $IN/*.pass.vcf.gz | parallel -j 12 --line-buffer annotate {}
echo done.
