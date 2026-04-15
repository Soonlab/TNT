#!/bin/bash
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate wes_somatic

IN=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/filtered
OUT=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/annotated
TSV=/mnt/sda1/data/TNT/analysis/02_wes_mutect2/variant_tables
LOG=/mnt/sda1/data/TNT/analysis/logs
mkdir -p $OUT $TSV
DATA=/mnt/sda1/data/TNT/refs/snpeff_data
export _JAVA_OPTIONS="-Xmx4g"

# Annotate sequentially but use 4 parallel, each getting own DB via symlinked dataDir (DB already present)
annotate_one() {
  local v=$1
  local s=$(basename $v .pass.vcf.gz)
  local A=$OUT/${s}.snpeff.vcf
  # Re-annotate if .gz is empty or missing
  local size=0; [ -s "${A}.gz" ] && size=$(stat -c%s "${A}.gz")
  if [ "$size" -lt 5000 ]; then
    echo "[$(date +%T)] annotate $s"
    snpEff -dataDir $DATA -noStats -nodownload GRCh38.99 $v > $A 2> $LOG/snpeff_${s}.log
    bgzip -f $A
    tabix -p vcf ${A}.gz
  fi
  # Extract fields
  SnpSift extractFields ${A}.gz CHROM POS REF ALT "ANN[0].EFFECT" "ANN[0].IMPACT" "ANN[0].GENE" "ANN[0].HGVS_P" "ANN[0].FEATUREID" FILTER "GEN[0].AD" "GEN[0].DP" "GEN[0].AF" > $TSV/${s}.tsv 2> $LOG/snpsift_${s}.log
  local nl=$(wc -l < $TSV/${s}.tsv)
  echo "  $s: $nl lines"
}
export -f annotate_one
export IN OUT TSV LOG DATA

# Run sequentially (one java process at a time) to avoid DB locking issues
for v in $IN/*.pass.vcf.gz; do annotate_one "$v"; done
echo "=== done ==="
ls $TSV | wc -l
wc -l $TSV/*.tsv | tail -5
