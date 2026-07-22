#!/bin/bash
# TRUST4 TCR/BCR repertoire reconstruction from RNA-seq BAMs
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate immune_tcr

RNA_BAM=/mnt/sda1/data/TNT/TNT_RNAseq/BAM_files
OUT=/mnt/sda1/data/TNT/analysis/06_rna_immune/trust4
mkdir -p $OUT
REF=/mnt/sda1/data/TNT/refs/trust4/hg38_bcrtcr.nochr.fa
IMGT=/mnt/sda1/data/TNT/refs/trust4/human_IMGT+C.fa

run_one() {
  local s=$1
  local B=$RNA_BAM/${s}/${s}_sorted.bam
  local OD=$OUT/${s}
  [ -s "$OD/${s}_report.tsv" ] && return 0
  mkdir -p $OD
  cd $OD
  run-trust4 -b $B -f $REF --ref $IMGT -o ${s} -t 4 2>&1 | tail -2
}
export -f run_one
export RNA_BAM OUT REF IMGT

ls $RNA_BAM | parallel -j 8 --line-buffer run_one {}
echo "=== TRUST4 done ==="
ls $OUT | wc -l
