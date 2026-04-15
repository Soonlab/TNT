#!/bin/bash
# MSI calling pipeline with msisensor-pro
# 1) scan reference for microsatellites
# 2) T-N calling for 41 matched tumors (msi mode)
# 3) Build baseline from 28 normals + tumor-only 'pro' mode for 8 unmatched
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate wes_somatic

REF=/mnt/sda1/data/TNT/refs/hg38_nochr/genome.fa
BED=/mnt/sda1/data/TNT/refs/hg38_nochr/exome_targets.bed
WES=/mnt/sda1/data/TNT/TNT_WES
OUT=/mnt/sda1/data/TNT/analysis/02_wes_tmb_msi/msi
MS=$OUT/microsatellites.list
BAS=$OUT/baseline
mkdir -p $OUT/paired $OUT/tumor_only $OUT/baseline_tmp $BAS

# Step 1: scan
if [ ! -s $MS ]; then
  echo "[$(date +%T)] scan microsatellites"
  msisensor-pro scan -d $REF -o $MS -b 8 2>&1 | tail -3
fi
wc -l $MS

bam() { echo "$WES/${1}_DNA/${1}_DNA.recal.bam"; }

# Step 2: T-N for 41 matched
msi_tn() {
  local tumor=$1
  local subj=$(echo $tumor | cut -d- -f1)
  local NBAM=$(bam "${subj}-N")
  local TBAM=$(bam $tumor)
  local O=$OUT/paired/${tumor}
  [ -s "${O}" ] && return 0
  [ ! -s "$NBAM" ] && return 1
  msisensor-pro msi -d $MS -n $NBAM -t $TBAM -e $BED -o $O -b 4 2> $OUT/paired/${tumor}.log
  echo "tn $tumor done"
}
export -f msi_tn bam
export MS OUT BED WES

awk -F'\t' 'NR>1 && $3!="normal" {print $1}' /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | \
  parallel -j 12 --line-buffer msi_tn {}

# Step 3: Build baseline from normals
# msisensor-pro baseline -d <ms_list> -i <config> -o <baseline_dir>
if [ ! -s $BAS/baseline.list ]; then
  # config: tab-separated case_name, normal_bam
  CFG=$OUT/baseline_config.tsv
  : > $CFG
  awk -F'\t' 'NR>1 && $3=="normal" {print $1}' /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | while read s; do
    B=$(bam $s); [ -s "$B" ] && printf "%s\t%s\n" "$s" "$B" >> $CFG
  done
  wc -l $CFG
  msisensor-pro baseline -d $MS -i $CFG -o $BAS -b 8 -e $BED 2>&1 | tail -3
fi
ls $BAS | head

# Step 4: tumor-only 'pro' for 8 unmatched
msi_to() {
  local tumor=$1
  local TBAM=$(bam $tumor)
  local O=$OUT/tumor_only/${tumor}
  [ -s "${O}" ] && return 0
  msisensor-pro pro -d $BAS/Homo_sapiens_assembly38_baseline -t $TBAM -e $BED -o $O -b 4 2> $OUT/tumor_only/${tumor}.log
  echo "to $tumor done"
}
export -f msi_to
export BAS

for t in 13-PR 13-PO 15-P 16-P 17-P 18-P 19-P 33-P; do
  msi_to $t &
done
wait
echo "=== MSI pipeline done ==="
ls $OUT/paired/ | head
ls $OUT/tumor_only/ | head
