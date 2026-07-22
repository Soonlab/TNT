#!/bin/bash
# HLA class I typing using OptiType on normal WES BAMs (28 samples)
# Extracts reads mapping to HLA class I region (chr6:28-34Mb), converts to FASTQ,
# OptiType maps via razers3 against IMGT/HLA reference.
set -euo pipefail
source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate hla_typing

WES=/mnt/sda1/data/TNT/TNT_WES
OUT=/mnt/sda1/data/TNT/analysis/03_hla/optitype
TMP=/mnt/sda1/data/TNT/analysis/03_hla/tmp
mkdir -p $OUT $TMP

# HLA class I region on GRCh38 (no chr prefix) — MHC region
# The full MHC region is 6:28,510,120-33,480,577 (hg38)
REGION="6:28510120-33480577"

OPTI_DIR=$(dirname $(dirname $(which OptiTypePipeline.py)))
REF_HLA=$(find $OPTI_DIR -name "hla_reference_dna.fasta" 2>/dev/null | head -1)
CFG_TMPL=$(find $OPTI_DIR -name "config.ini" 2>/dev/null | head -1)
if [ -z "$REF_HLA" ]; then
  echo "HLA reference not found; check OptiType install"; exit 1
fi
echo "HLA ref: $REF_HLA"
echo "Config: $CFG_TMPL"

# Prepare config
CFG=$OUT/config.ini
if [ ! -s $CFG ]; then
  cp $CFG_TMPL $CFG
  # adjust paths
  sed -i 's|^razers3=.*|razers3='"$(which razers3)"'|' $CFG
  sed -i 's|^solver=.*|solver=cbc|' $CFG || true
  sed -i 's|^threads=.*|threads=4|' $CFG || true
fi

type_sample() {
  local s=$1  # e.g. 1-N
  local BAM=$WES/${s}_DNA/${s}_DNA.recal.bam
  local OD=$OUT/${s}
  [ -s "$OD/${s}_result.tsv" ] && { echo "skip $s"; return 0; }
  mkdir -p $OD
  local PFX=$TMP/${s}
  echo "[$(date +%T)] HLA $s"
  # Extract reads in MHC region + unmapped-with-mapped-mate (both useful)
  samtools view -b $BAM 6:28510120-33480577 -o $PFX.mhc.bam
  samtools sort -n -@ 4 $PFX.mhc.bam -o $PFX.nsort.bam
  samtools fastq -@ 4 -1 $PFX.1.fq -2 $PFX.2.fq -0 /dev/null -s /dev/null -N $PFX.nsort.bam 2>/dev/null
  local N=$(wc -l < $PFX.1.fq)
  if [ "$N" -lt 40 ]; then echo "  too few reads for $s ($N)"; return 1; fi
  OptiTypePipeline.py -i $PFX.1.fq $PFX.2.fq --dna -o $OD -c $CFG -v 2>&1 | tail -5
  # Rename result
  mv $OD/*/*_result.tsv $OD/${s}_result.tsv 2>/dev/null || true
  mv $OD/*/*_coverage_plot.pdf $OD/${s}_coverage_plot.pdf 2>/dev/null || true
  rm -f $PFX.mhc.bam $PFX.nsort.bam $PFX.1.fq $PFX.2.fq
  echo "  done $s"
}
export -f type_sample
export WES OUT TMP CFG

# HLA typing from N (normal) samples where available, otherwise from tumor (just for patients w/o normal)
awk -F'\t' 'NR>1 && $3=="normal" {print $1}' /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | \
  parallel -j 8 --line-buffer type_sample {}

# For subjects without normal (13,15-19,33), use tumor sample
for t in 13-PR 15-P 16-P 17-P 18-P 19-P 33-P; do
  type_sample $t || true
done
echo "=== HLA typing done ==="
