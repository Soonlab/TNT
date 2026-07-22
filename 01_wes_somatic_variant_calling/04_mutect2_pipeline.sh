#!/bin/bash
# TNT Mutect2 somatic variant calling pipeline
# Runs: (1) PoN build from 28 N samples, (2) T-N calling for matched tumors, (3) tumor-only+PoN for unmatched
set -euo pipefail

source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate wes_somatic

REF=/mnt/sda1/data/TNT/refs/hg38_nochr/genome.fa
GNOMAD=/mnt/sda1/data/TNT/refs/hg38_nochr/af-only-gnomad.vcf.gz
SMALL_EXAC=/mnt/sda1/data/TNT/refs/hg38_nochr/small_exac_common_3.vcf.gz
INTERVALS=/mnt/sda1/data/TNT/refs/hg38_nochr/exome_targets.bed
WES=/mnt/sda1/data/TNT/TNT_WES
OUT=/mnt/sda1/data/TNT/analysis/02_wes_mutect2
LOG=$OUT/logs
mkdir -p $OUT/{pon_tumoronly,pon,tn_calls,tumor_only_calls,filtered,contamination} $LOG

# Parallel workers
N_JOBS=${N_JOBS:-12}
# Per-process memory
MEM="-Xmx6g"

bam() { echo "$WES/${1}_DNA/${1}_DNA.recal.bam"; }
sn() {
  # extract SM from BAM header (tolerant to missing SM) — fallback to sample_id
  local b="$1"
  samtools view -H "$b" 2>/dev/null | awk -F'\t' '/^@RG/{for(i=1;i<=NF;i++) if($i ~ /^SM:/){sub(/^SM:/,"",$i); print $i; exit}}' | head -1
}

# --- Stage 1: PoN tumor-only runs on 28 normals ---
pon_call() {
  local sid=$1  # e.g. 1-N
  local B=$(bam $sid)
  local O=$OUT/pon_tumoronly/${sid}.vcf.gz
  [ -s "$O" ] && { echo "skip $sid"; return 0; }
  local smsm=$(sn $B); [ -z "$smsm" ] && smsm=$sid
  echo "[$(date +%T)] PoN call: $sid (SM=$smsm)"
  gatk --java-options "$MEM" Mutect2 \
    -R $REF -I $B -tumor "$smsm" \
    --germline-resource $GNOMAD \
    -L $INTERVALS --interval-padding 100 \
    --max-mnp-distance 0 \
    -O $O 2> $LOG/pon_${sid}.log
}
export -f pon_call bam sn
export WES REF GNOMAD INTERVALS OUT LOG MEM

stage_pon_runs() {
  echo "=== Stage 1: PoN tumor-only calls ==="
  awk -F'\t' 'NR>1 && $3=="normal" {print $1}' /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | \
    parallel -j $N_JOBS --line-buffer pon_call {}
}

stage_pon_build() {
  echo "=== Stage 2: CreateSomaticPanelOfNormals ==="
  local ARGS=""
  for v in $OUT/pon_tumoronly/*.vcf.gz; do ARGS+=" -V $v"; done
  gatk --java-options "-Xmx32g" GenomicsDBImport \
    -R $REF -L $INTERVALS --interval-padding 100 \
    --genomicsdb-workspace-path $OUT/pon/gdb \
    --merge-input-intervals --reader-threads 8 \
    $ARGS 2> $LOG/gdb_import.log
  gatk --java-options "-Xmx16g" CreateSomaticPanelOfNormals \
    -R $REF -V gendb://$OUT/pon/gdb \
    --germline-resource $GNOMAD \
    -O $OUT/pon/tnt_cohort_pon.vcf.gz 2> $LOG/pon_build.log
  echo "PoN: $OUT/pon/tnt_cohort_pon.vcf.gz"
}

# --- Stage 3: T-N calling (matched) ---
tn_call() {
  local tumor=$1  # e.g. 1-PR
  local subj=$(echo $tumor | cut -d- -f1)
  local normal="${subj}-N"
  local NBAM=$(bam $normal)
  local TBAM=$(bam $tumor)
  if [ ! -s "$NBAM" ]; then echo "  no-N $tumor: run tumor-only branch"; return 1; fi
  local O=$OUT/tn_calls/${tumor}.vcf.gz
  [ -s "$O" ] && { echo "skip $tumor"; return 0; }
  local tsm=$(sn $TBAM); [ -z "$tsm" ] && tsm=$tumor
  local nsm=$(sn $NBAM); [ -z "$nsm" ] && nsm=$normal
  echo "[$(date +%T)] T-N: $tumor vs $normal"
  gatk --java-options "$MEM" Mutect2 \
    -R $REF -I $TBAM -I $NBAM \
    -tumor "$tsm" -normal "$nsm" \
    --germline-resource $GNOMAD \
    --panel-of-normals $OUT/pon/tnt_cohort_pon.vcf.gz \
    -L $INTERVALS --interval-padding 100 \
    --f1r2-tar-gz $OUT/tn_calls/${tumor}.f1r2.tar.gz \
    -O $O 2> $LOG/tn_${tumor}.log
}

# --- Stage 4: tumor-only + PoN for unmatched ---
to_call() {
  local tumor=$1
  local TBAM=$(bam $tumor)
  local O=$OUT/tumor_only_calls/${tumor}.vcf.gz
  [ -s "$O" ] && { echo "skip $tumor"; return 0; }
  local tsm=$(sn $TBAM); [ -z "$tsm" ] && tsm=$tumor
  echo "[$(date +%T)] TO+PoN: $tumor"
  gatk --java-options "$MEM" Mutect2 \
    -R $REF -I $TBAM -tumor "$tsm" \
    --germline-resource $GNOMAD \
    --panel-of-normals $OUT/pon/tnt_cohort_pon.vcf.gz \
    -L $INTERVALS --interval-padding 100 \
    --f1r2-tar-gz $OUT/tumor_only_calls/${tumor}.f1r2.tar.gz \
    -O $O 2> $LOG/to_${tumor}.log
}

export -f tn_call to_call

stage_calls() {
  echo "=== Stage 3: T-N calling (matched tumors) ==="
  awk -F'\t' 'NR>1 && $3!="normal" {print $1}' /mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv | \
    parallel -j $N_JOBS --line-buffer --joblog $LOG/joblog_tn.txt tn_call {} || true

  echo "=== Stage 4: Tumor-only + PoN (unmatched) ==="
  # Find tumors whose subject has no matched normal, or whose tn_call failed
  python3 -c "
import pandas as pd
inv=pd.read_csv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', sep='\t')
has_normal=set(inv[inv.timepoint=='normal'].subject_id)
tumors=inv[inv.timepoint!='normal']
print('\n'.join(tumors[~tumors.subject_id.isin(has_normal)].sample_id.tolist()))
" | parallel -j $N_JOBS --line-buffer --joblog $LOG/joblog_to.txt to_call {} || true
}

# --- Stage 5: Filter ---
filter_sample() {
  local tumor=$1
  local branch=$2  # tn or to
  local VCF=$OUT/${branch}_calls/${tumor}.vcf.gz
  local F1R2=$OUT/${branch}_calls/${tumor}.f1r2.tar.gz
  local ROM=$OUT/${branch}_calls/${tumor}.rom.tar.gz
  local CONTAM=$OUT/contamination/${tumor}.contam.tsv
  local TBAM=$(bam $tumor)
  local OUTV=$OUT/filtered/${tumor}.filtered.vcf.gz
  local PASS=$OUT/filtered/${tumor}.pass.vcf.gz
  [ -s "$PASS" ] && { echo "skip filter $tumor"; return 0; }
  [ ! -s "$VCF" ] && { echo "no VCF $tumor"; return 1; }

  # LearnReadOrientationModel
  gatk --java-options "$MEM" LearnReadOrientationModel \
    -I $F1R2 -O $ROM 2> $LOG/rom_${tumor}.log || true

  # GetPileupSummaries + CalculateContamination
  gatk --java-options "$MEM" GetPileupSummaries \
    -I $TBAM -V $SMALL_EXAC -L $SMALL_EXAC \
    -O $OUT/contamination/${tumor}.pileup.tsv 2> $LOG/pileup_${tumor}.log || true
  gatk --java-options "$MEM" CalculateContamination \
    -I $OUT/contamination/${tumor}.pileup.tsv \
    -O $CONTAM \
    --tumor-segmentation $OUT/contamination/${tumor}.segments.tsv 2> $LOG/contam_${tumor}.log || true

  CONTARG=""
  [ -s "$CONTAM" ] && CONTARG="--contamination-table $CONTAM --tumor-segmentation $OUT/contamination/${tumor}.segments.tsv"
  ROMARG=""
  [ -s "$ROM" ] && ROMARG="--ob-priors $ROM"

  gatk --java-options "$MEM" FilterMutectCalls \
    -R $REF -V $VCF $CONTARG $ROMARG \
    --filtering-stats $OUT/filtered/${tumor}.stats \
    -O $OUTV 2> $LOG/filter_${tumor}.log

  bcftools view -f PASS $OUTV -Oz -o $PASS
  bcftools index -t $PASS
  echo "filtered $tumor"
}
export -f filter_sample

stage_filter() {
  echo "=== Stage 5: Filter + LROM + Contamination ==="
  (for v in $OUT/tn_calls/*.vcf.gz; do s=$(basename $v .vcf.gz); echo -e "$s\ttn"; done
   for v in $OUT/tumor_only_calls/*.vcf.gz; do s=$(basename $v .vcf.gz); echo -e "$s\tto"; done) | \
    parallel -j $N_JOBS --line-buffer --colsep '\t' --joblog $LOG/joblog_filter.txt filter_sample {1} {2} || true
}

# Main
case "${1:-all}" in
  pon)    stage_pon_runs ;;
  ponbuild) stage_pon_build ;;
  calls)  stage_calls ;;
  filter) stage_filter ;;
  all)    stage_pon_runs && stage_pon_build && stage_calls && stage_filter ;;
  *) echo "usage: $0 [pon|ponbuild|calls|filter|all]"; exit 1 ;;
esac
