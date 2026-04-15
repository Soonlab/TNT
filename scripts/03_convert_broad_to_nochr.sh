#!/bin/bash
# Convert Broad hg38 resources (chr-prefixed) to NCBI/Ensembl style (no chr prefix)
# to match Macrogen's BAM naming convention.
# Ensures MT/chrM name reconciliation.
set -euo pipefail

REF=/mnt/sda1/data/TNT/refs/hg38
OUT=/mnt/sda1/data/TNT/refs/hg38_nochr
mkdir -p "$OUT"

source /home/soon/miniconda3/etc/profile.d/conda.sh
conda activate wes_somatic

# 1. FASTA: strip 'chr' prefix; rename chrM -> MT
if [ ! -s "$OUT/genome.fa" ]; then
  echo "[1/6] Renaming FASTA..."
  awk 'BEGIN{OFS=""} /^>/{
    name=$1; sub(/^>chr/,">",name);
    if(name==">M") name=">MT";
    # keep only primary + MT; drop alt/decoy contigs
    if(name ~ /^>([0-9]+|X|Y|MT)$/) { keep=1; print name }
    else { keep=0 }
    next
  } { if(keep) print }' "$REF/Homo_sapiens_assembly38.fasta" > "$OUT/genome.fa"
  samtools faidx "$OUT/genome.fa"
  gatk CreateSequenceDictionary -R "$OUT/genome.fa" 2>&1 | tail -3
fi

# 2. Build chromosome rename map
cat > "$OUT/chr_rename.txt" <<'EOF'
chr1	1
chr2	2
chr3	3
chr4	4
chr5	5
chr6	6
chr7	7
chr8	8
chr9	9
chr10	10
chr11	11
chr12	12
chr13	13
chr14	14
chr15	15
chr16	16
chr17	17
chr18	18
chr19	19
chr20	20
chr21	21
chr22	22
chrX	X
chrY	Y
chrM	MT
EOF

# 3. Rename gnomAD
if [ ! -s "$OUT/af-only-gnomad.vcf.gz" ]; then
  echo "[3/6] Renaming gnomAD VCF..."
  bcftools annotate --rename-chrs "$OUT/chr_rename.txt" "$REF/af-only-gnomad.hg38.vcf.gz" \
    --threads 8 -Oz -o "$OUT/af-only-gnomad.vcf.gz"
  bcftools index -t "$OUT/af-only-gnomad.vcf.gz"
fi

# 4. Rename 1000g PoN
if [ ! -s "$OUT/1000g_pon.vcf.gz" ]; then
  echo "[4/6] Renaming 1000g PoN..."
  bcftools annotate --rename-chrs "$OUT/chr_rename.txt" "$REF/1000g_pon.hg38.vcf.gz" \
    --threads 8 -Oz -o "$OUT/1000g_pon.vcf.gz"
  bcftools index -t "$OUT/1000g_pon.vcf.gz"
fi

# 5. Rename small_exac_common (for contamination)
if [ ! -s "$OUT/small_exac_common_3.vcf.gz" ]; then
  echo "[5/6] Renaming small_exac_common..."
  bcftools annotate --rename-chrs "$OUT/chr_rename.txt" "$REF/small_exac_common_3.hg38.vcf.gz" \
    --threads 8 -Oz -o "$OUT/small_exac_common_3.vcf.gz"
  bcftools index -t "$OUT/small_exac_common_3.vcf.gz"
fi

# 6. Rename exome interval list -> BED (no chr)
if [ ! -s "$OUT/exome_targets.bed" ]; then
  echo "[6/6] Converting exome interval list to BED (no chr prefix)..."
  awk 'BEGIN{OFS="\t"} /^@/{next} {
    chr=$1; sub(/^chr/,"",chr); if(chr=="M") chr="MT";
    if(chr ~ /^([0-9]+|X|Y|MT)$/) print chr, $2-1, $3
  }' "$REF/exome_calling_regions.v1.interval_list" | \
    sort -k1,1 -k2,2n > "$OUT/exome_targets.bed"
  wc -l "$OUT/exome_targets.bed"
fi

echo "=== Done ==="
ls -lh "$OUT/"
