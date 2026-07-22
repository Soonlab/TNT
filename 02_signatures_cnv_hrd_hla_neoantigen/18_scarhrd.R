# scarHRD genomic scar score (LOH + LST + TAI)
# Input: CNVkit .call.cns converted to sequenza-style format
# Output: per-sample HRD score + sum
suppressPackageStartupMessages({
  library(scarHRD); library(dplyr); library(readr)
})
IN <- '/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/cnvkit'
OUT <- '/mnt/sda1/data/TNT/analysis/04_wes_cnv_clonal/hrd'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# scarHRD expects a sequenza-like input: SampleID, Chromosome, Start_position, End_position, total_cn, A_cn, B_cn
# From CNVkit: cn (total), depth, log2; allele-specific requires BAF but CNVkit doesn't provide natively.
# Approximate: assume heterozygous split cn/2 floor and ceil; if cn==1 treat as 1,0; if cn==2 as 1,1; if >=3 assume 2,cn-2
files <- list.files(IN, pattern='\\.call\\.cns$', full.names=TRUE)
cat('Processing', length(files), 'samples\n')
res <- data.frame()
for (f in files) {
  sid <- sub('_DNA\\.call\\.cns$','', basename(f))
  d <- read.table(f, header=TRUE, sep='\t', stringsAsFactors=FALSE)
  d$SampleID <- sid
  d$Chromosome <- d$chromosome
  d$Start_position <- d$start
  d$End_position <- d$end
  d$total_cn <- d$cn
  # Naive A/B allele split (major/minor)
  d$A_cn <- pmax(d$cn - pmin(1, d$cn), 0) # major
  d$B_cn <- pmin(1, d$cn)  # minor (0 if cn=0, 1 if cn>=1)
  # Proper: if cn==0 → 0,0; cn==1 → 1,0; cn==2 → 1,1; cn>=3 → cn-1,1 (assume minor remains 1)
  d$A_cn <- ifelse(d$cn==0, 0, ifelse(d$cn==1, 1, d$cn-1))
  d$B_cn <- ifelse(d$cn<=1, 0, 1)
  d2 <- d[, c('SampleID','Chromosome','Start_position','End_position','total_cn','A_cn','B_cn')]
  tmp <- tempfile(fileext='.txt')
  write.table(d2, tmp, sep='\t', quote=FALSE, row.names=FALSE)
  out <- tryCatch(scar_score(tmp, reference='grch38', seqz=FALSE), error=function(e) NA)
  if (is.null(out) || all(is.na(out))) {
    cat(' ', sid, ': FAILED\n'); next
  }
  rr <- as.data.frame(t(out)); rr$SampleID <- sid
  res <- rbind(res, rr)
  cat(' ', sid, ': HRD sum =', ifelse(!is.null(out[['HRD-sum']]), out[['HRD-sum']], 'NA'), '\n')
}
write_tsv(res, file.path(OUT, 'hrd_scores.tsv'))
cat('Saved', file.path(OUT,'hrd_scores.tsv'), '\n')

# Response association
inv <- readr::read_tsv('/mnt/sda1/data/TNT/analysis/00_cohort/wes_inventory.tsv', show_col_types=FALSE)
m <- merge(res, inv[,c('sample_id','subject_id','timepoint','response_bin')],
           by.x='SampleID', by.y='sample_id')
cat('\n=== HRD good vs bad (pre, matched) ===\n')
pre <- m[m$timepoint=='pre' & !m$subject_id %in% c(13,15,16,17,18,19,33),]
if (nrow(pre)>=8) {
  g <- pre[pre$response_bin=='good','HRD-sum']; b <- pre[pre$response_bin=='bad','HRD-sum']
  cat('good n=', length(g), ' median=', median(g, na.rm=TRUE),
      ' | bad n=', length(b), ' median=', median(b, na.rm=TRUE),
      ' | p=', wilcox.test(g,b)$p.value, '\n', sep='')
}
write_tsv(m, file.path(OUT, 'hrd_with_meta.tsv'))
