# CMS subtype classification via CMScaller
suppressPackageStartupMessages({
  library(CMScaller); library(dplyr); library(readr)
})
RNA_X <- '/mnt/sda1/data/TNT/TNT_RNAseq/result_RNAseq/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.txt'
META  <- '/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv'
OUT   <- '/mnt/sda1/data/TNT/analysis/07_rna_cms'; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

ex <- readr::read_tsv(RNA_X, show_col_types=FALSE)
# Use TPM for CMScaller (it expects log-transformed, but accepts normalized expression)
tpm_cols <- grep('_TPM$', colnames(ex), value=TRUE)
sym <- ex$Gene_Symbol
mat <- as.matrix(ex[, tpm_cols])
colnames(mat) <- sub('_TPM$','', colnames(mat))

# CMScaller requires Entrez IDs by default. Convert via org.Hs.eg.db
# Actually CMScaller() rowNames=TRUE argument + symbol: try "rowNames='symbol'"
# Alternatively use package's built-in entrez conversion
# Use emat with Entrez via AnnotationDbi
suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
})
# Collapse to max by symbol
rownames(mat) <- make.unique(as.character(sym))
# Get Entrez
entrez <- mapIds(org.Hs.eg.db, keys=as.character(sym), column='ENTREZID', keytype='SYMBOL', multiVals='first')
mat2 <- mat[!is.na(entrez),]
entrez2 <- entrez[!is.na(entrez)]
rownames(mat2) <- entrez2
# Collapse duplicates by max
mat2 <- mat2[!duplicated(rownames(mat2)),]
cat('Matrix for CMS:', dim(mat2), '\n')

# CMScaller expects log2 transformed usually
emat <- log2(mat2 + 1)
res <- CMScaller(emat, RNAseq=TRUE, doPlot=FALSE, templates=CMScaller::templates.CMS, rowNames='entrez')
print(head(res))

# Merge with meta
meta <- readr::read_tsv(META, show_col_types=FALSE)
res$sample_id <- rownames(res)
out <- merge(res, meta[, c('sample_id','subject_id','timepoint','response_bin','response_num')], by='sample_id')
write_tsv(out, file.path(OUT,'cms_assignments.tsv'))

# Summary per timepoint
cat('\n=== CMS frequency overall ===\n'); print(table(out$prediction))
for (tp in c('pre','post','normal')) {
  sub <- out[out$timepoint==tp,]
  if (nrow(sub)<5) next
  cat('\n=== ',tp,' CMS x response ===\n')
  print(table(sub$prediction, sub$response_bin, useNA='ifany'))
  t <- table(sub$prediction, sub$response_bin)
  if (all(dim(t)>=2)) {
    cat('Fisher p=', fisher.test(t, simulate.p.value=TRUE)$p.value, '\n')
  }
}
cat('\nSaved to', OUT, '\n')
