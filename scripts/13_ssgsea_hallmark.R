# ssGSEA per-sample pathway scores via GSVA (Hallmark + key curated sets)
suppressPackageStartupMessages({
  library(GSVA); library(msigdbr); library(readr); library(dplyr)
})
RNA_X <- '/mnt/sda1/data/TNT/TNT_RNAseq/result_RNAseq/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.txt'
META  <- '/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv'
OUT   <- '/mnt/sda1/data/TNT/analysis/08_rna_pathway'; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

ex <- readr::read_tsv(RNA_X, show_col_types=FALSE)
tpm_cols <- grep('_TPM$', colnames(ex), value=TRUE)
sym <- ex$Gene_Symbol
mat <- as.matrix(ex[, tpm_cols])
colnames(mat) <- sub('_TPM$','', colnames(mat))
rownames(mat) <- make.unique(as.character(sym))
mat <- mat[!is.na(rownames(mat)) & rownames(mat)!='', ]
# log2 TPM
lmat <- log2(mat + 1)
cat('matrix:', dim(lmat), '\n')

# Gene sets
hm <- msigdbr(species='Homo sapiens', collection='H')
gsH <- split(hm$gene_symbol, hm$gs_name)

# Also import key immune/stromal/proliferation curated sets
keep <- c('REACTOME_CELL_CYCLE_CHECKPOINTS','REACTOME_DNA_REPAIR','REACTOME_HOMOLOGY_DIRECTED_REPAIR',
          'REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR','REACTOME_S_PHASE','REACTOME_M_PHASE',
          'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_INTERFERON_ALPHA_BETA_SIGNALING',
          'REACTOME_INTERFERON_GAMMA_SIGNALING','REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION',
          'REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION')
re <- msigdbr(species='Homo sapiens', collection='C2', subcollection='CP:REACTOME')
gsR <- split(re$gene_symbol, re$gs_name)
gsR <- gsR[names(gsR) %in% keep]

gs <- c(gsH, gsR)
cat('pathways:', length(gs), '\n')

par <- gsvaParam(lmat, gs, kcdf='Gaussian')
scores <- gsva(par)
cat('score matrix:', dim(scores), '\n')
df <- as.data.frame(t(scores))
df$sample_id <- rownames(df)
df <- df[, c(ncol(df), 1:(ncol(df)-1))]
write_tsv(df, file.path(OUT,'ssgsea_scores.tsv'))

# Response association
meta <- readr::read_tsv(META, show_col_types=FALSE)
m <- merge(df, meta[, c('sample_id','timepoint','response_bin')], by='sample_id')

rows <- list()
for (tp in c('pre','post')) {
  sub <- m[m$timepoint==tp,]
  if (nrow(sub)<8) next
  for (p in names(gs)) {
    g <- sub[[p]][sub$response_bin=='good']
    b <- sub[[p]][sub$response_bin=='bad']
    if (length(g)<2 || length(b)<2) next
    u <- wilcox.test(g, b, exact=FALSE)
    rows[[length(rows)+1]] <- data.frame(timepoint=tp, pathway=p,
      n_good=length(g), n_bad=length(b),
      mean_good=mean(g), mean_bad=mean(b),
      delta=mean(g)-mean(b), pvalue=u$p.value)
  }
}
res <- do.call(rbind, rows)
res$qvalue <- p.adjust(res$pvalue, 'fdr')
res <- res[order(res$pvalue),]
write_tsv(res, file.path(OUT,'ssgsea_response_stats.tsv'))
cat('\n=== ssGSEA pre top ===\n')
print(head(res[res$timepoint=='pre',], 15))
cat('\n=== ssGSEA post top ===\n')
print(head(res[res$timepoint=='post',], 10))
