# DESeq2 DEG good vs bad (pre-treatment) + fgsea Hallmark/Reactome
# Covariate-adjusted for sex, clinical T stage
suppressPackageStartupMessages({
  library(DESeq2); library(dplyr); library(readr)
  library(fgsea); library(msigdbr)
})

RNA_X <- '/mnt/sda1/data/TNT/TNT_RNAseq/result_RNAseq/Expression_profile/StringTie/Expression_Profile.GRCh38.gene.txt'
META  <- '/mnt/sda1/data/TNT/analysis/00_cohort/rna_inventory.tsv'
OUT   <- '/mnt/sda1/data/TNT/analysis/05_rna_deg_gsea'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# Build count matrix from xlsx via the txt export
ex <- readr::read_tsv(RNA_X, show_col_types=FALSE)
cnt_cols <- grep('_Read_Count$', colnames(ex), value=TRUE)
sym <- ex$Gene_Symbol
cts <- ex[, cnt_cols]
colnames(cts) <- sub('_Read_Count$','', colnames(cts))
cts <- as.data.frame(cts); rownames(cts) <- make.unique(as.character(sym))
cts <- cts[!is.na(rownames(cts)) & rownames(cts)!='', ]
# drop non-integer / collapse by max
cts <- round(as.matrix(cts))
storage.mode(cts) <- 'integer'
cat('count matrix:', dim(cts), '\n')

meta <- readr::read_tsv(META, show_col_types=FALSE)
meta <- meta[match(colnames(cts), meta$sample_id),]
meta$response_bin <- factor(meta$response_bin, levels=c('bad','good'))
meta$sex <- factor(meta$sex, levels=c('M','F'))
meta$cT_simple <- factor(ifelse(meta$cT %in% c('T2','T2/T3'), 'T2-3low', meta$cT), levels=c('T2-3low','T3','T4'))
meta$timepoint <- factor(meta$timepoint, levels=c('normal','pre','post'))

# Subset to PRE samples
keep <- which(meta$timepoint=='pre')
md <- meta[keep,]
co <- cts[, md$sample_id]
cat('pre samples:', ncol(co), ' good:', sum(md$response_bin=='good'), ' bad:', sum(md$response_bin=='bad'), '\n')

# Filter low-count genes
co <- co[rowSums(co>=10) >= 3, ]
cat('genes after filter:', nrow(co), '\n')

dds <- DESeqDataSetFromMatrix(countData=co, colData=md,
                              design= ~ sex + cT_simple + response_bin)
dds <- DESeq(dds, parallel=TRUE)
res <- results(dds, name='response_bin_good_vs_bad')
res_df <- as.data.frame(res) %>% tibble::rownames_to_column('gene') %>%
  arrange(pvalue)
write_tsv(res_df, file.path(OUT,'DEG_good_vs_bad_pre.tsv'))
cat('top 10 DEG:\n'); print(head(res_df, 10))

# fgsea with Hallmark + Reactome
ranks <- setNames(res_df$log2FoldChange, res_df$gene)
ranks <- ranks[!is.na(ranks) & !is.na(res_df$pvalue)]
ranks <- sort(ranks, decreasing=TRUE)

for (cat1 in list(c('H','Hallmark'), c('C2','Reactome'))) {
  if (cat1[1]=='H') {
    gs <- msigdbr(species='Homo sapiens', collection='H')
  } else {
    gs <- msigdbr(species='Homo sapiens', collection='C2', subcollection='CP:REACTOME')
  }
  gs_list <- split(gs$gene_symbol, gs$gs_name)
  fg <- fgsea(pathways=gs_list, stats=ranks, minSize=10, maxSize=500)
  fg <- fg[order(pval)]
  write_tsv(as.data.frame(fg)[, !colnames(fg) %in% 'leadingEdge'], file.path(OUT, paste0('GSEA_',cat1[2],'_pre.tsv')))
  cat('\n=== GSEA ',cat1[2],' top ===\n'); print(head(fg[, .(pathway, pval, padj, NES, size)], 10))
}
cat('\nDONE\n')
