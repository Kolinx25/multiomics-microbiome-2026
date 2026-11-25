# Exercise 2, step 5: gene table by outcome (week 18). Saves genes_by_outcome.csv
library(tidyverse)
library(here)
library(DESeq2)

source(here("R/count_tables.R"))
d <- read_counts("instructional")
metadata <- d$metadata

week18 <- metadata$sample_id[metadata$time_point == "week_18"]
counts_w18 <- d$genes[, week18]
meta_w18 <- metadata[match(week18, metadata$sample_id), ]
stopifnot(identical(colnames(counts_w18), meta_w18$sample_id))

keep <- rowSums(counts_w18 >= 10) >= ncol(counts_w18) / 4

dds <- DESeqDataSetFromMatrix(countData = counts_w18[keep, ],
                              colData = meta_w18, design = ~ outcome)
dds <- DESeq(dds)
# Positive log2FoldChange = higher in preterm.
res <- results(dds, contrast = c("outcome", "preterm", "term"))
res_g <- as.data.frame(res) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)

write_csv(res_g, here("data/genes_by_outcome.csv"))
cat("Genes tested:", sum(!is.na(res_g$pvalue)),
    "| padj < 0.05:", sum(res_g$padj < 0.05, na.rm = TRUE),
    "| padj < 0.05 & >2-fold:",
    sum(res_g$padj < 0.05 & abs(res_g$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("Saved data/genes_by_outcome.csv\n")
