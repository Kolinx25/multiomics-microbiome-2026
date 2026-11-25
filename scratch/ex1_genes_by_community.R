# Exercise 1: gene table by community (week 18), DESeq2
# Follows notebook 04, step 1. Raw counts in, DESeq2 does its own normalisation.
library(tidyverse)
library(here)
library(DESeq2)

source(here("R/count_tables.R"))

d <- read_counts("instructional")
sample_groups <- read_csv(here("data/my_sample_groups.csv"))
metadata <- d$metadata |> left_join(sample_groups, by = "sample_id")

dplyr::count(metadata, time_point, community)

# ---- sample order check --------------------------------------------------
# Keep only week 18, then force the count-matrix columns to match the
# metadata rows in the same order. DESeq2 does NOT check this for you.
week18 <- metadata$sample_id[metadata$time_point == "week_18"]
counts_w18 <- d$genes[, week18]
meta_w18 <- metadata[match(week18, metadata$sample_id), ]

cat("Order check (should be TRUE):",
    identical(colnames(counts_w18), meta_w18$sample_id), "\n")
cat("First 3 sample ids:", paste(head(week18, 3), collapse = ", "), "\n")
cat("Samples per community:\n")
print(table(meta_w18$community))

# ---- gene filter ---------------------------------------------------------
keep <- rowSums(counts_w18 >= 10) >= ncol(counts_w18) / 4
cat("Genes kept by filter:", sum(keep), "of", nrow(counts_w18), "\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_w18[keep, ],
  colData  = meta_w18,
  design   = ~ community
)
dds <- DESeq(dds)

cat("DESeq2 used this level order as the reference:",
    paste(levels(dds$community), collapse = " < "), "\n")

# Contrast: anaerobe against lactobacillus.
# POSITIVE log2FoldChange = higher in anaerobe-dominant samples.
# NEGATIVE = higher in lactobacillus-dominant samples.
res <- results(dds, contrast = c("community", "anaerobe", "lactobacillus"))
res_df <- as.data.frame(res) |>
  rownames_to_column("feature") |>
  left_join(d$gene_info, by = c("feature" = "ensembl_gene_id")) |>
  select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
         symbol, gene_name)

# ---- direction check: does the sign match the raw medians? ---------------
# Independent check on a handful of top genes: median CPM in each group.
top10 <- res_df |> arrange(padj) |> head(10) |> pull(feature)
totals <- colSums(counts_w18)
cpm <- t(t(counts_w18) / totals * 1e6)
dir_check <- sapply(top10, function(g) {
  a <- median(cpm[g, meta_w18$community == "anaerobe"])
  l <- median(cpm[g, meta_w18$community == "lactobacillus"])
  log2(a / l)   # rough log2FC; sign should match the DESeq2 log2FoldChange
})
check_tbl <- res_df |> filter(feature %in% top10) |>
  arrange(padj) |>
  select(feature, symbol, log2FoldChange) |>
  mutate(median_cpm_log2fc = round(dir_check[feature], 2),
         signs_agree = sign(log2FoldChange) == sign(dir_check[feature]))
# Print as plain text so NaN/Inf values cannot break the console output
cat("Direction check, top 10 genes (does the DESeq2 sign match the raw medians?)\n")
print(as.data.frame(lapply(check_tbl, function(x) format(x, trim = TRUE))),
      row.names = FALSE)

cat("Genes tested (non-NA p):", sum(!is.na(res_df$pvalue)), "\n")
saveRDS(res_df, here("scratch/res_genes_by_community.rds"))
saveRDS(list(counts = counts_w18, meta = meta_w18),
        here("scratch/w18_counts_meta.rds"))
write_csv(res_df, here("data/genes_by_community.csv"))
cat("Saved data/genes_by_community.csv\n")
