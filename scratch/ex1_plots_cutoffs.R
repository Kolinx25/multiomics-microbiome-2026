# Exercise 1, steps 2-3: summary numbers, four plots, cutoff grid, filters
# Reads the saved result so DESeq2 does not need to run again.
library(tidyverse)
library(here)
library(DESeq2)

sp <- function(x) print(as.data.frame(lapply(x, function(c) format(c, trim = TRUE))), row.names = FALSE)

res_df <- readRDS(here("scratch/res_genes_by_community.rds"))
w <- readRDS(here("scratch/w18_counts_meta.rds"))
counts_w18 <- w$counts; meta_w18 <- w$meta

# ---- short report --------------------------------------------------------
res_df |>
  summarise(
    pvalue_lt_0.05 = sum(pvalue < 0.05, na.rm = TRUE),
    padj_lt_0.05   = sum(padj < 0.05, na.rm = TRUE),
    more_than_2fold = sum(padj < 0.05 & abs(log2FoldChange) > 1, na.rm = TRUE),
    higher_in_anaerobe      = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
    higher_in_lactobacillus = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE)
  ) |> sp()

cat("\nTop 20 by padj:\n")
res_df |> arrange(padj) |> head(20) |>
  select(symbol, gene_name, log2FoldChange, padj) |> sp()

# ---- four plots ----------------------------------------------------------
com_col <- c("anaerobe" = "#B2182B", "lactobacillus" = "#2166AC")
sig_line <- sum(!is.na(res_df$pvalue)) / 10
png(here("scratch/ex1_four_plots.png"), width = 1600, height = 1200)
par(mar = c(1, 1, 1, 1))

# 1. Volcano. Positive log2FC = higher in anaerobe.
p1 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(colour = padj < 0.05), alpha = 0.4, size = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_colour_manual(values = c("grey70", "black"), name = "padj < 0.05") +
  labs(title = "Volcano", x = "log2 fold change (positive = higher in anaerobe)",
       y = "-log10(p-value)") +
  theme_bw()

# 2. MA plot
p2 <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(colour = padj < 0.05), alpha = 0.4, size = 0.6) +
  scale_x_log10() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "MA", x = "mean expression (log scale)", y = "log2 fold change") +
  theme_bw()

# 3. p-value histogram with the flat no-difference line
p3 <- ggplot(res_df, aes(x = pvalue)) +
  geom_histogram(binwidth = 0.05, boundary = 0, colour = "white") +
  geom_hline(yintercept = sig_line, linetype = "dashed") +
  labs(title = "Raw p-values", subtitle = "dashed = genes tested / 10",
       x = "raw p-value", y = "number of genes") +
  theme_bw()

# 4. Boxplots of top 6 genes by padj (CPM then log2)
top6 <- res_df |> arrange(padj) |> head(6) |> pull(feature)
totals <- colSums(counts_w18)
box_df <- t(t(counts_w18[top6, ]) / totals * 1e6) |>
  as_tibble(rownames = "feature") |>
  pivot_longer(-feature, names_to = "sample_id", values_to = "cpm") |>
  mutate(log_cpm = log2(cpm + 1)) |>
  left_join(select(meta_w18, sample_id, community), by = "sample_id") |>
  left_join(select(res_df, feature, symbol), by = "feature")
p4 <- ggplot(box_df, aes(x = community, y = log_cpm, fill = community)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ symbol) +
  scale_fill_manual(values = com_col) +
  labs(title = "Top 6 genes", y = "log2(cpm + 1)") +
  theme_bw()

print(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2))
dev.off()
cat("Saved scratch/ex1_four_plots.png\n")

# ---- step 4: cutoffs crossed --------------------------------------------
# padj 0.01 / 0.05 / 0.1  x  |log2FC| 0 / 1.5 / 2 / 4 (0 = fold change ignored)
grid <- expand.grid(padj_cut = c(0.01, 0.05, 0.10), lfc_cut = c(0, 1.5, 2, 4))
cutoff_tbl <- grid |>
  rowwise() |>
  mutate(
    n = sum(res_df$padj < padj_cut & abs(res_df$log2FoldChange) > lfc_cut,
            na.rm = TRUE),
    up   = sum(res_df$padj < padj_cut & res_df$log2FoldChange >  lfc_cut,
               na.rm = TRUE),
    down = sum(res_df$padj < padj_cut & res_df$log2FoldChange < -lfc_cut,
               na.rm = TRUE)
  ) |>
  ungroup()
sp(cutoff_tbl)

# ---- gene filters: looser and stricter -----------------------------------
# 10% of samples: gene with count >= 10 in 10% of samples
# 50% of samples: the stricter version (used 25% in the main run)
keep10 <- rowSums(counts_w18 >= 10) >= ncol(counts_w18) / 10
keep50 <- rowSums(counts_w18 >= 10) >= ncol(counts_w18) / 2
cat("Genes kept: 10% filter", sum(keep10),
    "| 25% filter (main run)", sum(keep10 & keep50),
    "| 50% filter", sum(keep50), "\n")

filter_run <- function(keep_vec) {
  dds_f <- DESeqDataSetFromMatrix(countData = counts_w18[keep_vec, ],
                                  colData = meta_w18, design = ~ community)
  dds_f <- DESeq(dds_f, quiet = TRUE)
  as.data.frame(results(dds_f, contrast = c("community", "anaerobe",
                                            "lactobacillus")))
}
res10 <- filter_run(keep10)
res25 <- res_df  # main run, already computed
res50 <- filter_run(keep50)

filter_tbl <- tibble(
  filter = c("10% of samples", "25% (main run)", "50% of samples"),
  genes_tested = c(sum(!is.na(res10$pvalue)), sum(!is.na(res25$pvalue)),
                   sum(!is.na(res50$pvalue))),
  padj_lt_0.05 = c(sum(res10$padj < 0.05, na.rm = TRUE),
                   sum(res25$padj < 0.05, na.rm = TRUE),
                   sum(res50$padj < 0.05, na.rm = TRUE)),
  `padj_lt_0.05_&_2fold` = c(sum(res10$padj < 0.05 & abs(res10$log2FoldChange) > 1,
                                 na.rm = TRUE),
                             sum(res25$padj < 0.05 & abs(res25$log2FoldChange) > 1,
                                 na.rm = TRUE),
                             sum(res50$padj < 0.05 & abs(res50$log2FoldChange) > 1,
                                 na.rm = TRUE))
)
sp(filter_tbl)
saveRDS(filter_tbl, here("scratch/filter_tbl.rds"))
