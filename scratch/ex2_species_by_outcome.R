# Exercise 2, steps 1-3: species table by outcome (week 18), and size factors
# Default run first, then poscounts, then compare. Saves species_by_outcome.csv
library(tidyverse)
library(here)
library(DESeq2)

source(here("R/count_tables.R"))
sp <- function(x) print(as.data.frame(lapply(x, function(c) format(c, trim = TRUE))),
                        row.names = FALSE)

d <- read_counts("instructional")
metadata <- d$metadata   # outcome comes from the sample table itself

week18 <- metadata$sample_id[metadata$time_point == "week_18"]
counts_w18 <- d$species[, week18]
meta_w18 <- metadata[match(week18, metadata$sample_id), ]
stopifnot(identical(colnames(counts_w18), meta_w18$sample_id))
cat("Order check:", identical(colnames(counts_w18), meta_w18$sample_id), "\n")
cat("Outcome groups:\n"); print(table(meta_w18$outcome))

# Species present (count > 0) in at least 10% of the week 18 samples
keep10 <- rowSums(counts_w18 > 0) >= ceiling(ncol(counts_w18) * 0.1)
keep_all <- rowSums(counts_w18 > 0) == ncol(counts_w18)   # no zeros anywhere
cat("Species surviving 10% filter:", sum(keep10),
    "| present in all 157:", sum(keep_all), "\n")

run_dds <- function(sf_type = c("default", "poscounts")) {
  dds <- DESeqDataSetFromMatrix(countData = counts_w18[keep10, ],
                                colData = meta_w18, design = ~ outcome)
  if (sf_type == "poscounts") {
    sizeFactors(dds) <- estimateSizeFactorsForMatrix(counts(dds),
                                                     type = "poscounts")
  } else {
    dds <- estimateSizeFactors(dds)   # default: geometric-mean, needs a control species
  }
  DESeq(dds, quiet = TRUE)
}

fmt_res <- function(dds) {
  as.data.frame(results(dds, contrast = c("outcome", "preterm", "term"))) |>
    rownames_to_column("feature") |>
    left_join(d$taxonomy, by = c("feature" = "taxon")) |>
    select(feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj,
           Genus, Family)
}

dds_def <- run_dds("default")
res_def <- fmt_res(dds_def)
cat("\n== DEFAULT run ==\n")
cat("Species tested:", sum(!is.na(res_def$pvalue)),
    "| padj < 0.05:", sum(res_def$padj < 0.05, na.rm = TRUE), "\n")
sp(res_def |> filter(padj < 0.05) |> arrange(padj))

# ---- step 2: what are the size factors doing? ----------------------------
sf_def <- sizeFactors(dds_def)
lib_size <- colSums(counts_w18[keep10, ])
iners <- counts_w18["lactobacillus_iners", ]
iners_share <- iners / lib_size

# Which species have no zero in any sample? Those are the only ones the
# default (geometric-mean) size factors can actually "see".
cat("\nSpecies with no zero in any week-18 sample:", sum(keep_all), "\n")
cat("Of the 53 tested, how many have no zeros:",
    sum(keep_all[keep10]), "\n")

gm <- exp(rowMeans(log(counts_w18[keep10, ])))   # geometric mean per species
cat("Species among the 53 whose geometric mean is 0",
    "(invisible to default size factors):",
    sum(gm == 0), "\n")

sf_df <- tibble(sample = names(sf_def), size_factor = sf_def,
                log_lib_size = log(lib_size[sample]),
                log_iners_share = log(iners_share[sample] + 1e-9))
png(here("scratch/ex2_size_factors.png"), width = 1100, height = 500)
p1 <- ggplot(sf_df, aes(x = log_iners_share, y = log(size_factor))) +
  geom_point(alpha = 0.6, colour = "#B2182B") +
  geom_smooth(method = "lm") +
  labs(title = "Size factor vs L. iners share", x = "log(L. iners share)",
       y = "log(size factor)") + theme_bw()
p2 <- ggplot(sf_df, aes(x = log_lib_size, y = log(size_factor))) +
  geom_point(alpha = 0.6, colour = "#2166AC") +
  geom_smooth(method = "lm") +
  labs(title = "Size factor vs library size", x = "log(library size)",
       y = "log(size factor)") + theme_bw()
print(gridExtra::grid.arrange(p1, p2, nrow = 1))
dev.off()
cat("Cor(size factor, L. iners share):",
    round(cor(sf_df$size_factor, iners_share[sf_df$sample]), 3), "\n")
cat("Cor(size factor, library size):",
    round(cor(sf_df$size_factor, lib_size[sf_df$sample]), 3), "\n")

# ---- step 3: poscounts rerun, side by side -------------------------------
dds_pos <- run_dds("poscounts")
res_pos <- fmt_res(dds_pos)

cat("\n== POSCOUNTS run ==\n")
cat("Species tested:", sum(!is.na(res_pos$pvalue)),
    "| padj < 0.05:", sum(res_pos$padj < 0.05, na.rm = TRUE), "\n")

comp <- res_def |>
  select(feature, Genus, Family, lfc_def = log2FoldChange, padj_def = padj) |>
  left_join(res_pos |>
              select(feature, lfc_pos = log2FoldChange, padj_pos = padj),
            by = "feature") |>
  mutate(sig_def = padj_def < 0.05, sig_pos = padj_pos < 0.05,
         agree = sig_def == sig_pos)
cat("Significant in default only:", sum(comp$sig_def & !comp$sig_pos, na.rm = TRUE),
    "| poscounts only:", sum(comp$sig_pos & !comp$sig_def, na.rm = TRUE),
    "| both:", sum(comp$sig_def & comp$sig_pos, na.rm = TRUE), "\n")
sp(comp |> filter(sig_def | sig_pos) |> arrange(desc(padj_def)))

# Save the poscounts version for tomorrow
write_csv(res_pos, here("data/species_by_outcome.csv"))
cat("\nSaved data/species_by_outcome.csv\n")
