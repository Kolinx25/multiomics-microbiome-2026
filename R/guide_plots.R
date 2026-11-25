# Plots for measurements-guide.html, coloured by sequencing batch.
# Run: Rscript R/guide_plots.R  (from the project root)

library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(ragg)

source(here("R/count_tables.R"))

out_dir <- here("measurements-guide_files")
dir.create(out_dir, showWarnings = FALSE)

d <- read_counts("group_a")
metadata <- d$metadata

batch_cols <- c(B1 = "#0f766e", B2 = "#d97706")
arm_shapes <- c(probiotic = 16, placebo = 17)
save_plot <- function(p, file, w = 2400, h = 1150) {
  ragg::agg_png(file = file.path(out_dir, file), width = w, height = h,
          units = "px", res = 220, background = "white")
  print(p)
  dev.off()
}
base_theme <- theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

## 1. Library sizes by batch, with the low-read screen line -----------------
library_sizes <- bind_rows(
  tibble(table = "Microbiome", sample_id = colnames(d$species),
         library_size = colSums(d$species)),
  tibble(table = "RNA-seq",   sample_id = colnames(d$genes),
         library_size = colSums(d$genes))
) |>
  left_join(select(metadata, sample_id, time_point, arm, rnaseq_batch),
            by = "sample_id") |>
  group_by(table) |>
  mutate(lower_screen = quantile(library_size, 0.25) - 1.5 * IQR(library_size)) |>
  ungroup()

cutoffs <- library_sizes |> distinct(table, lower_screen)

p1 <- ggplot(library_sizes,
             aes(x = time_point, y = library_size, colour = rnaseq_batch)) +
  geom_hline(data = cutoffs, aes(yintercept = lower_screen), linetype = "dashed",
             colour = "#92400e", linewidth = 0.4) +
  geom_boxplot(aes(group = time_point), fill = NA, colour = "grey55",
               outlier.shape = NA, width = 0.35) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.7) +
  scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  scale_colour_manual(values = batch_cols, name = "Batch") +
  facet_wrap(~ table, scales = "free_y") +
  labs(title = "Library sizes by sequencing batch",
       subtitle = "Dashed line = low-read screen (25th percentile − 1.5 × IQR)",
       x = NULL, y = "Total counts per sample (log scale)") +
  base_theme
save_plot(p1, "plot-library-batch.png")

low_flag <- library_sizes |> group_by(table) |> summarise(flagged = sum(library_size < lower_screen))

## 2. Batch vs arm balance ---------------------------------------------------
p2 <- metadata |>
  count(rnaseq_batch, arm, time_point) |>
  ggplot(aes(x = rnaseq_batch, y = n, fill = arm)) +
  geom_col(position = position_dodge(), width = 0.7, colour = "white") +
  geom_text(aes(label = n), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3.4) +
  scale_fill_manual(values = c(probiotic = "#0f766e", placebo = "#94a3b8"),
                    name = "Arm") +
  facet_wrap(~ time_point) +
  labs(title = "Is batch mixed up with treatment group?",
       subtitle = "Want: both arms well represented in both batches, at each time point",
       x = "Sequencing batch", y = "Samples") +
  base_theme
save_plot(p2, "plot-batch-balance.png", w = 2000, h = 1000)

## 3. Microbiome ordination (Bray–Curtis PCoA) by batch ---------------------
comm <- d$species[, metadata$sample_id]
comm <- comm[rowSums(comm > 0) >= 3, ]          # drop near-absent features
pcoa <- cmdscale(vegdist(t(comm), method = "bray"), k = 2, eig = TRUE)
var_expl <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

ord <- tibble(
  PC1 = pcoa$points[, 1], PC2 = pcoa$points[, 2],
  metadata |> select(time_point, arm, rnaseq_batch)
)

p3 <- ggplot(ord, aes(PC1, PC2, colour = rnaseq_batch, shape = arm)) +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_colour_manual(values = batch_cols, name = "Batch") +
  scale_shape_manual(values = arm_shapes, name = "Arm") +
  facet_wrap(~ time_point) +
  labs(title = "Microbiome communities: does batch structure them?",
       subtitle = "Bray–Curtis PCoA; points close together have similar species profiles",
       x = paste0("PCoA 1 (", var_expl[1], "%)"),
       y = paste0("PCoA 2 (", var_expl[2], "%)")) +
  base_theme
save_plot(p3, "plot-pcoa-batch.png")

## 4. Baseline (pre) clinical measures by arm -------------------------------
measure_labels <- c(ph = "Vaginal pH", nugent = "Nugent score (0–10)",
                    clue_cells = "Clue cells (%)")
p4 <- metadata |>
  filter(time_point == "pre") |>
  select(arm, ph, nugent, clue_cells) |>
  pivot_longer(-arm, names_to = "measure", values_to = "value") |>
  ggplot(aes(arm, value, fill = arm)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.7) +
  scale_fill_manual(values = c(probiotic = "#0f766e", placebo = "#94a3b8"),
                    name = NULL) +
  facet_wrap(~ measure, scales = "free_y",
             labeller = labeller(measure = measure_labels)) +
  labs(title = "Baseline differences: the arms before the intervention",
       subtitle = "Pre-visit samples only; big gaps here mean adjust, don't credit the probiotic",
       x = NULL, y = NULL) +
  base_theme
save_plot(p4, "plot-baseline.png", w = 2000, h = 1000)

## 5. Cytokine PCA by batch -------------------------------------------------
cyto <- d$cytokines |> left_join(select(metadata, sample_id, time_point, arm, rnaseq_batch),
                                 by = "sample_id")
cyto_mat <- cyto |> select(-sample_id, -time_point, -arm, -rnaseq_batch) |>
  as.matrix()
cyto_mat <- log10(replace_na_mat <- apply(cyto_mat, 2, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE); x
}))
cat("cytokine NAs imputed with molecule medians:", sum(is.na(d$cytokines[,-1])), "\n")
pca <- prcomp(cyto_mat, center = TRUE, scale. = TRUE)
pvar <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)

p5 <- tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
             cyto |> select(time_point, arm, rnaseq_batch)) |>
  ggplot(aes(PC1, PC2, colour = rnaseq_batch, shape = arm)) +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_colour_manual(values = batch_cols, name = "Batch") +
  scale_shape_manual(values = arm_shapes, name = "Arm") +
  facet_wrap(~ time_point) +
  labs(title = "Cytokine profiles: does batch structure them?",
       subtitle = "PCA on log10 cytokine values",
       x = paste0("PC1 (", pvar[1], "%)"),
       y = paste0("PC2 (", pvar[2], "%)")) +
  base_theme
save_plot(p5, "plot-cytokine-pca.png")

## 6. Missing values per metadata column -----------------------------------
p6 <- metadata |>
  summarise(across(everything(), \(x) sum(is.na(x)))) |>
  pivot_longer(everything(), names_to = "column", values_to = "missing") |>
  ggplot(aes(reorder(column, -missing), missing, fill = missing > 0)) +
  geom_col() +
  geom_text(aes(label = missing), vjust = -0.3, size = 3.4) +
  scale_fill_manual(values = c(`FALSE` = "#0f766e", `TRUE` = "#d97706"), guide = "none") +
  labs(title = "Missing values in the sample table",
       subtitle = "Required design fields must be complete; clinical gaps need a look",
       x = NULL, y = "Number of NA values (out of 352 samples)") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p6, "plot-missing.png", w = 2000, h = 950)

## Numbers needed for the captions -----------------------------------------
cat("---- samples per batch x arm x time ----\n")
print(count(metadata, rnaseq_batch, arm, time_point))
cat("---- low-read flags ----\n"); print(low_flag)
cat("---- pre means by arm ----\n")
print(metadata |> filter(time_point == "pre") |>
        group_by(arm) |> summarise(across(c(ph, nugent, clue_cells, age), mean)))
cat("---- PCoA/PCA variance ----\n")
cat("pcoa:", var_expl, " pca:", pvar, "\n")
cat("---- adonis: ordination ~ batch + arm + time ----\n")
set.seed(1)
print(adonis2(vegdist(t(comm), "bray") ~ rnaseq_batch + arm + time_point,
              data = as.data.frame(metadata), permutations = 999))
cat("---- cytokine: same ----\n")
set.seed(1)
print(adonis2(vegdist(cyto_mat, "euclid") ~ rnaseq_batch + arm + time_point,
              data = as.data.frame(cyto), permutations = 999))
