# Three-way tissue comparison: Andes virus index-case ORFs vs Lung HKG vs
# Kidney HKG, for CAI and mFE - hantaviruses cause both pulmonary (HPS) and
# renal (HFRS/AKI) disease, so both tissue panels are relevant benchmarks.

suppressMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
base <- normalizePath(file.path(dirname(script_path), ".."))

results_dir <- file.path(base, "results")
fig_dir <- file.path(base, "figures")

orf_metrics <- read_tsv(file.path(results_dir, "predecessor_orf_metrics.tsv"), show_col_types = FALSE)
lung_hkg <- read_tsv(file.path(results_dir, "hkg_metrics.tsv"), show_col_types = FALSE)
kidney_hkg <- read_tsv(file.path(results_dir, "kidney_hkg_metrics.tsv"), show_col_types = FALSE)
mfe_lung <- read_tsv(file.path(results_dir, "mfe_predecessor_hkg.tsv"), show_col_types = FALSE)
mfe_kidney <- read_tsv(file.path(results_dir, "mfe_kidney_hkg.tsv"), show_col_types = FALSE)

cai_df <- bind_rows(
  orf_metrics %>% transmute(group = "Andes virus (index-case ORFs)", CAI, ENC),
  lung_hkg %>% transmute(group = "Lung HKG", CAI, ENC),
  kidney_hkg %>% transmute(group = "Kidney HKG", CAI, ENC)
)
mfe_df <- bind_rows(
  mfe_lung %>% filter(group == "Andes_virus_predecessor") %>% transmute(group = "Andes virus (index-case ORFs)", mFE = mFE_kcal_per_mol),
  mfe_lung %>% filter(group == "Lung_HKG") %>% transmute(group = "Lung HKG", mFE = mFE_kcal_per_mol),
  mfe_kidney %>% transmute(group = "Kidney HKG", mFE = mFE_kcal_per_mol)
)

sink(file.path(results_dir, "tissue_comparison_stats.txt"))

cat("=== CAI by group (Andes virus ORFs / Lung HKG / Kidney HKG) ===\n")
print(cai_df %>% group_by(group) %>% summarise(n = n(), mean_CAI = mean(CAI), mean_ENC = mean(ENC)))

cat("\n=== One-way ANOVA: CAI ~ group ===\n")
print(summary(aov(CAI ~ group, data = cai_df)))

cat("\n=== Andes virus ORFs vs Kidney HKG: CAI t-test ===\n")
print(t.test(orf_metrics$CAI, kidney_hkg$CAI))

cat("\n=== Andes virus ORFs vs Kidney HKG: ENC t-test ===\n")
print(t.test(orf_metrics$ENC, kidney_hkg$ENC))

cat("\n=== mFE by group ===\n")
print(mfe_df %>% group_by(group) %>% summarise(n = n(), mean_mFE = mean(mFE)))

cat("\n=== One-way ANOVA: mFE ~ group ===\n")
print(summary(aov(mFE ~ group, data = mfe_df)))

cat("\n=== Andes virus ORFs vs Kidney HKG: mFE t-test ===\n")
print(t.test(mfe_df$mFE[mfe_df$group == "Andes virus (index-case ORFs)"], mfe_df$mFE[mfe_df$group == "Kidney HKG"]))

sink()
cat("Wrote results/tissue_comparison_stats.txt\n")

# wrap the long "Andes virus (index-case ORFs)" label onto two lines so it
# isn't clipped by the plot panel at these figure widths
relabel <- function(g) factor(
  recode(g, "Andes virus (index-case ORFs)" = "Andes virus\n(index-case ORFs)"),
  levels = c("Andes virus\n(index-case ORFs)", "Kidney HKG", "Lung HKG")
)
cai_df$group <- relabel(cai_df$group)
mfe_df$group <- relabel(mfe_df$group)

p_cai <- ggplot(cai_df, aes(x = group, y = CAI, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal(base_size = 13) +
  labs(title = "CAI: Andes virus vs Lung HKG vs Kidney HKG", x = NULL, y = "CAI") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig6_tissue_CAI_boxplot.png"), p_cai, width = 7.5, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig6_tissue_CAI_boxplot.pdf"), p_cai, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig6_tissue_CAI_boxplot.tiff"), p_cai, width = 7.5, height = 5.5, dpi = 300, compression = "lzw")

p_mfe <- ggplot(mfe_df, aes(x = group, y = mFE, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal(base_size = 13) +
  labs(title = "Folding energy: Andes virus vs Lung HKG vs Kidney HKG", x = NULL, y = expression(Delta*G~(kcal/mol))) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig7_tissue_mFE_boxplot.png"), p_mfe, width = 7.5, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig7_tissue_mFE_boxplot.pdf"), p_mfe, width = 7.5, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig7_tissue_mFE_boxplot.tiff"), p_mfe, width = 7.5, height = 5.5, dpi = 300, compression = "lzw")

cat("Wrote 2 figures (png/pdf/tiff each) to figures/\n")
