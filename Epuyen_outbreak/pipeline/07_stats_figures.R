# Statistics + figures analogous to paper sections 3.1-3.4, Figures 1-3.
# Run with Rscript (see README for required R packages).

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(stringr)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
base <- normalizePath(file.path(dirname(script_path), ".."))

results_dir <- file.path(base, "results")
fig_dir <- file.path(base, "figures")
dir.create(fig_dir, showWarnings = FALSE)

orf_metrics <- read_tsv(file.path(results_dir, "predecessor_orf_metrics.tsv"), show_col_types = FALSE)
hkg_metrics <- read_tsv(file.path(results_dir, "hkg_metrics.tsv"), show_col_types = FALSE)
mfe <- read_tsv(file.path(results_dir, "mfe_predecessor_hkg.tsv"), show_col_types = FALSE)
pnps <- read_tsv(file.path(results_dir, "pnps_per_patient_segment.tsv"), show_col_types = FALSE)
rand_summary <- read_tsv(file.path(results_dir, "randomization_summary.tsv"), show_col_types = FALSE)
rand_all <- read_tsv(file.path(results_dir, "randomization_all_replicates.tsv"), show_col_types = FALSE)
strain_metrics <- read_tsv(file.path(results_dir, "strain_genome_metrics.tsv"), show_col_types = FALSE)

sink(file.path(results_dir, "stats_summary.txt"))

cat("=== CAI vs mFE correlation: Lung HKG panel (n=12) ===\n")
hkg_mfe <- mfe %>% filter(group == "Lung_HKG") %>% rename(gene_symbol = name)
hkg_joined <- inner_join(hkg_metrics, hkg_mfe, by = "gene_symbol")
ct_hkg <- cor.test(hkg_joined$CAI, hkg_joined$mFE_kcal_per_mol, method = "pearson")
print(ct_hkg)

cat("\n=== CAI vs mFE correlation: Andes virus predecessor ORFs (n=3, CAVEAT: severely underpowered vs paper's n=12) ===\n")
andes_mfe <- mfe %>% filter(group == "Andes_virus_predecessor") %>%
  mutate(segment = str_extract(name, "^[A-Z]")) %>% rename(mFE_virus = mFE_kcal_per_mol)
andes_joined <- inner_join(orf_metrics, andes_mfe, by = "segment")
if (nrow(andes_joined) >= 3) {
  ct_andes <- cor.test(andes_joined$CAI, andes_joined$mFE_virus, method = "pearson")
  print(ct_andes)
} else {
  cat("Not enough ORFs for correlation.\n")
}

cat("\n=== Two-sample t-test: CAI, Andes virus predecessor ORFs (n=3) vs Lung HKG (n=12) ===\n")
tt_cai <- t.test(orf_metrics$CAI, hkg_metrics$CAI)
print(tt_cai)

cat("\n=== Two-sample t-test: mFE, Andes virus predecessor ORFs (n=3) vs Lung HKG (n=12) ===\n")
tt_mfe <- t.test(andes_mfe$mFE_virus, hkg_mfe$mFE_kcal_per_mol)
print(tt_mfe)

cat("\n=== Randomization one-tailed paired t-test (native vs best-of-100 randomized mFE, n=3 ORFs) ===\n")
rand_tt <- read_tsv(file.path(results_dir, "randomization_ttest.tsv"), show_col_types = FALSE)
print(rand_tt)

cat("\n=== pn/ps summary across all successor patients x segments ===\n")
print(summary(pnps$pn_ps_ratio))
cat(sprintf("n(pn/ps > 1, diversifying/positive selection): %d\n", sum(pnps$pn_ps_ratio > 1, na.rm = TRUE)))
cat(sprintf("n(pn/ps < 1, purifying selection): %d\n", sum(pnps$pn_ps_ratio < 1, na.rm = TRUE)))

sink()
cat("Wrote results/stats_summary.txt\n")

# ---- Figures ----

# Figure 2 analog: CAI boxplot, HKG vs Andes virus
cai_df <- bind_rows(
  hkg_metrics %>% transmute(group = "Lung HKG", CAI),
  orf_metrics %>% transmute(group = "Andes virus (predecessor ORFs)", CAI)
)
p_cai <- ggplot(cai_df, aes(x = group, y = CAI, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal(base_size = 13) +
  labs(title = "CAI: Lung HKG vs Andes virus predecessor ORFs", x = NULL, y = "CAI") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig2_analog_CAI_boxplot.png"), p_cai, width = 6, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "fig2_analog_CAI_boxplot.pdf"), p_cai, width = 6, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig2_analog_CAI_boxplot.tiff"), p_cai, width = 6, height = 5, dpi = 300, compression = "lzw")

# Figure 3 analog: mFE boxplot, HKG vs Andes native vs Andes randomized
mfe_df <- bind_rows(
  hkg_mfe %>% transmute(group = "Lung HKG", mFE = mFE_kcal_per_mol),
  andes_mfe %>% transmute(group = "Andes virus (native)", mFE = mFE_virus),
  rand_all %>% transmute(group = "Andes virus (randomized, n=300)", mFE = mFE_kcal_per_mol)
)
p_mfe <- ggplot(mfe_df, aes(x = group, y = mFE, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  theme_minimal(base_size = 13) +
  labs(title = "Folding energy: HKG vs Andes virus (native & randomized)", x = NULL, y = expression(Delta*G~(kcal/mol))) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(fig_dir, "fig3_analog_mFE_boxplot.png"), p_mfe, width = 7, height = 5, dpi = 150)
ggsave(file.path(fig_dir, "fig3_analog_mFE_boxplot.pdf"), p_mfe, width = 7, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "fig3_analog_mFE_boxplot.tiff"), p_mfe, width = 7, height = 5, dpi = 300, compression = "lzw")

# Figure 1 analog: pn/ps bar charts per segment, across successor patients
p_pnps <- ggplot(pnps %>% filter(!is.na(pn_ps_ratio)), aes(x = patient_num, y = pn_ps_ratio)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~segment, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(title = "pn/ps per successor patient, by ORF segment", x = "Patient", y = "pn/ps") +
  theme(axis.text.x = element_text(angle = 90, size = 6))
ggsave(file.path(fig_dir, "fig1_analog_pnps_bars.png"), p_pnps, width = 10, height = 6, dpi = 150)
ggsave(file.path(fig_dir, "fig1_analog_pnps_bars.pdf"), p_pnps, width = 10, height = 6, dpi = 300)
ggsave(file.path(fig_dir, "fig1_analog_pnps_bars.tiff"), p_pnps, width = 10, height = 6, dpi = 300, compression = "lzw")

cat("Wrote 3 figures (png/pdf/tiff each) to figures/\n")
