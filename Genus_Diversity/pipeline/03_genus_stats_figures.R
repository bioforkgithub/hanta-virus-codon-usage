# Genus-wide statistics + figures: pairwise dN/dS by comparison-type
# (pathogenic-pathogenic vs environmental-environmental vs cross), and
# per-species CAI/ENC by clinical category.
#
# CAVEAT (stated explicitly in the manuscript too): pairwise dN/dS values
# are not independent observations - each of the 22 species appears in 21
# pairs, so a naive ANOVA/t-test here violates the independence assumption
# a full phylogenetic comparative method would respect. Treat these as
# descriptive/exploratory statistics, not a phylogenetically corrected test.

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
dir.create(fig_dir, showWarnings = FALSE)

dnds <- read_tsv(file.path(results_dir, "pairwise_dnds_matrix.tsv"), show_col_types = FALSE) %>%
  filter(dNdS != "NA") %>% mutate(dNdS = as.numeric(dNdS), ps = as.numeric(ps), pn = as.numeric(pn))
cai_enc <- read_tsv(file.path(results_dir, "genus_diversity_cai_enc.tsv"), show_col_types = FALSE)

sink(file.path(results_dir, "genus_stats_summary.txt"))

cat("=== Pairwise dN/dS by comparison type (n pairs, mean, sd) ===\n")
print(dnds %>% group_by(comparison_type) %>% summarise(n = n(), mean_dNdS = mean(dNdS), sd_dNdS = sd(dNdS), mean_ps = mean(ps), mean_pn = mean(pn)))

cat("\n=== Kruskal-Wallis test: dN/dS ~ comparison_type (descriptive only - see caveat in script header) ===\n")
print(kruskal.test(dNdS ~ comparison_type, data = dnds))

cat("\n=== Pairwise Wilcoxon rank-sum tests between comparison types ===\n")
print(pairwise.wilcox.test(dnds$dNdS, dnds$comparison_type, p.adjust.method = "BH"))

cat("\n=== Per-species CAI (vs human high-expression reference) and ENC, by clinical category ===\n")
print(cai_enc %>% group_by(category) %>% summarise(n = n(), mean_CAI = mean(CAI_vs_human_ref), mean_ENC = mean(ENC)))

cat("\n=== Two-sample t-test: CAI by category ===\n")
print(t.test(CAI_vs_human_ref ~ category, data = cai_enc))

cat("\n=== Two-sample t-test: ENC by category ===\n")
print(t.test(ENC ~ category, data = cai_enc))

sink()
cat("Wrote results/genus_stats_summary.txt\n")

# ---- Figures ----

p_dnds <- ggplot(dnds, aes(x = comparison_type, y = dNdS, fill = comparison_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  theme_minimal(base_size = 13) +
  labs(title = "Pairwise nucleoprotein dN/dS across the Hantavirus genus",
       subtitle = "By clinical category of the species pair (22 species, 231 pairs)",
       x = NULL, y = "dN/dS (Nei-Gojobori)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.png"), p_dnds, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.pdf"), p_dnds, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.tiff"), p_dnds, width = 7, height = 5.5, dpi = 300, compression = "lzw")

p_cai <- ggplot(cai_enc, aes(x = category, y = CAI_vs_human_ref, fill = category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal(base_size = 13) +
  labs(title = "Nucleoprotein CAI by Hantavirus clinical category",
       subtitle = "vs. human high-expression reference, genus-wide (22 species)",
       x = NULL, y = "CAI") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.png"), p_cai, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.pdf"), p_cai, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.tiff"), p_cai, width = 7, height = 5.5, dpi = 300, compression = "lzw")

cat("Wrote 2 figures (png/pdf/tiff each) to figures/\n")
