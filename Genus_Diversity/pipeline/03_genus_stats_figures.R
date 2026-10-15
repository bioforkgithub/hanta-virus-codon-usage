# Across-species statistics + figures: pairwise pn/ps by comparison-type
# (pathogenic-pathogenic vs environmental-environmental vs cross), and
# per-species CAI/ENC by clinical category.
#
# The quantity is pn/ps, uncorrected, not dN/dS - see the header of
# 02_pairwise_dnds.py for the relabelling applied at revision.
#
# CAVEAT (stated explicitly in the manuscript too): pairwise pn/ps values
# are not independent observations - each of the 22 species appears in 21
# pairs, so a naive ANOVA/t-test here violates the independence assumption
# a full phylogenetic comparative method would respect. Treat these as
# descriptive/exploratory statistics, not a phylogenetically corrected test.
# The p-values reported in the manuscript come from 05_permutation_test.py.
#
# All t.test() calls below are Welch tests on unequal variances, which is
# R's default for two samples and is what the manuscript's Methods state.

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

dnds <- read_tsv(file.path(results_dir, "pairwise_pnps_matrix.tsv"), show_col_types = FALSE) %>%
  filter(pnps != "NA") %>% mutate(pnps = as.numeric(pnps), ps = as.numeric(ps), pn = as.numeric(pn))
cai_enc <- read_tsv(file.path(results_dir, "genus_diversity_cai_enc.tsv"), show_col_types = FALSE)

sink(file.path(results_dir, "genus_stats_summary.txt"))

cat("=== Pairwise pn/ps by comparison type (n pairs, mean, sd) ===\n")
print(dnds %>% group_by(comparison_type) %>% summarise(n = n(), mean_pnps = mean(pnps), sd_pnps = sd(pnps), mean_ps = mean(ps), mean_pn = mean(pn)))

cat("\n=== Kruskal-Wallis test: pn/ps ~ comparison_type (descriptive only - see caveat in script header) ===\n")
print(kruskal.test(pnps ~ comparison_type, data = dnds))

cat("\n=== Pairwise Wilcoxon rank-sum tests between comparison types ===\n")
print(pairwise.wilcox.test(dnds$pnps, dnds$comparison_type, p.adjust.method = "BH"))

cat("\n=== Per-species CAI (vs human high-expression reference) and ENC, by clinical category ===\n")
print(cai_enc %>% group_by(category) %>% summarise(n = n(), mean_CAI = mean(CAI_vs_human_ref), mean_ENC = mean(ENC)))

cat("\n=== Welch two-sample t-test: CAI by category (manuscript Table 7) ===\n")
print(t.test(CAI_vs_human_ref ~ category, data = cai_enc))

cat("\n=== Welch two-sample t-test: ENC by category ===\n")
print(t.test(ENC ~ category, data = cai_enc))

cat("\n=== CAI by category, adjusted for GC3 (analysis of covariance) ===\n")
# Added at revision. CAI and GC3 are strongly collinear, so the unadjusted
# comparison above cannot separate a codon-preference effect from a base-
# composition one. The adjusted test is the clinical-category term of
# CAI ~ category + GC3. Because the term carries one degree of freedom and is
# entered last, Type I, II and III sums of squares all give the same p-value.
#
# The submitted version adjusted in two steps instead - fit CAI ~ GC3 across
# the genes of the comparison, then t-test the two groups' residuals - which
# does not propagate the first-stage fit's uncertainty into the second-stage
# test. That procedure is reproduced here for comparability only; the value
# reported in the revised manuscript is the covariance-analysis one. The
# conclusion is the same under both.
if (!"GC3" %in% names(cai_enc)) {
  cat("GC3 column absent - rerun 02_pairwise_dnds.py to emit it.\n")
} else {
  fit <- lm(CAI_vs_human_ref ~ category + GC3, data = cai_enc)
  print(summary(fit))
  cat("\nType II sums of squares on the same model:\n")
  print(anova(fit))
  cat("\nCAI ~ GC3 alone, for the collinearity it implies:\n")
  print(summary(lm(CAI_vs_human_ref ~ GC3, data = cai_enc)))
  cat("\nSubmitted two-step residual procedure, for comparison only:\n")
  resid_two_step <- residuals(lm(CAI_vs_human_ref ~ GC3, data = cai_enc))
  print(t.test(resid_two_step ~ cai_enc$category))
}

sink()
cat("Wrote results/genus_stats_summary.txt\n")

# ---- Figures ----

p_dnds <- ggplot(dnds, aes(x = comparison_type, y = pnps, fill = comparison_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  theme_minimal(base_size = 13) +
  labs(title = "Pairwise nucleoprotein pn/ps across 22 hantavirus species",
       subtitle = "By clinical category of the species pair (22 species, 231 pairs)",
       x = NULL, y = "pn/ps (Nei-Gojobori, uncorrected)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.png"), p_dnds, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.pdf"), p_dnds, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig4_genus_dnds_boxplot.tiff"), p_dnds, width = 7, height = 5.5, dpi = 300, compression = "lzw")

p_cai <- ggplot(cai_enc, aes(x = category, y = CAI_vs_human_ref, fill = category)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  theme_minimal(base_size = 13) +
  labs(title = "Nucleoprotein CAI by hantavirus clinical category",
       subtitle = "vs. human high-expression reference, 22 species",
       x = NULL, y = "CAI") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.png"), p_cai, width = 7, height = 5.5, dpi = 150)
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.pdf"), p_cai, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_dir, "fig5_genus_CAI_by_category.tiff"), p_cai, width = 7, height = 5.5, dpi = 300, compression = "lzw")

cat("Wrote 2 figures (png/pdf/tiff each) to figures/\n")
