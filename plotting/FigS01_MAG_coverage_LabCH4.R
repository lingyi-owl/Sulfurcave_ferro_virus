# =============================================================================
# Fig. S1 — mean DNA coverage of MAGs, coloured by Lab CH4 coverage
# Figure numbering follows the manuscript (Cave_Ferro_Virus, Nov 2026).
# Split from the Fig. 1 script; data preparation is identical.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
setwd(REPO)

library(dplyr)
library(ggplot2)
##########################################
# MAG plot
##########################################

library(ggrepel)

# --- Load data ---
mean_coverage <- read.delim("metagenomics/bins_across_samples/mean_coverage.txt")
abundance     <- read.delim("metagenomics/bins_across_samples/abundance.txt")
bin_sum       <- read.delim("metagenomics/bins_summary.txt")

# --- Transform and merge ---
mean_coverage_log10 <- apply(mean_coverage[,-1], 2, log10)
colnames(mean_coverage_log10) <- paste0("log10_", colnames(mean_coverage_log10))

relab <- apply(abundance[,-1], 2, function(x) x / sum(x))
colnames(relab) <- paste0("relab_", colnames(relab))

mean_cov_all <- merge(
  cbind(mean_coverage, mean_coverage_log10),
  bin_sum, by = "bins"
)

# --- Derived variables ---
mean_cov_all$MAG_quality  <- ifelse(grepl("_MAG_", mean_cov_all$bins), "High", "Low or medium")
mean_cov_all$phylum       <- factor(ifelse(mean_cov_all$t_phylum == "", "unclassified", mean_cov_all$t_phylum))
mean_cov_all$in_culture   <- mean_cov_all$Sample_ERR10036470 > 0
mean_cov_all$border_color <- ifelse(mean_cov_all$in_culture, "black", "transparent")
mean_cov_all$size_var     <- ifelse(mean_cov_all$log10_Sample_ERR10036470 > 0,
                                    mean_cov_all$log10_Sample_ERR10036470, 0)

# --- Relative abundance & genus fill label ---
samples <- c("Sample_ERR10036468", "Sample_ERR10036469", "Sample_ERR10036470")
for (s in samples) {
  mean_cov_all[[paste0("rel_abund_", s)]] <- mean_cov_all[[s]] / sum(mean_cov_all[[s]]) * 100
}

mean_cov_all$abundant <- !(mean_cov_all$rel_abund_Sample_ERR10036468 < 1 &
                             mean_cov_all$rel_abund_Sample_ERR10036469 < 1 &
                             mean_cov_all$rel_abund_Sample_ERR10036470 < 1)

mean_cov_all$genus_fill <- ifelse(
  mean_cov_all$abundant & !is.na(mean_cov_all$t_genus) & mean_cov_all$t_genus != "",
  mean_cov_all$t_genus,
  ifelse(!mean_cov_all$abundant, "MAG below 1%", "Unclassified")
)

# --- Factor levels and colors ---
genus_levels <- c("Ferroplasma","Mycobacterium","Acidithiobacillus",
                  "Cuniculiplasma","Cutibacterium","Lawsonella",
                  "MAG below 1%","Unclassified")
mean_cov_all$genus_fill <- factor(mean_cov_all$genus_fill, levels = genus_levels)

colors_genus_fill <- c(
  "Ferroplasma"       = '#ffd92f',
  "Mycobacterium"     = '#66c2a5',
  "Acidithiobacillus" = '#fc8d62',
  "Cuniculiplasma"    = '#e5c494',
  "Cutibacterium"     = '#a6d854',
  "Lawsonella"        = '#80c1e3',
  "MAG below 1%"      = '#8da0cb',
  "Unclassified"      = '#b3b3b3'
)

# --- Labels and nudges ---
label_map <- data.frame(
  bins = c("SC_MAG_00004","SC_MAG_00005","SC_MAG_00006","SC_MAG_00008",
           "SC_MAG_00013","SC_MAG_00016","SC_MAG_00019","SC_Bin_00030"),
  marker_label = c("MAG 4","MAG 5","MAG 6","MAG 8",
                   "MAG 13","MAG 16","MAG 19","Bin 30"),
  nudge_x = c(0, 0, 0, -0.6, 0, 0, 0, 0),
  nudge_y = c(0, 0, 0, -0.01, 0, 0, 0, 0)
)
mean_cov_all <- merge(mean_cov_all, label_map, by = "bins", all.x = TRUE)
mean_cov_all$nudge_x[is.na(mean_cov_all$nudge_x)] <- 0
mean_cov_all$nudge_y[is.na(mean_cov_all$nudge_y)] <- 0

# --- Plot ---
theme_set(theme_classic(base_size = 20, base_family = "Arial"))

# Log-scale grid: major breaks at each decade (10^-2 ... 10^4), minor breaks
# log-positioned within each decade (at 2,3,...,9) so the grid tightens toward
# the top of every decade, matching the log-axis look in fig3c.
log_major <- -2:4
log_minor <- as.vector(sapply(min(log_major):(max(log_major) - 1),
                              function(k) k + log10(2:9)))
log_labels <- parse(text = paste0("10^", log_major))

# color the markers by abundance in the Lab CH4 sample
# --- handle -Inf (from log10(0)) before plotting ---
mean_cov_all$log10_Sample_ERR10036470[is.infinite(mean_cov_all$log10_Sample_ERR10036470)] <- NA

mag_plot_color_abundance <- ggplot(mean_cov_all, aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469)) +
  annotate("path",
           x = 1.3 + 0.6 * cos(seq(0, 2*pi, length.out = 100)),
           y = -1.3 + 0.6 * sin(seq(0, 2*pi, length.out = 100)),
           color = "blue", linetype = "dashed", linewidth = 1.0
  ) +
  geom_point(aes(fill = log10_Sample_ERR10036470, shape = MAG_quality),
             color = mean_cov_all$border_color, size = 8, stroke = 1.5
  ) +
  geom_point(aes(color = in_culture), size = 3, alpha = 0) +
  geom_text(aes(label = marker_label),
            hjust = -0.3, vjust = 0.5, size = 6, family = "Arial", na.rm = TRUE,
            nudge_x = mean_cov_all$nudge_x, nudge_y = mean_cov_all$nudge_y
  ) +
  annotate("segment",
           x = 0, y = 0, xend = 4, yend = 4,
           linetype = "dashed", color = "red",
           linewidth = 1.0   # increase this value for thicker line
  ) +
  scale_fill_viridis_c(
    name = "Lab CH4\nmean DNA coverage (log10)",
    na.value = "grey80"
  ) +
  scale_color_manual(
    values = c("TRUE" = "transparent", "FALSE" = "black"),
    name = "MAG in Lab CH4", labels = c("TRUE" = "Absent", "FALSE" = "Present")
  ) +
  scale_shape_manual(values = c("High" = 21, "Low or medium" = 24), name = "MAG quality") +
  labs(x = "Biofilm 1\n mean DNA coverage (log10)", y = "Biofilm 2\n  mean DNA coverage (log10)", title = "Metagenome-assembled genome") +
  theme(
    axis.text.x          = element_text(vjust = 0.5),
    panel.grid.major     = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor     = element_line(color = "grey92", linewidth = 0.5),
    panel.border         = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
    # legend.position      = "none"
  ) +
  scale_x_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  scale_y_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  coord_equal() +
  guides(
    shape = guide_legend(override.aes = list(fill = "white", color = "black", size = 6)),
    size  = guide_legend(override.aes = list(fill = "grey50", color = "black", shape = 21)),
    color = guide_legend(override.aes = list(shape = 21, fill = "grey50", size = 5, alpha = 1, stroke = 1.5))
  )

mag_plot_color_abundance