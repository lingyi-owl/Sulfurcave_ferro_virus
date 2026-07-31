##########################################
# KAIJU plots
##########################################

library(dplyr)
library(ggplot2)

setwd("metagenomics/KAIJU_reads")
kaiju_reads <- list.files(pattern = "kaiju_ERR")

# ── Helper functions ───────────────────────────────────────────────────────────

# Rename sample files to readable names
rename_samples <- function(fac) {
  levels(fac)[levels(fac) == "ERR10036468.out"] <- "Biofilm 1"
  levels(fac)[levels(fac) == "ERR10036469.out"] <- "Biofilm 2"
  levels(fac)[levels(fac) == "ERR10036470.out"] <- "Lab CH4"
  as.character(fac)
}

# Label taxa below 1% only if they never exceed 1% in any sample
label_below1 <- function(df, label) {
  abundant <- df$taxon_name[df$percent >= 1]
  df$taxon_final <- df$taxon_name
  df$taxon_final[df$percent < 1 & !df$taxon_name %in% abundant] <- label
  df
}

# ── 1. Phylum ──────────────────────────────────────────────────────────────────
kaiju_phylum <- do.call(rbind, lapply(kaiju_reads[grep("phylum", kaiju_reads)], read.delim))
kaiju_phylum$Sample_name <- rename_samples(factor(kaiju_phylum$file))
kaiju_phylum <- label_below1(kaiju_phylum, "Phylum below 1%")
kaiju_phylum$taxon_final[kaiju_phylum$taxon_final == "cannot be assigned to a (non-viral) phylum"] <- "Viral"
kaiju_phylum$taxon_final[kaiju_phylum$taxon_name  == "unclassified"]                               <- "Unclassified"
kaiju_phylum$taxon_final[kaiju_phylum$taxon_final == "Candidatus Thermoplasmatota"]                <- "Thermoplasmatota"

kaiju_phylum_agg <- kaiju_phylum %>%
  group_by(Sample_name, taxon_final) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  mutate(taxon_final = factor(taxon_final, levels = c(
    "Thermoplasmatota", "Actinobacteria", "Proteobacteria",
    "Firmicutes", "Viral", "Phylum below 1%", "Unclassified"
  )))

colors_p <- c(
  "Thermoplasmatota" = '#ffd92f',
  "Actinobacteria"   = '#66c2a5',
  "Proteobacteria"   = '#fc8d62',
  "Firmicutes"       = '#a6d854',
  "Viral"            = '#e78ac3',
  "Phylum below 1%"  = '#8da0cb',
  "Unclassified"     = '#b3b3b3'
)

kaiju_phylum_ggplot <- ggplot(kaiju_phylum_agg, aes(fill = taxon_final, y = percent, x = Sample_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors_p) +
  labs(x = "", y = "% mapped reads", fill = "Phylum", title = "Phylum") +
  theme_classic(base_size = 20, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

kaiju_phylum_ggplot

# ── 2. Genus ───────────────────────────────────────────────────────────────────
kaiju_genus <- do.call(rbind, lapply(kaiju_reads[grep("genus", kaiju_reads)], read.delim))
kaiju_genus$Sample_name <- rename_samples(factor(kaiju_genus$file))
kaiju_genus <- label_below1(kaiju_genus, "Genus below 1%")
kaiju_genus$taxon_final[kaiju_genus$taxon_final == "cannot be assigned to a (non-viral) genus"] <- "Viral"
kaiju_genus$taxon_final[kaiju_genus$taxon_name  == "unclassified"]                              <- "Unclassified"

kaiju_genus_agg <- kaiju_genus %>%
  group_by(Sample_name, taxon_final) %>%
  summarise(percent = sum(percent), .groups = "drop") %>%
  mutate(taxon_final = factor(taxon_final, levels = c(
    "Ferroplasma", "Mycobacterium", "Acidithiobacillus",
    "Acidiplasma", "Viral", "Genus below 1%", "Unclassified"
  )))

colors_g <- c(
  "Ferroplasma"       = '#ffd92f',
  "Mycobacterium"     = '#66c2a5',
  "Acidithiobacillus" = '#fc8d62',
  "Acidiplasma"       = '#377eb8',
  "Viral"             = '#e78ac3',
  "Genus below 1%"    = '#8da0cb',
  "Unclassified"      = '#b3b3b3'
)

genus_labels <- c(
  expression(italic("Ferroplasma")),
  expression(italic("Mycobacterium")),
  expression(italic("Acidithiobacillus")),
  expression(italic("Acidiplasma")),
  "Viral", "Genus below 1%", "Unclassified"
)

kaiju_genus_ggplot <- ggplot(kaiju_genus_agg, aes(fill = taxon_final, y = percent, x = Sample_name)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = colors_g, labels = genus_labels) +
  labs(x = "", y = "% mapped reads", fill = "Taxonomy (genus)", title = "Genus") +
  theme_classic(base_size = 20, base_family = "Arial") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

kaiju_genus_ggplot

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

mag_plot <- ggplot(mean_cov_all, aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469)) +
  annotate("path",
           x = 1.3 + 0.6 * cos(seq(0, 2*pi, length.out = 100)),
           y = -1.3 + 0.6 * sin(seq(0, 2*pi, length.out = 100)),
           color = "blue", linetype = "dashed", linewidth = 1.0
  ) +
  geom_point(aes(fill = genus_fill, shape = MAG_quality),
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
  scale_fill_manual(
    values = colors_genus_fill, name = "Genus",
    labels = c(
      expression(italic("Ferroplasma")),
      expression(italic("Mycobacterium")),
      expression(italic("Acidithiobacillus")),
      expression(italic("Cuniculiplasma")),
      expression(italic("Cutibacterium")),
      expression(italic("Lawsonella")),
      "MAG below 1%",
      "Unclassified"
    )
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
    legend.position      = "none"
  ) +
  scale_x_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  scale_y_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  coord_equal() +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21, color = "white", size = 6)),
    shape = guide_legend(override.aes = list(fill = "white", color = "black", size = 6)),
    size  = guide_legend(override.aes = list(fill = "grey50", color = "black", shape = 21)),
    color = guide_legend(override.aes = list(shape = 21, fill = "grey50", size = 5, alpha = 1, stroke = 1.5))
  )

mag_plot

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

library(cowplot)

plot_for_color_legend <- ggplot(mean_cov_all, aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469)) +
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
    name = "MAG in Lab CH4", labels = c("TRUE" = "Absent", "FALSE" = "Present"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("High" = 21, "Low or medium" = 24), name = "MAG quality",
    guide = "none"
  ) +
  labs(x = "Biofilm 1\n mean DNA coverage (log10)", y = "Biofilm 2\n  mean DNA coverage (log10)", title = "Metagenome-assembled genome") +
  theme(
    axis.text.x          = element_text(vjust = 0.5),
    panel.grid.major     = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor     = element_line(color = "grey92", linewidth = 0.5),
    panel.border         = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
    legend.position      = "right"
  ) +
  scale_x_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  scale_y_continuous(breaks = log_major, minor_breaks = log_minor, labels = log_labels) +
  coord_equal()

plot_for_color_legend 
# legend_plot = your mag plot with legend.position = "right" and
# guide = "none" on scale_color_manual / scale_shape_manual
g <- ggplotGrob(plot_for_color_legend)
legend_only <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
grid::grid.newpage(); grid::grid.draw(legend_only)

##########################################
# Combine plots
##########################################

library(patchwork)
library(cowplot)

# --- Build a dummy data frame for MAG quality legend ---
quality_legend_df <- data.frame(
  quality = factor(c("High", "Low or medium"), levels = c("High", "Low or medium")),
  x = 1, y = 1:2
)

# --- Build MAG quality legend ---
legend_quality_combined <- get_legend(
  ggplot(quality_legend_df, aes(x, y, shape = quality)) +
    geom_point(size = 8, fill = "white", color = "grey", stroke = 1.5) +
    scale_shape_manual(
      values = c("High" = 21, "Low or medium" = 24),
      name   = "MAG quality"
    ) +
    theme_void(base_size = 20, base_family = "Arial") +
    theme(
      legend.position      = "right",
      legend.key.spacing.y = unit(6, "pt")
    )
)

# CH4 legend
# --- Build a dummy data frame for Lab CH4 legend ---
ch4_legend_df <- data.frame(
  in_culture = factor(c("Present", "Absent"), levels = c("Present", "Absent")),
  x = 1, y = 1:2
)
# --- Build Lab CH4 legend ---
legend_ch4_combined <- get_legend(
  ggplot(ch4_legend_df, aes(x, y, color = in_culture)) +
    geom_point(shape = 22, size = 8, fill = "grey50", stroke = 1.5) +
    scale_color_manual(
      values = c("Present" = "black", "Absent" = "transparent"),
      name   = "MAG in Lab CH4"
    ) +
    theme_void(base_size = 20, base_family = "Arial") +
    theme(
      legend.position      = "right",
      legend.key.spacing.y = unit(6, "pt")
    )
)

# phylum legend
# --- Build a dummy data frame for phylum legend ---
phylum_legend_df <- data.frame(
  phylum = factor(names(colors_p), levels = names(colors_p)),
  x = 1, y = seq_along(colors_p)
)

# --- Build phylum legend using same square style ---
legend_phylum_combined <- get_legend(
  ggplot(phylum_legend_df, aes(x, y, fill = phylum)) +
    geom_point(shape = 22, size = 8, color = "white") +
    scale_fill_manual(
      values = colors_p,
      name   = "Phylum",
      labels = c(
        "Thermoplasmatota",
        "Actinobacteria",
        "Proteobacteria",
        "Firmicutes",
        "Viral",
        "Phylum below 1%",
        "Unclassified"
      )
    ) +
    theme_void(base_size = 20, base_family = "Arial") +
    theme(
      legend.position      = "right",
      legend.key.spacing.y = unit(-4, "pt")
    )
)

# genus legend
# --- Combined color map for shared genus legend ---
colors_genus_combined <- c(
  "Ferroplasma"       = '#ffd92f',
  "Mycobacterium"     = '#66c2a5',
  "Acidithiobacillus" = '#fc8d62',
  "Cuniculiplasma"    = '#e5c494',   # from mag_plot only
  "Cutibacterium"     = '#a6d854',   # from mag_plot only
  "Lawsonella"        = '#80c1e3',   # from mag_plot only
  "Acidiplasma"       = '#377eb8',   # from kaiju only
  "Viral"             = '#e78ac3',   # from kaiju only
  "Below 1% / MAG below 1%" = '#8da0cb',
  "Unclassified"      = '#b3b3b3'
)

# --- Build a dummy data frame covering all genus levels ---
genus_combined_df <- data.frame(
  genus = factor(names(colors_genus_combined), levels = names(colors_genus_combined)),
  x = 1, y = seq_along(colors_genus_combined)
)

# --- Build legend using kaiju_genus_ggplot point shape (filled circle = 16) ---
legend_genus_combined <- get_legend(
  ggplot(genus_combined_df, aes(x, y, fill = genus)) +
    geom_point(shape = 22, size = 8, color = "white") +
    scale_fill_manual(
      values = colors_genus_combined,
      name   = "Genus",
      labels = c(
        expression(italic("Ferroplasma")),
        expression(italic("Mycobacterium")),
        expression(italic("Acidithiobacillus")),
        expression(italic("Cuniculiplasma")),
        expression(italic("Cutibacterium")),
        expression(italic("Lawsonella")),
        expression(italic("Acidiplasma")),
        "Viral",
        "Genus/MAG below 1%",
        "Unclassified"
      )
    ) +
    theme_void(base_size = 20, base_family = "Arial") +
    theme(
      legend.position    = "right",
      legend.key.spacing.y = unit(-4, "pt")  # negative value reduces spacing
    )
)

# --- Stack legends ---
legends_combined <- plot_grid(
  legend_phylum_combined,
  legend_genus_combined,
  legend_quality_combined,
  legend_ch4_combined,
  ncol = 1,
  align = "v",
  rel_heights = c(1.5, 2.5, 1, 1)
)

legends_combined
# --- Combine plots (no legends) ---

plots_ab <- (
  (kaiju_phylum_ggplot + theme(legend.position = "none")) +
    (kaiju_genus_ggplot  + theme(legend.position = "none"))
)
plots_ab

leftbottom_mags <- subset(mean_cov_all,
                          log10_Sample_ERR10036468 >= 0.5 & log10_Sample_ERR10036468 <= 2 &
                            log10_Sample_ERR10036469 >= -2   & log10_Sample_ERR10036469 <= -1)
