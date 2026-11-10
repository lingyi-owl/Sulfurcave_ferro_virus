# =============================================================================
# Fig. 2 — KEGG-module beeswarm (a), amino-acid KO dumbbell (b), Ferroplasma KO UpSet (c)
# Figure numbering follows the manuscript (Cave_Ferro_Virus, Nov 2026).
# Panel d (tree + ANI heatmap) has no script in this repository.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

library(dplyr)
data_dir <- file.path(REPO, "pangenomics")

# =============================================================================
# Figure 2a: KEGG/KOfam functional enrichment scatterplot
# =============================================================================

library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

result_SD_pw <-  read.delim(file.path(data_dir, "result_SD_pw.csv"), sep = ',')

result_SD_pw$is_ferroplasma <- factor(result_SD_pw$is_ferroplasma, 
                                      levels = c(TRUE, FALSE),
                                      labels = c("Ferroplasma", "Not Ferroplasma"))

f2a <- ggplot(result_SD_pw, aes(x = SD_unique_KOs, y = class)) +
  geom_beeswarm(aes(color = is_ferroplasma, size = avg_unique_KOs)) +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(
    text             = element_text(color = "black", size = 20),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.5),
    panel.border     = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    axis.text        = element_text(color = "black")
  ) +
  labs(
    size = "Average KOs",
    x = "Standard deviation of KOs",
    y = "KEGG module functional lass",
    color = "Genome"
  ) +
  
  # ── MODIFIED: Dynamic Legend Formatting & NA Translation ────────────────────
  scale_color_manual(
    values = c("Ferroplasma"     = "#ffd92f",
               "Not Ferroplasma" = "#7570b3"),
    na.value = "grey50",
    na.translate = TRUE,      # Forces ggplot to draw the missing data entry
    breaks = c("Ferroplasma", "Not Ferroplasma", NA), # Sets the explicit order
    labels = list(
      expression(italic("Ferroplasma")),
      expression("Not " * italic("Ferroplasma")), # Keeps "Not " regular and "Ferroplasma" italic
      expression("NA")
    )
  ) +
  
  # ── MODIFIED: Removed the box container and bumped label text size ──
  geom_text_repel(
    data = subset(result_SD_pw, SD_unique_KOs > 3 & is_ferroplasma == "Ferroplasma"),  
    aes(label = KEGG_Module_ACC),
    size = 5,             
    max.overlaps = Inf    
  ) 

f2a
# =============================================================================
# Figure 2b: unique amino acid KOs per KEGG module,
# Ferroplasma vs all other genomes, for one functional class.
# =============================================================================

## ---- load & subset --------------------------------------------------------
aa_df <- result_SD_pw %>%
  mutate(
    class = trimws(class),
    fc    = suppressWarnings(as.numeric(avg_unique_KOs_ferroplasma)),
    ot    = suppressWarnings(as.numeric(avg_unique_KOs_other))
  ) %>%
  filter(class == "Amino acid metabolism", !is.na(fc), !is.na(ot))

# order rows by the gap (other - Ferroplasma): largest deficit at the top
aa_df <- aa_df %>%
  mutate(gap = ot - fc) %>%
  arrange(gap) %>%
  mutate(row = row_number())

n_mod  <- nrow(aa_df)
n_lower<- sum(aa_df$fc < aa_df$ot)
mean_fc<- mean(aa_df$fc)
mean_ot<- mean(aa_df$ot)

colnames(aa_df)
keep <- c("row", "NAME", "CLASS", "PATHWAY.ID", "PATHWAY", "avg_unique_KOs_ferroplasma", "avg_unique_KOs_other", "gap")
aa_sub <- aa_df[, keep]
aa_sub <- aa_sub[order(aa_sub$row, decreasing = TRUE), ]
num_cols <- sapply(aa_sub, is.numeric)
aa_sub[num_cols] <- round(aa_sub[num_cols], 2)
write.csv(aa_sub, 
          file.path(REPO, "pangenomics/ferro_non_ferro_aa_module_ko_gap.csv"), 
          row.names = F,
          quote = F)
unique(aa_sub$PATHWAY)
## ---- colours (threaded to match the other figures) -----------------------
C_FC  <- "#ffd92f"   # Ferroplasma
C_OT  <- "#7570b3"   # all other genomes
C_SEG <- "#cfd3d9"   # connector

## ---- plot -----------------------------------------------------------------
# long form for the point legend
pts <- data.frame(
  row   = rep(aa_df$row, 2),
  value = c(aa_df$fc, aa_df$ot),
  grp   = rep(c("Ferroplasma", "All other genomes*"), each = nrow(aa_df))
)

pts$grp <- factor(pts$grp, levels = c("Ferroplasma", "All other genomes*"))

f2b <- ggplot() +
  geom_segment(data = aa_df,
               aes(x = fc, xend = ot, y = row, yend = row),
               colour = C_SEG, linewidth = 1.2) +
  geom_point(data = pts,
             aes(x = value, y = row, colour = grp), size = 3) +
  scale_colour_manual(
    values = c("Ferroplasma" = C_FC, "All other genomes*" = C_OT),
    breaks = c("Ferroplasma", "All other genomes*"),
    labels = c(expression(italic("Ferroplasma")), expression("Not " * italic("Ferroplasma"))),
    name   = "Genome"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  labs(
    x = "Mean distinct KOs per genome",
    y = sprintf("%d amino acid metabolism modules\n(sorted by gap)", n_mod)
  ) +
  theme_minimal(base_size = 20) +
  theme(
    text               = element_text(size = 20, colour = "black"),
    axis.line          = element_line(colour = "black", linewidth = 1),  # L-shaped x + y axis lines
    axis.ticks.x       = element_line(colour = "black", linewidth = 1),  # x tick marks
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 1),
    axis.text          = element_text(size = 20, colour = "black"),
    axis.title         = element_text(size = 20, colour = "black"),
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.title       = element_text(size = 20, colour = "black"),
    legend.text        = element_text(size = 20, colour = "black"),
    legend.background  = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(10, 14, 8, 10)
  )

f2b

# =============================================================================
# Figure 2c: KO-sharing UpSet plot with recalculated enrichment significance
# =============================================================================
# Run this section independently. Uses only base R and grid.
# Groups: MAG 8 (1), MAG 6 (1), non-cave Ferroplasma (10).
# Acidiplasma S12936484 is excluded. Original gene clusters and their functional
# assignments are retained; occurrence is NOT inferred from associated_groups.
# For these data, the verified anvi'o qvalue null-proportion estimate is 1,
# making its q-values identical to BH adjustment (across all 920 retained KOs).
# Optional overrides before running: fig2c_data_dir and fig2c_output_dir.

if (!exists("fig2c_data_dir")) {
  fig2c_data_dir <- if (exists("data_dir") &&
    file.exists(file.path(data_dir, "Ferro_cave_un_enriched-KOfam.txt"))) {
    data_dir
  } else {
    file.path(REPO, "pangenomics")
  }
}
if (!exists("fig2c_output_dir")) {
  fig2c_output_dir <- file.path(dirname(fig2c_data_dir), "figures", "main", "Fig2", "panels")  # Fig. 2c PDF + tables
}
dir.create(fig2c_output_dir, recursive = TRUE, showWarnings = FALSE)

fig2c_ko <- read.delim(file.path(fig2c_data_dir,
  "Ferro_cave_un_enriched-KOfam.txt"), check.names = FALSE, quote = "")
fig2c_summary_path <- file.path(fig2c_data_dir,
  "Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
fig2c_columns <- names(read.delim(fig2c_summary_path, nrows = 0, quote = ""))
fig2c_classes <- setNames(rep("NULL", length(fig2c_columns)), fig2c_columns)
fig2c_classes[c("gene_cluster_id", "genome_name")] <- "character"
fig2c_genes <- read.delim(fig2c_summary_path, colClasses = fig2c_classes,
  quote = "", check.names = FALSE)
fig2c_members <- split(fig2c_genes$genome_name, fig2c_genes$gene_cluster_id)
fig2c_mag8 <- "SFerroplasmacircular"
fig2c_mag6 <- "SMAG00006"
fig2c_references <- setdiff(unique(fig2c_genes$genome_name),
  c(fig2c_mag8, fig2c_mag6, "S12936484"))
stopifnot(length(fig2c_references) == 10L)
fig2c_occurrence <- t(vapply(fig2c_ko$gene_clusters_ids, function(ids) {
  ids <- trimws(strsplit(ids, ",", fixed = TRUE)[[1]])
  stopifnot(all(ids %in% names(fig2c_members)))
  genomes <- unique(unlist(fig2c_members[ids], use.names = FALSE))
  c(MAG8 = as.integer(fig2c_mag8 %in% genomes),
    MAG6 = as.integer(fig2c_mag6 %in% genomes),
    references = sum(fig2c_references %in% genomes))
}, numeric(3)))
fig2c_keep <- rowSums(fig2c_occurrence) > 0
fig2c_ko <- fig2c_ko[fig2c_keep, , drop = FALSE]
fig2c_occurrence <- fig2c_occurrence[fig2c_keep, , drop = FALSE]

# Refit the binomial-logit GLM / Rao score test for each distinct count pattern.
fig2c_pattern <- apply(fig2c_occurrence, 1, paste, collapse = ":")
fig2c_tests <- lapply(unique(fig2c_pattern), function(key) {
  counts <- fig2c_occurrence[match(key, fig2c_pattern), ]
  dat <- data.frame(x = as.numeric(counts), N = c(1, 1, 10),
    group = c("Ferroplasma_c", "MAG00006", "rest"))
  fit <- suppressWarnings(glm(cbind(x, N - x) ~ group,
    data = dat, family = binomial(link = "logit")))
  test <- suppressWarnings(anova(fit, test = "Rao"))
  score <- test$Rao[2]
  p <- test$`Pr(>Chi)`[2]
  if (is.na(p)) {
    stopifnot(score < 1e-3)
    p <- 1
  }
  stopifnot(score > -1e-5)
  c(score = max(score, 0), p = p)
})
fig2c_tests <- do.call(rbind, fig2c_tests)
fig2c_index <- match(fig2c_pattern, unique(fig2c_pattern))
fig2c_p <- fig2c_tests[fig2c_index, "p"]
fig2c_q <- p.adjust(fig2c_p, method = "BH")
fig2c_sig <- fig2c_p < 0.05 & fig2c_q < 0.1
fig2c_results <- data.frame(KO = fig2c_ko$accession,
  function_name = fig2c_ko$KOfam, fig2c_occurrence,
  enrichment_score = fig2c_tests[fig2c_index, "score"],
  p_value = fig2c_p, q_value = fig2c_q, significant = fig2c_sig,
  row.names = NULL, check.names = FALSE)

# Seven exact intersections; reference membership means present in >=1 reference.
fig2c_sets <- expand.grid(MAG8 = 0:1, MAG6 = 0:1, references = 0:1)
fig2c_sets <- fig2c_sets[rowSums(fig2c_sets) > 0, ]
fig2c_keys <- apply(fig2c_occurrence > 0, 1, function(x) paste(as.integer(x), collapse = ""))
fig2c_set_keys <- apply(fig2c_sets, 1, paste, collapse = "")
fig2c_sets$total <- vapply(fig2c_set_keys, function(k) sum(fig2c_keys == k), integer(1))
fig2c_sets$significant <- vapply(fig2c_set_keys,
  function(k) sum(fig2c_keys == k & fig2c_sig), integer(1))
fig2c_sets$nonsignificant <- fig2c_sets$total - fig2c_sets$significant
fig2c_sets <- fig2c_sets[order(-fig2c_sets$total), ]
rownames(fig2c_sets) <- NULL
stopifnot(sum(fig2c_sets$total) == nrow(fig2c_results))

# Native R grid grob: uniform Arial 44 pt; genus names italic only.
fig2c_grobs <- list()
fig2c_add <- function(g) fig2c_grobs[[length(fig2c_grobs) + 1L]] <<- g
fig2c_text <- function(x, y, label, colour = "#24323B",
                       centre = FALSE, rotation = 0) {
  fig2c_add(grid::textGrob(label, x, y, default.units = "native",
    just = if (centre) "centre" else "left", rot = rotation,
    gp = grid::gpar(fontsize = 44, fontfamily = "Arial",
      fontface = "plain", col = colour)))
}
fig2c_rect <- function(x, y, width, height, fill) {
  if (height <= 0) return(invisible(NULL))
  fig2c_add(grid::rectGrob(x, y, width, height, default.units = "native",
    just = c("left", "bottom"), gp = grid::gpar(fill = fill, col = NA)))
}
fig2c_line <- function(x0, y0, x1, y1, colour = "#24323B", width = 1) {
  fig2c_add(grid::segmentsGrob(x0, y0, x1, y1, default.units = "native",
    gp = grid::gpar(col = colour, lwd = width)))
}
fig2c_red <- "#FF0000"
fig2c_grey <- "#C5CDD4"
fig2c_dot_colours <- c(MAG8 = "#FFD646", MAG6 = "#FF8865", references = "#0974AC")
fig2c_max <- ceiling(max(fig2c_sets$total) / 200) * 200
for (tick in seq(0, fig2c_max, 200)) {
  y <- 320 + 850 * tick / fig2c_max
  fig2c_line(1362, y, 1370, y, "#707070", 1)
  fig2c_text(1300, y, tick, centre = TRUE)
}
fig2c_text(1225, 745, "Shared KOs between Ferroplasma group", centre = TRUE, rotation = 90)
# Visible axes for the intersection-size bars.
fig2c_line(1370, 320, 2080, 320, "#707070", 1.2)
fig2c_line(1370, 320, 1370, 1190, "#707070", 1.2)
# Alternating backgrounds for the three presence/absence matrix rows.
# Add them before circles and connecting lines so they remain behind the data.
for (j in 1:3) {
  fig2c_rect(1370, c(260, 205, 150)[j] - 27.5, 710, 55,
    if (j %% 2 == 1) "#F2F2F2" else "#FFFFFF")
}

for (i in seq_len(nrow(fig2c_sets))) {
  x <- 1393 + (i - 1) * 99
  h <- fig2c_sets$total[i] / fig2c_max * 850
  hs <- fig2c_sets$significant[i] / fig2c_max * 850
  fig2c_rect(x, 320, 62, h - hs, fig2c_grey)
  fig2c_rect(x, 320 + h - hs, 62, hs, fig2c_red)
  fig2c_text(x + 31, 320 + h + 44, fig2c_sets$total[i], centre = TRUE)
  flags <- unlist(fig2c_sets[i, c("MAG8", "MAG6", "references")])
  ys <- c(260, 205, 150)
  if (sum(flags) > 1) fig2c_line(x + 31, min(ys[flags == 1]),
    x + 31, max(ys[flags == 1]), width = 1.3)
  for (j in 1:3) fig2c_add(grid::circleGrob(x + 31, ys[j], r = 5,
    default.units = "native", gp = grid::gpar(col = NA,
      fill = if (flags[j]) fig2c_dot_colours[j] else "#E0E5E9")))
  fig2c_text(x + 31, 85, fig2c_sets$significant[i], fig2c_red, TRUE)
  fig2c_text(x + 31, 30, fig2c_sets$nonsignificant[i], "#60707D", TRUE)
}
fig2c_labels <- list(expression(italic("Ferroplasma")~"c. MAG 8"),
  expression(italic("Ferroplasma")~"MAG 6"),
  expression("Non-cave"~italic("Ferroplasma")))
for (i in 1:3) fig2c_text(650, c(260, 205, 150)[i], fig2c_labels[[i]])
fig2c_text(650, 85, "Significant", fig2c_red)
fig2c_text(650, 30, "Non-significant", "#60707D")
# Marginal totals: KO presence for genome rows; significance across all retained KOs.
# Non-cave Ferroplasma counts the union of the ten reference genomes.
# These overlapping row totals share one linear bar-length scale.
fig2c_row_totals <- c(colSums(fig2c_occurrence > 0),
  significant = sum(fig2c_sig), nonsignificant = sum(!fig2c_sig))
fig2c_row_y <- c(260, 205, 150, 85, 30)
# Reference-style marginal bars: black, leftward from a common zero baseline.
# Row labels sit to their right, between the bars and intersection matrix.
fig2c_row_max <- ceiling(max(fig2c_row_totals) / 200) * 200
fig2c_bar_zero <- 620
fig2c_bar_span <- 560
for (i in seq_along(fig2c_row_totals)) {
  bar_width <- fig2c_bar_span * fig2c_row_totals[i] / fig2c_row_max
  fig2c_rect(fig2c_bar_zero - bar_width, fig2c_row_y[i] - 6,
    bar_width, 12, "black")
}
fig2c_line(fig2c_bar_zero - fig2c_bar_span, -10, fig2c_bar_zero, -10,
  colour = "black", width = 0.7)
for (tick in seq(0, fig2c_row_max, 200)) {
  tick_x <- fig2c_bar_zero - fig2c_bar_span * tick / fig2c_row_max
  fig2c_line(tick_x, -10, tick_x, -16, colour = "black", width = 0.7)
  fig2c_text(tick_x, -45, tick, colour = "black", centre = TRUE)
}
fig2c_text(30, -104,
  "Total number of KOs per Ferroplasma group", colour = "black")
fig2c_rect(30, 976, 18, 18, fig2c_red)
fig2c_text(60, 985, sprintf("Significant (%s)", sum(fig2c_sig)))
fig2c_rect(30, 916, 18, 18, fig2c_grey)
fig2c_text(60, 925, sprintf("Non-significant (%s)", sum(!fig2c_sig)))
# Circle-color legend, aligned below the significance legend.
for (i in 1:3) {
  legend_y <- c(860, 800, 740)[i]
  fig2c_add(grid::circleGrob(39, legend_y, r = 5,
    default.units = "native", gp = grid::gpar(col = NA,
      fill = fig2c_dot_colours[i])))
  fig2c_text(60, legend_y, fig2c_labels[[i]])
}
f2c <- grid::gTree(children = do.call(grid::gList, fig2c_grobs),
  vp = grid::viewport(xscale = c(0, 2140), yscale = c(-140, 1230)), name = "Figure02c_upset")

# Cairo PDF embeds Arial fonts, including italic genus labels.
if (!capabilities("cairo")) stop("Cairo support is required to export the Arial PDF.")
grDevices::cairo_pdf(file.path(fig2c_output_dir, "Figure02c_upset_R.pdf"),
  width = 2140 / 72, height = 1370 / 72, family = "Arial", pointsize = 44)
grid::grid.newpage()
grid::grid.draw(f2c)
grDevices::dev.off()
write.csv(fig2c_results, file.path(fig2c_output_dir, "Figure02c_KO_enrichment_R.csv"), row.names = FALSE)
write.csv(fig2c_sets, file.path(fig2c_output_dir, "Figure02c_intersections_R.csv"), row.names = FALSE)
if (interactive()) {
  grid::grid.newpage()
  grid::grid.draw(f2c)
}
