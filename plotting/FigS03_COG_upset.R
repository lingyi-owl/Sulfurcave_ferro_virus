# =============================================================================
# Fig. S3: UpSet plot of COG intersections (cave pangenome)
# =============================================================================
# Builds a genome x gene_cluster_id presence/absence matrix from the anvi'o
# pangenome summary table, then plots set/intersection sizes with UpSetR.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

library(dplyr)
library(UpSetR)
library(ggplotify)
library(grid)
library(cowplot)

# -----------------------------------------------------------------------------
# 1. Load pangenome summary table
# -----------------------------------------------------------------------------
data_dir <- file.path(REPO, "pangenomics")

CAVE <- read.delim(file.path(data_dir, "CAVE_R_Pan_gene_clusters_summary.txt.gz"))

# -----------------------------------------------------------------------------
# 2. Build genome x gene_cluster_id presence/absence matrix
# -----------------------------------------------------------------------------
genomes <- unique(CAVE$genome_name)
gene_clusters <- unique(CAVE$gene_cluster_id)

# For each genome, the set of gene cluster IDs present in it
OG_list <- lapply(genomes, function(g) {
  unique(CAVE$gene_cluster_id[CAVE$genome_name == g])
})
names(OG_list) <- genomes

# Binary matrix: rows = genomes, columns = gene clusters
OG_matrix <- t(sapply(OG_list, function(present) as.integer(gene_clusters %in% present)))
colnames(OG_matrix) <- gene_clusters
rownames(OG_matrix) <- genomes

# Cache to disk (uncomment to regenerate; otherwise read cached version below)
# write.csv(OG_matrix, file.path(data_dir, "OG_matrix_cave.csv"))
OG_matrix <- read.csv(file.path(data_dir, "OG_matrix_cave.csv"), row.names = "X")

# -----------------------------------------------------------------------------
# 3. Relabel genomes of interest for plotting
# -----------------------------------------------------------------------------
# Sanity check before relabeling by position - confirm row order matches
# the expected genome identities:
#   [2]  -> "F. c."                  (Ferroplasma circular contig)
#   [8]  -> "F. MAG6"                (Ferroplasma MAG6)
#   [19] -> "M. MAG2"                (Mycobacterium MAG2)
#   [20] -> "M.MAG3"                 (Mycobacterium MAG3)
#   [21] -> "C. M. methanotrophicum"
stopifnot(length(rownames(OG_matrix)) >= 21)
print(rownames(OG_matrix)[c(2, 8, 19, 20, 21)])  # verify before relabeling

rownames(OG_matrix)[c(2, 8, 19, 20, 21)] <- c(
  "F. c.", "F. MAG6", "M. MAG2", "M.MAG3", "C. M. methanotrophicum"
)

# -----------------------------------------------------------------------------
# 4. UpSet plot figS2
# -----------------------------------------------------------------------------
sets_of_interest <- c("F. c.", "F. MAG6", "M. MAG2", "M.MAG3", "C. M. methanotrophicum")

# Bar colors, in plotting order (order.by = "freq"):
#   teal   (#66c2a5) = Mycobacterium-only intersections
#   yellow (#ffd92f) = Ferroplasma-only intersections
#   gray   (#a3a3a3) = mixed Ferroplasma + Mycobacterium intersections
bar_colors <- c(
  "#66c2a5", "#66c2a5", "#66c2a5", "#ffd92f", "#66c2a5",
  "#ffd92f", "#ffd92f", "#66c2a5", "#66c2a5", "#66c2a5",
  "#a3a3a3", "#a3a3a3", "#a3a3a3", "#a3a3a3", "#a3a3a3",
  "#a3a3a3", "#a3a3a3", "#a3a3a3", "#a3a3a3", "#a3a3a3"
)

# NOTE on font/size: UpSetR is built on grid/base graphics, not a single
# ggplot object, so there is no one-line ggplot theme() equivalent like
# theme_classic(base_size = 20, base_family = "Arial"). Two separate settings
# stand in for it:
#   - par(family = ...) sets the font family on the graphics device, used by
#     all text UpSetR draws (axis titles, labels, bar counts, set names).
#   - text.scale inside upset() scales all text sizes together (not an exact
#     point size like base_size=20; treat it as a relative multiplier -
#     2-2.2 is roughly equivalent to size-20-ish text for this plot's bar
#     count, so left at 2 here, adjust if it still looks small/large).
par(family = "Arial")

upset_plot <- upset(
  as.data.frame(t(OG_matrix)),
  main.bar.color = bar_colors,
  sets.bar.color = "black",
  order.by = "freq",
  nintersects = NA,
  nsets = length(sets_of_interest),
  sets = sets_of_interest,
  keep.order = TRUE,     # <- keeps Set Size bars in 'sets' order instead of
  #    being dropped/resorted; restores the bottom-left
  #    "Set Size" side bars matching the 5 named genomes
  text.scale = 2,
  mainbar.y.label = "shared COGs between MAGs",
  sets.x.label = "total number of COGs per MAG"
)

figS2 <- upset_plot
figS2

# -----------------------------------------------------------------------------
# 5. Composite a manual color legend into the top-left corner
# -----------------------------------------------------------------------------
# UpSetR has no native legend for main.bar.color, so it's added here as a
# separate small ggplot legend, overlaid on the upset plot with cowplot.
figS2_legend_df <- data.frame(
  label = factor(
    c("Ferroplasma unique", "Mycobacterium unique", "Shared"),
    levels = c("Ferroplasma unique", "Mycobacterium unique", "Shared")
  ),
  x = 0, y = 0
)

fig2a_legend_colors <- c(
  "Ferroplasma unique"   = "#ffd92f",
  "Mycobacterium unique" = "#66c2a5",
  "Shared"  = "#a3a3a3"
)

fig2a_legend_plot <- ggplot(legend_df, aes(x, y, color = label)) +
  geom_point(size = 4, shape = 15) +
  scale_color_manual(values = legend_colors, name = NULL) +
  theme_void(base_size = 20, base_family = "Arial") +
  theme(
    legend.position = "left",
    legend.justification = "left",
    legend.text = element_text(size = 20),
    legend.key = element_blank()
  )

figS2_legend_plot
