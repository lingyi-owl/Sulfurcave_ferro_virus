# =============================================================================
# rTCA / TCA marker-protein heatmap  --  ComplexHeatmap version
#
# Same data and row/column semantics as analysis/marker_protein_heatmap_data.R (the ggplot2
# supplementary version); this script only swaps the renderer to
# ComplexHeatmap so the columns can be split by pathway block and the grey
# "gene present / protein not detected" state kept distinct from white
# "gene absent".
#
# Cell convention:
#   colour (viridis) = marker protein DETECTED & quantified in the cave biofilm
#                      metaproteome (average log10 LFQ over the two replicates)
#   grey             = marker KO gene PRESENT in the genome, protein NOT detected
#   white            = genome does not carry that marker KO
#
# It re-uses the data-building half of analysis/marker_protein_heatmap_data.R (sourced below),
# which produces `df_cell` (tax_id, KO_role factor, avg_LFQ, status, row_label
# factor with markdown italics). Set `repo` if the path differs.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

suppressPackageStartupMessages({
  library(ComplexHeatmap); library(circlize)
  library(grid); library(gridtext); library(viridisLite); library(svglite)
})

repo <- REPO
outdir <- file.path(repo, "plotting")
FIGOUT_S04 <- file.path(repo, "figures/supplementary/FigS04")  # figure output folder
dir.create(FIGOUT_S04, recursive = TRUE, showWarnings = FALSE)

# ---- global text style: Arial, size 20 -------------------------------------
FONT <- "Arial"; FS <- 20

# ---- build df_cell / row_label / KO_role from the ggplot2 script ------------
# (source the pipeline; it also writes the ggplot2 PDF/PNG, which is harmless)
source(file.path(repo, "analysis", "marker_protein_heatmap_data.R"), echo = FALSE)

# ---- value + status matrices (rows top->bottom = MAG 1 ... MAG 33) ----------
rl_top <- rev(levels(df_cell$row_label))      # ggplot puts level[1] at bottom; reverse for top-down
cl     <- levels(df_cell$KO_role)
mat  <- matrix(NA_real_,  length(rl_top), length(cl), dimnames = list(rl_top, cl))
stat <- matrix("absent",  length(rl_top), length(cl), dimnames = list(rl_top, cl))
for (i in seq_len(nrow(df_cell))) {
  a <- as.character(df_cell$row_label[i]); b <- as.character(df_cell$KO_role[i])
  mat[a, b]  <- df_cell$avg_LFQ[i]
  stat[a, b] <- df_cell$status[i]             # measured / not_detected
}

# ---- row grouping by genus --------------------------------------------------
# genus of each row (via tax_id -> MAG number -> Table 2 genus); blank -> "NA"
row_tax <- sapply(rl_top, function(lbl) df_cell$tax_id[match(lbl, as.character(df_cell$row_label))])
row_gen <- sapply(row_tax, function(t){ g <- genus_by_mag[[as.character(mag_num(t))]]
                                        if (is.null(g) || g == "") "NA" else g })
# slice order: ascending minimum MAG number within each genus, unassigned "NA" last
mn  <- tapply(sapply(row_tax, mag_num), row_gen, min)
gen_levels <- names(sort(mn)); gen_levels <- c(setdiff(gen_levels, "NA"), if ("NA" %in% gen_levels) "NA")
row_split <- factor(unname(row_gen), levels = gen_levels)
# slice titles: italic for the genera Table 2 italicises
# genus group titles: italic; bold-italic if the genus is cited in the rTCA manuscript text
# (only Ferroplasma is named in the rTCA paragraph)
mentioned_genera <- c("Ferroplasma")
gen_title <- vapply(gen_levels, function(g) {
  if (g == "NA") return("NA")
  if (g %in% mentioned_genera) paste0("***", g, "***") else paste0("*", g, "*")
}, character(1))

# ---- column split by pathway block + KO-id column labels --------------------
col_block <- factor(sub(" \\(K.*$", "", cl), levels = unique(sub(" \\(K.*$", "", cl)))
col_ko    <- sub("^.*\\((K[0-9]+)\\)$", "\\1", cl)
# KO IDs cited in the manuscript text -> drawn bold on the x-axis
mentioned_ko <- c("K01006","K01007","K00024","K00031","K01681","K01902","K01903","K00174","K00175",
                  "K03841","K01624","K00600","K00058","K00831","K16157","K16158","K16159","K16161",
                  "K08093","K08094","K21071","K01595","K01689")
ko_face <- ifelse(col_ko %in% mentioned_ko, "bold", "plain")

# ---- markdown row labels (mixed italics) via gridtext -----------------------
row_md <- gt_render(rl_top)

# ---- colour ramp: viridis over the observed LFQ range -----------------------
lfq_max <- max(mat, na.rm = TRUE)
col_fun <- colorRamp2(seq(0, lfq_max, length.out = 5), viridis(5))

# ---- three-state cell painter ----------------------------------------------
cell_fun <- function(j, i, x, y, w, h, fill) {
  s <- stat[i, j]
  fillcol <- if (s == "measured") col_fun(mat[i, j]) else if (s == "not_detected") "grey80" else "white"
  grid.rect(x, y, w, h, gp = gpar(fill = fillcol, col = "white", lwd = 0.6))
}

ttl_gp <- gpar(fontfamily = FONT, fontsize = FS, fontface = "plain", col = "black")
lab_gp <- gpar(fontfamily = FONT, fontsize = FS)

lgd_val  <- Legend(col_fun = col_fun, title = "Mean protein intensity (LFQ; log 10)",
                   at = round(seq(0, lfq_max, length.out = 5), 1),
                   title_gp = ttl_gp, labels_gp = lab_gp,
                   grid_width = unit(6, "mm"), legend_height = unit(45, "mm"))
lgd_grey <- Legend(labels = c("Yes", "No"), title = "Gene presence",
                   legend_gp = gpar(fill = c("grey80", "white")), border = "grey40",
                   title_gp = ttl_gp, labels_gp = lab_gp,
                   grid_height = unit(8, "mm"), grid_width = unit(8, "mm"),
                   row_gap = unit(6, "mm"))   # extra space between the Yes / No boxes

# row labels: genus portion italic (already encoded as *...* markdown in rl_top),
# everything else regular; render at Arial 20 via gridtext
row_md <- gt_render(rl_top, gp = gpar(fontfamily = FONT, fontsize = FS))
# column-block titles at the BOTTOM, 0-degree rotation, Arial 20
ht <- Heatmap(
  mat, name = "LFQ", col = col_fun, na_col = "white",
  rect_gp = gpar(type = "none"), cell_fun = cell_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_order = rl_top, column_order = cl,
  row_split = row_split, cluster_row_slices = FALSE,
  row_title = gt_render(gen_title, gp = gpar(fontfamily = FONT, fontsize = FS)),
  row_title_rot = 0, row_title_side = "right", row_gap = unit(1.2, "mm"),
  column_split = col_block, cluster_column_slices = FALSE,
  column_title = levels(col_block), column_title_side = "bottom",
  column_title_gp = gpar(fontfamily = FONT, fontsize = FS), column_title_rot = 0,
  column_labels = col_ko,
  column_names_gp = gpar(fontfamily = FONT, fontsize = FS, fontface = ko_face), column_names_rot = 45,
  row_labels = row_md,
  row_names_gp = gpar(fontfamily = FONT, fontsize = FS), row_names_side = "left",
  show_heatmap_legend = FALSE, border = TRUE, column_gap = unit(1.5, "mm"))

# extra left padding so the long dual reference-genome labels are not clipped
draw_it <- function() draw(ht, annotation_legend_list = list(lgd_val, lgd_grey),
                           annotation_legend_side = "right", legend_gap = unit(10, "mm"),
                           padding = unit(c(2, 165, 2, 2), "mm"))  # bottom, left, top, right

W <- 30; H <- 13
cairo_pdf(file.path(FIGOUT_S04, "FigS04_rTCA_heatmap.pdf"), width = W, height = H); draw_it(); dev.off()
png(file.path(FIGOUT_S04, "FigS04_rTCA_heatmap.png"), width = W, height = H, units = "in", res = 300, type = "cairo"); draw_it(); dev.off()
svglite(file.path(FIGOUT_S04, "FigS04_rTCA_heatmap.svg"), width = W, height = H); draw_it(); dev.off()
