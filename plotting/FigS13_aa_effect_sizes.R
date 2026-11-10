# =====================================================================
# Amino-acid composition heatmap: Cave vs IMG/VR effect sizes.
# Remakes aa_heatmap.pdf from the recomputed raw feature table (the same
# data used for FigS10_aa_percentage_radar_paired), instead of a
# precomputed diff table.
#
# For each amino acid x protein type this script computes, from the raw
# per-protein %(AA) values:
#   - effect size = Cohen's d (mean_cave - mean_reference)/pooled SD, and
#   - one-tailed permutation p, using the EXACT scheme of FigS9/FigS10/FigS11:
#     set.seed(42), 1000 iterations, each drawing n = floor(0.90 * cave_count)
#     WITH replacement from BOTH the cave and reference pools; statistic =
#     per-iteration difference of subsample means (cave - reference);
#     p = fraction of iterations on the side OPPOSITE the observed mean
#     (resolution 1/1000 -> p < 0.001).
#
# Rows (20 amino acids) and columns (12 protein types) are both clustered,
# columns split Structural / Non-Structural -- matching the attached figure.
#
# Input : ../viromics/cave_imgvr_features_raw.csv
# Output: figures/supplementary/FigS13/FigS13_aa_effect_size_heatmap.svg (Fig. S13)
# =====================================================================
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_S13 <- file.path(REPO, "figures/supplementary/FigS13")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_S13, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(ComplexHeatmap); library(circlize); library(grid); library(data.table)})

data_csv <- file.path(REPO, "viromics/cave_imgvr_features_raw.csv")
out_svg  <- file.path(FIGOUT_S13, "FigS13_aa_effect_size_heatmap.svg")

feat <- fread(data_csv)
aa_cols <- grep("^[A-Z]_Percent$", names(feat), value = TRUE)   # the 20 single-AA percent columns

protein_order <- c("major_head_protein","minor_head_protein","major_tail_protein",
                   "minor_tail_protein","tail_fiber_protein","tail_spike_protein",
                   "dna_polymerase","endolysin","head_maturation_protease",
                   "integrase","terminase_large_subunit","terminase_small_subunit")
nice <- function(x){ s<-gsub("_"," ",x); s<-paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s))); sub("^Dna ","DNA ",s) }

# --- effect size (Cohen's d) + one-tailed permutation p (FigS9 method) ---
set.seed(42); N_ITER <- 1000; FRAC <- 0.90
resample_means <- function(cave, ref) {                 # identical to FigS9/FigS10/FigS11
  n <- floor(FRAC * length(cave))
  list(cave = vapply(seq_len(N_ITER), function(i) mean(sample(cave, n, replace = TRUE)), numeric(1)),
       ref  = vapply(seq_len(N_ITER), function(i) mean(sample(ref,  n, replace = TRUE)), numeric(1)))
}
cohen_d <- function(a, b){
  na<-length(a); nb<-length(b)
  sp <- sqrt(((na-1)*var(a) + (nb-1)*var(b)) / (na+nb-2))
  if (is.na(sp) || sp==0) return(NA_real_)
  (mean(a)-mean(b))/sp
}

mat   <- matrix(NA_real_, length(aa_cols), length(protein_order), dimnames=list(aa_cols, protein_order))
p_mat <- matrix(NA_real_, length(aa_cols), length(protein_order), dimnames=list(aa_cols, protein_order))
for (fc in aa_cols) for (pr in protein_order){
  cv <- feat[origin=="sulfurcave" & protein==pr, get(fc)]; cv <- cv[is.finite(cv)]
  rf <- feat[origin=="imgvr"      & protein==pr, get(fc)]; rf <- rf[is.finite(rf)]
  if (length(cv) < 2 || length(rf) < 2) next
  m   <- resample_means(cv, rf)
  d   <- m$cave - m$ref; obs <- mean(d)
  p_mat[fc, pr] <- if (obs > 0) mean(d <= 0) else mean(d >= 0)
  mat[fc, pr]   <- cohen_d(cv, rf)
}
colnames(mat) <- nice(colnames(mat)); colnames(p_mat) <- colnames(mat)

# rows in FigS10 order (grouped by amino-acid physicochemical property, from A),
# labelled with the single-letter amino-acid code
aa_order <- c("A","V","L","I","M","F","W","Y","H","R","K","D","E","S","T","N","Q","G","P","C")
ord <- paste0(aa_order, "_Percent")
mat   <- mat[ord, , drop = FALSE]; p_mat <- p_mat[ord, , drop = FALSE]
rownames(mat) <- aa_order;         rownames(p_mat) <- aa_order

# --- heatmap (styling as in the attached original) -------------------
col_fun <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
dynamic_split <- factor(ifelse(protein_order %in% protein_order[1:6], "Structural", "Non-Structural"),
                        levels = c("Structural", "Non-Structural"))
BASE_FAMILY <- "Arial"; BASE_SIZE <- 20

# amino-acid property groups: order + Paul-Tol colours as in FigS10
aa_groups  <- list(Aliphatic=c("A","V","L","I","M"), Aromatic=c("F","W","Y","H"),
                   Positive=c("R","K"), Negative=c("D","E"),
                   Polar=c("S","T","N","Q"), Special=c("G","P","C"))
grp_colors <- c(Aliphatic="#117733", Aromatic="#332288", Positive="#88CCEE",
                Negative="#AA4499", Polar="#44AA99", Special="#DDCC77")
row_grp <- factor(rep(names(aa_groups), lengths(aa_groups)), levels = names(aa_groups))

# right side (left->right): AA letters, then colour bar per group, then group
# name (coloured to match), like FigS10. Letters live in the annotation so the
# group bars/labels sit to the RIGHT of the amino-acid letters.
right_ann <- rowAnnotation(
  aa   = anno_text(rownames(mat), rot = 0, just = "centre", location = unit(0.5, "npc"),
                   gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY)),
  bar  = anno_block(gp = gpar(fill = grp_colors, col = "white"), width = unit(0.55, "cm")),
  labs = anno_block(gp = gpar(fill = NA, col = NA),
                    labels = names(aa_groups), labels_rot = 0,
                    labels_gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY, col = grp_colors),
                    width = unit(4.6, "cm")),
  show_annotation_name = FALSE, gap = unit(2, "mm"))

hp_aa <- Heatmap(mat,
  name = "Effect Size", col = col_fun,
  width  = unit(20, "cm"), height = unit(24, "cm"),   # fixed body -> room for size-20 labels
  column_title = c("", ""),                            # no title (and no slice titles)
  cluster_columns = TRUE, cluster_column_slices = FALSE,
  cluster_rows = FALSE, cluster_row_slices = FALSE, show_column_dend = TRUE,  # FigS10 row order
  row_split = row_grp, row_title = NULL, row_gap = unit(1.5, "mm"),          # group by property, no left titles
  column_split = dynamic_split,
  column_names_rot = 45, column_names_side = "bottom",
  column_names_gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY),
  show_row_names = FALSE,                 # letters drawn in the right annotation instead
  right_annotation = right_ann,
  heatmap_legend_param = list(direction = "horizontal", legend_width = unit(5, "cm"),
                              title_gp  = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY),
                              labels_gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY)),
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = function(j, i, x, y, width, height, fill) {
    p <- p_mat[i, j]
    if (!is.na(p) && p < 0.05) {
      lbl <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))
      grid.text(lbl, x, y, gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY))
    }
  })

svglite::svglite(out_svg, width = 15, height = 14)
draw(hp_aa, heatmap_legend_side = "bottom", padding = unit(c(6, 8, 6, 8), "mm"))
dev.off()
