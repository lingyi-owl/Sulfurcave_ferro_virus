# =====================================================================
# General physiochemical feature heatmap: Cave vs IMG/VR effect sizes.
# Reproduces general_features_heatmap.pdf, but from the recomputed raw
# feature table (per-protein values) rather than a precomputed diff table.
#
# For each feature x protein type this script computes, from the raw values:
#   - effect size = Cohen's d (mean_cave - mean_reference)/pooled SD, and
#   - one-tailed permutation p, using the SAME scheme as FigS9
#     (physiochemical_features_paired): 1000 iterations, each drawing 90% of
#     the cave n with replacement from BOTH groups; p = fraction of per-iter
#     differences (cave - reference) that fell on the null side of zero
#     (resolution 1/1000 -> p < 0.001). seed = 42.
#
# Restricted to the 8 features shown in FigS9, and axis labels are matched
# to that figure (protein names + feature titles with units).
#
# Input : ../viromics/cave_imgvr_features_raw.csv
# Output: figures/supplementary/FigS12/FigS12_physicochemical_effect_size_heatmap.svg (Fig. S12)
# =====================================================================
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_S12 <- file.path(REPO, "figures/supplementary/FigS12")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_S12, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(ComplexHeatmap); library(circlize); library(dplyr)
                  library(tidyr); library(grid); library(data.table)})

data_csv <- file.path(REPO, "viromics/cave_imgvr_features_raw.csv")
out_svg  <- file.path(FIGOUT_S12, "FigS12_physicochemical_effect_size_heatmap.svg")

feat <- fread(data_csv)
# Reproduce FigS9's exact feature transforms so the resampled values are
# literally identical (Aromaticity = F+W+Y+H %, secondary-structure -> %).
# Scaling by a positive constant leaves both the one-tailed p and Cohen's d
# unchanged, but this keeps the resampled data byte-for-byte the same as FigS9.
feat[, Aromaticity := (F_Percent + W_Percent + Y_Percent + H_Percent) / 100]
feat[, `:=`(HelixFraction = HelixFraction * 100,
            SheetFraction = SheetFraction * 100,
            TurnFraction  = TurnFraction  * 100)]

# --- features kept as rows (+ FigS9 y-axis titles) -------------------
feature_labels <- c(
  ProteinLength    = "Length (aa)",
  IsoelectricPoint = "Isoelectric point (pH units)",
  FlexibilityMean  = "Mean flexibility index (a.u.)",
  Aromaticity      = "Aromatic content, F+W+Y+H (%)",
  Hydrophobicity   = "Hydrophobicity, GRAVY (a.u.)",
  HelixFraction    = "\u03b1-Helix fraction (%)",
  SheetFraction    = "\u03b2-Sheet fraction (%)",
  TurnFraction     = "Turn fraction (%)"
)
features <- names(feature_labels)

protein_order <- c("major_head_protein","minor_head_protein","major_tail_protein",
                   "minor_tail_protein","tail_fiber_protein","tail_spike_protein",
                   "dna_polymerase","endolysin","head_maturation_protease",
                   "integrase","terminase_large_subunit","terminase_small_subunit")
nice <- function(x){ s<-gsub("_"," ",x); s<-paste0(toupper(substr(s,1,1)),substr(s,2,nchar(s))); sub("^Dna ","DNA ",s) }

# --- effect size + one-tailed permutation p ---------------------------
# The permutation p reuses FigS9's EXACT resampling (same seed, same nested
# loop order over features x protein_order, same resample_means drawing all
# N cave means then all N reference means). Because the RNG stream, seed,
# data, and call order are identical, the per-iteration difference
# distributions -- and hence p_one_tailed -- match FigS9 cell for cell.
# Effect size is Cohen's d on the (raw) values: deterministic, no RNG.
set.seed(42); N_ITER <- 1000; FRAC <- 0.90

resample_means <- function(cave, ref) {           # identical to FigS9
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

mat   <- matrix(NA_real_, length(features), length(protein_order), dimnames=list(features, protein_order))
p_mat <- matrix(NA_real_, length(features), length(protein_order), dimnames=list(features, protein_order))
for (fc in features) for (pr in protein_order){
  cv <- feat[origin=="sulfurcave" & protein==pr, get(fc)]; cv <- cv[is.finite(cv)]
  rf <- feat[origin=="imgvr"      & protein==pr, get(fc)]; rf <- rf[is.finite(rf)]
  if (length(cv) < 2 || length(rf) < 2) next
  m   <- resample_means(cv, rf)           # draws first -> preserves FigS9 RNG order
  d   <- m$cave - m$ref; obs <- mean(d)
  p_mat[fc, pr] <- if (obs > 0) mean(d <= 0) else mean(d >= 0)
  mat[fc, pr]   <- cohen_d(cv, rf)
}
rownames(mat) <- feature_labels[rownames(mat)]; rownames(p_mat) <- rownames(mat)
colnames(mat) <- nice(colnames(mat));           colnames(p_mat) <- colnames(mat)

# --- heatmap (styling as in the original) ----------------------------
col_fun <- colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
dynamic_split <- factor(ifelse(protein_order %in% protein_order[1:6], "Structural", "Non-Structural"),
                        levels = c("Structural", "Non-Structural"))

BASE_FAMILY <- "Arial"; BASE_SIZE <- 20
hp_physio <- Heatmap(mat,
  name = "Effect Size", col = col_fun,
  width  = unit(22, "cm"), height = unit(13, "cm"),   # fixed body -> room for size-20 labels
  column_title = c("", ""),                     # no figure title (and no slice titles)
  cluster_columns = TRUE, cluster_column_slices = FALSE,
  cluster_rows = FALSE, show_column_dend = TRUE,   # keep FigS9 feature order (no row clustering)
  column_split = dynamic_split,
  column_names_rot = 45, column_names_side = "bottom",
  column_names_gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY),
  row_names_gp = gpar(fontsize = BASE_SIZE, fontfamily = BASE_FAMILY),
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

svglite::svglite(out_svg, width = 18, height = 10)
draw(hp_physio, heatmap_legend_side = "bottom", padding = unit(c(6, 8, 6, 8), "mm"))
dev.off()
