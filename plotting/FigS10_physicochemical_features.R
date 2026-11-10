# =====================================================================
# Physiochemical features: PAIRED resampled violins (Cave vs IMG/VR),
# one panel per feature, 12 phage protein types per panel.
#
# PERMUTATION (identical scheme to plotting/
# FigS14_myco_TM_pI.R): for each feature x protein type,
# 1000 iterations; each iteration draws n = floor(0.90 * cave_count) proteins
# WITH replacement from Cave and the same n from IMG/VR, and records the two
# subsample means. Violins/beeswarm show those per-iteration means. The
# one-tailed p is taken from the per-iteration difference (cave - reference):
# frac = P(delta <= 0) if the observed mean delta > 0 else P(delta >= 0);
# reported as "p < 0.001" when frac < 1/1000. seed = 42.
#
# Input : ../viromics/cave_imgvr_features_raw.csv  (origin, protein, features +
#         F/W/Y/H percent columns; *_Percent stored as fraction x 10000).
#         Cave = all proteins; IMG/VR = seeded 3000/protein-type subsample.
# Output: physiochemical_features_paired.svg
# House style: Arial, size 24, black text, house colours.
# =====================================================================
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
FIGOUT_S10 <- file.path(REPO, "figures/supplementary/FigS10")  # figure output folder (manuscript numbering)
dir.create(FIGOUT_S10, recursive = TRUE, showWarnings = FALSE)

suppressMessages({library(data.table); library(ggplot2); library(ggbeeswarm); library(patchwork)})

data_csv <- file.path(REPO, "viromics/cave_imgvr_features_raw.csv")
feat <- fread(data_csv)

# Aromaticity = aromatic content including histidine: (F+W+Y+H)/length, as a percentage.
# The *_Percent columns are fraction x 10000, so /100 gives percent.
feat[, Aromaticity := (F_Percent + W_Percent + Y_Percent + H_Percent) / 100]
# Secondary-structure fractions expressed as percentages (0-1 -> %).
feat[, `:=`(HelixFraction = HelixFraction * 100,
            SheetFraction = SheetFraction * 100,
            TurnFraction  = TurnFraction  * 100)]

# --- Features to keep (in panel order) + y-axis titles ---------------
features <- c(
  ProteinLength    = "Length (aa)",
  IsoelectricPoint = "Isoelectric point (pH units)",
  FlexibilityMean  = "Mean flexibility index (a.u.)",
  Aromaticity      = "Aromatic content, F+W+Y+H (%)",
  Hydrophobicity   = "Hydrophobicity, GRAVY (a.u.)",
  HelixFraction    = "\u03b1-Helix fraction (%)",
  SheetFraction    = "\u03b2-Sheet fraction (%)",
  TurnFraction     = "Turn fraction (%)"
)
protein_order <- c("major_head_protein","minor_head_protein","major_tail_protein",
                   "minor_tail_protein","tail_fiber_protein","tail_spike_protein",
                   "dna_polymerase","endolysin","head_maturation_protease",
                   "integrase","terminase_large_subunit","terminase_small_subunit")

# --- Resampling (same scheme as the reference script) ----------------
set.seed(42); N_ITER <- 1000; FRAC <- 0.90; P_RES <- 1/N_ITER; DISPLAY_PTS <- 200
fmt_p <- function(p) {
  base <- ifelse(p < P_RES, "p < 0.001", sprintf("p = %.3f", p))
  ifelse(p < 0.05, paste0(base, " *"), base)   # star significant comparisons
}

resample_means <- function(cave, ref) {
  n <- floor(FRAC * length(cave))
  data.table(cave = vapply(seq_len(N_ITER), function(i) mean(sample(cave, n, replace = TRUE)), numeric(1)),
             reference = vapply(seq_len(N_ITER), function(i) mean(sample(ref, n, replace = TRUE)), numeric(1)))
}

means_list <- list(); sig_list <- list()
for (fc in names(features)) for (pr in protein_order) {
  cv <- feat[origin == "sulfurcave" & protein == pr, get(fc)]
  rf <- feat[origin == "imgvr"      & protein == pr, get(fc)]
  cv <- cv[is.finite(cv)]; rf <- rf[is.finite(rf)]
  if (length(cv) < 2 || length(rf) < 2) next
  m <- resample_means(cv, rf); m[, `:=`(feature = fc, protein = pr, iter = seq_len(N_ITER))]
  means_list[[paste(fc, pr)]] <- m
  d <- m$cave - m$reference; obs <- mean(d)
  frac <- if (obs > 0) mean(d <= 0) else mean(d >= 0)
  sig_list[[paste(fc, pr)]] <- data.table(feature = fc, protein = pr,
                                          p_lab = fmt_p(frac), ytop = max(c(m$cave, m$reference)))
}
means_dt <- rbindlist(means_list); sig_dt <- rbindlist(sig_list)

longm <- melt(means_dt, id.vars = c("feature","protein","iter"),
              measure.vars = c("cave","reference"), variable.name = "origin", value.name = "mval")
longm[, Origin  := factor(fifelse(origin == "cave","Cave","IMG/VR"), levels = c("Cave","IMG/VR"))]
longm[, protein := factor(protein, levels = protein_order)]
sig_dt[, protein := factor(protein, levels = protein_order)]

# beeswarm display subsample (violins + p use ALL 1000 iterations)
set.seed(1); disp <- longm[, .SD[sample(.N, min(.N, DISPLAY_PTS))], by = .(feature, protein, Origin)]

# --- Style -----------------------------------------------------------
BASE_FAMILY <- "Arial"; BASE_SIZE <- 20
origin_colors <- c(Cave = "#e7298a", "IMG/VR" = "#0072B2")
nice <- function(x) { s <- gsub("_"," ", x); s <- paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s))); sub("^Dna ","DNA ", s) }

make_panel <- function(fc, show_x) {
  d  <- longm[feature == fc]; dd <- disp[feature == fc]; pv <- copy(sig_dt[feature == fc])
  span <- diff(range(d$mval)); pv[, y_pos := ytop + span * 0.06]
  ggplot(d, aes(protein, mval, fill = Origin)) +
    geom_violin(position = position_dodge(0.8), width = 0.8, alpha = 0.35,
                colour = NA, scale = "width", trim = FALSE) +
    geom_quasirandom(data = dd, aes(colour = Origin), dodge.width = 0.8,
                     size = 0.3, alpha = 0.4) +
    geom_text(data = pv, aes(x = protein, y = y_pos, label = p_lab), inherit.aes = FALSE,
              family = BASE_FAMILY, size = BASE_SIZE/.pt, angle = 90, hjust = 0, vjust = 0.5, colour = "black") +
    scale_fill_manual(values = origin_colors, name = NULL) +
    scale_colour_manual(values = origin_colors, name = NULL) +
    scale_x_discrete(labels = nice) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.55))) +
    labs(x = NULL, y = features[[fc]]) +
    theme_bw(base_family = BASE_FAMILY, base_size = BASE_SIZE) +
    theme(text = element_text(family = BASE_FAMILY, size = BASE_SIZE, colour = "black"),
          axis.text.y  = element_text(size = BASE_SIZE, colour = "black"),
          axis.title.y = element_text(size = BASE_SIZE, colour = "black"),
          axis.text.x  = if (show_x) element_text(size = BASE_SIZE, angle = 45, hjust = 1, colour = "black") else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(family = BASE_FAMILY, size = BASE_SIZE, colour = "black"),
          legend.position = "bottom")
}

codes  <- names(features); n <- length(codes)
panels <- lapply(seq_len(n), function(i) make_panel(codes[i], show_x = (i %in% c(n-1, n))))
final_plot <- wrap_plots(panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom",
        plot.tag = element_text(family = BASE_FAMILY, size = 20, face = "plain", colour = "black"))

ggsave(file.path(FIGOUT_S10, "FigS10_physicochemical_features_paired.svg"),
       plot = final_plot, device = svglite::svglite, width = 26, height = 32, limitsize = FALSE)
