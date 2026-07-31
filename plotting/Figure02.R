# =============================================================================
# Figure 2a: UpSet plot of gene cluster intersections (cave pangenome)
# =============================================================================
# Builds a genome x gene_cluster_id presence/absence matrix from the anvi'o
# pangenome summary table, then plots set/intersection sizes with UpSetR.
# =============================================================================

library(dplyr)
library(UpSetR)
library(ggplotify)
library(grid)
library(cowplot)

# -----------------------------------------------------------------------------
# 1. Load pangenome summary table
# -----------------------------------------------------------------------------
data_dir <- "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics"

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
# f2c unique amino acid KOs per KEGG module,
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
          "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/ferro_non_ferro_aa_module_ko_gap.csv", 
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
# Figure 2c: KEGG/KOfam functional enrichment barplot
# =============================================================================

Ferro_cave_un_enriched.KEGG_Module <- read.delim(file.path(data_dir, "Ferro_cave_un_enriched-KEGG_Module.txt"))
Ferro_cave_un_enriched.KEGG_Module$category <- "KEGG_Module"
colnames(Ferro_cave_un_enriched.KEGG_Module)[1] <- "Function"

# write.table(data_enrichment,
#             "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/data_enrichment.txt",
#             sep = '\t',
#             quote = F,
#             col.names = T,
#             row.names = F
#             )

data_enrichment <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/data_enrichment.txt")

## ---- classify each function by which set carries it ----------------------
# associated_groups lists the genome groups a function is enriched in.
# The two cave lineages are "Ferroplasma_c" and "MAG00006"; "rest" = reference.
library(dplyr)
library(stringr)
library(ggplot2)

classify_side <- function(g, sig) {
  in_fc   <- str_detect(g, "Ferroplasma_c")
  in_mag  <- str_detect(g, "MAG00006")
  in_ref  <- str_detect(g, "rest")
  in_cave <- in_fc | in_mag
  
  case_when(
    is.na(g)                           ~ "Shared: reference + both cave",
    in_cave & in_ref & in_fc & !in_mag ~ "Shared: reference + Ferroplasma_c",
    in_cave & in_ref & in_mag & !in_fc ~ "Shared: reference + MAG00006",
    in_cave & in_ref                   ~ "Shared: reference + both cave",
    in_cave & sig == "Yes"             ~ "Significant cave-only",
    in_cave & sig == "No"              ~ "Cave-only",
    in_ref                             ~ "Reference-only",
    TRUE                               ~ NA_character_
  )
}
# Note: Using your streamlined version here to avoid rowwise vapply issues
df_enrichment <- data_enrichment %>%
  mutate(side = classify_side(associated_groups, Significance))

## ---- counts + proportions per layer --------------------------------------
df_enrichment_counts <- df_enrichment %>%
  count(category, side, name = "n") %>%
  group_by(category) %>%
  mutate(prop = n / sum(n),
         total = sum(n)) %>%
  ungroup()

# readable layer labels carrying the per-layer n
layer_labels <- df_enrichment_counts %>%
  distinct(category, total) %>%
  mutate(label = dplyr::recode(category,
                               KOfam       = paste0("KOs \n(n=", total, ")"),
                               KEGG_Module = paste0("KEGG modules\n(n=", total, ")")))

# stack order (bottom -> top): Significant cave-only, Cave-only, Shared..., Reference-only
df_enrichment_counts <- df_enrichment_counts %>%
  left_join(layer_labels %>% select(category, label), by = "category") %>%
  mutate(
    label = factor(label, levels = layer_labels$label[order(-layer_labels$total)]),
    side  = factor(side, levels = c("Significant cave-only",
                                    "Cave-only", 
                                    "Shared: reference + Ferroplasma_c",
                                    "Shared: reference + MAG00006",
                                    "Shared: reference + both cave",
                                    "Reference-only")),
    
    # ── FIX 1: Only keep labels for Cave-only and Reference-only ──────────────
    seg_label = ifelse(
      side %in% c("Cave-only", "Reference-only"),
      paste0(round(prop * 100), "%\n(n=", n, ")"),
      "" # Leaves the shared segments blank
    )
  )

# Compute label y-positions explicitly
df_enrichment_counts <- df_enrichment_counts %>%
  arrange(category, desc(side)) %>%      
  group_by(category) %>%
  mutate(ymax = cumsum(prop),
         ymid = ymax - prop / 2) %>%
  ungroup()

## ---- palette -------------------------------------------------------------
pal <- c("Significant cave-only" = "#dd1c77", 
         "Cave-only" = "#ffd92f", 
         "Shared: reference + Ferroplasma_c" = "#66c2a5",
         "Shared: reference + MAG00006" = "#fc8d62",
         "Shared: reference + both cave" = "#b8bfc9", 
         "Reference-only" = "#7570b3"
)

# ── FIX 2: Correct string matching for text coloring ────────────────────────
df_enrichment_counts <- df_enrichment_counts %>% 
  mutate(txt_col = ifelse(str_detect(side, "Shared"), "#444444", "white"))

## ---- plot -----------------------------------------------------------------
library(ggrepel)

# ── DATA PREP (assuming context from previous steps) ─────────────────────────

# 1. Update seg_label to include "Significant cave-only"
# We update the upstream mutate call where seg_label was created:
# The user already added "Significant cave-only" in previous interactions.
# We modify the logical check so labels are generated for these three groups.

df_enrichment_counts <- df_enrichment_counts %>%
  mutate(
    # Updated logical check to include the new category
    seg_label = ifelse(
      side %in% c("Cave-only", "Reference-only", "Significant cave-only"),
      paste0(round(prop * 100), "%\n(n=", n, ")"),
      "" # Leaves other segments blank
    ),
    
    # Text coloring rule (Shared is light, others are dark)
    txt_col = ifelse(str_detect(side, "Shared"), "#444444", "white")
  )

# 2. Re-compute explicit positions (ensure sorting matches factor levels)
# Stack order (bottom -> top): Significant -> Cave-only -> Shared... -> Reference-only
df_enrichment_counts <- df_enrichment_counts %>%
  left_join(layer_labels %>% select(category, label), by = "category") %>%
  mutate(
    label = factor(label, levels = layer_labels$label[order(-layer_labels$total)]),
    side  = factor(side, levels = c(
      "Significant cave-only",
      "Cave-only", 
      "Shared: reference + Ferroplasma_c",
      "Shared: reference + MAG00006",
      "Shared: reference + both cave",
      "Reference-only"
    ))
  ) %>%
  arrange(category, desc(side)) %>% # forces cumulative sum logic to match factor stack
  group_by(category) %>%
  mutate(ymax = cumsum(prop),
         ymid = ymax - prop / 2) %>%
  ungroup()

# ── DEFINE PALETTE ────────────────────────────────────────────────────────────
pal <- c(
  "Significant cave-only" = "#e31a1c", # Distinct RED for significance
  "Cave-only" = "#ffd92f", 
  "Shared: reference + Ferroplasma_c" = "#66c2a5",
  "Shared: reference + MAG00006" = "#fc8d62",
  "Shared: reference + both cave" = "#b8bfc9", 
  "Reference-only" = "#2c7fb8"
)

# ── THE NEW PLOT IMPLEMENTATION ──────────────────────────────────────────────
f2c <- ggplot(df_enrichment_counts, aes(x = label, y = prop, fill = side)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.6) +
  
  # ── Layer 1: Internal Labels (Standard Cave and Reference-Only) ──────────────
  geom_text(
    data = subset(df_enrichment_counts, side %in% c("Cave-only", "Reference-only")),
    aes(y = ymid, label = seg_label, color = txt_col),
    size = 5, fontface = "bold", lineheight = 0.9,
    show.legend = FALSE
  ) +
  
  # ── Layer 2: External Linking Arrows (Significant Cave-Only) ────────────────
  # Draws a physical link line from the thin segment out into clear canvas space
  geom_segment(
    data = subset(df_enrichment_counts, side == "Significant cave-only" & seg_label != ""),
    aes(
      x = as.numeric(label),        # Start line right at the bar's x-coordinate
      xend = as.numeric(label) + 0.43, # Extend the line out to the right side
      y = ymid,                     # Track exactly centered within the segment
      yend = ymid
    ),
    color = "black",
    linewidth = 0.5,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed", ends = "first"), # Arrow heads pointing at the slice
    inherit.aes = FALSE            # Prevents aesthetic inheritance conflicts
  ) +
  
  # ── Layer 3: Fixed External Text Labels (Significant Cave-Only) ──────────────
  # Places the label safely directly next to the end tip of our custom linking line
  geom_text(
    data = subset(df_enrichment_counts, side == "Significant cave-only" & seg_label != ""),
    aes(
      x = as.numeric(label) + 0.45, # Positions text just past the arrow pointer line termination
      y = ymid, 
      label = seg_label
    ),
    color = "black",                # Solid crisp black text for out-of-bar canvas regions
    size = 4.5, 
    fontface = "bold", 
    lineheight = 0.9,
    hjust = 0,                      # Left-align text box start for neat alignments
    inherit.aes = FALSE
  ) +
  
  # ── Remainder of your standard plotting layers remain identical ──────────────
  scale_fill_manual(
    values = pal, 
    breaks = c(
      "Significant cave-only",
      "Cave-only", 
      "Shared: reference + Ferroplasma_c", 
      "Shared: reference + MAG00006",
      "Shared: reference + both cave",
      "Reference-only"
    ),
    labels = list(
      expression("Significant cave-only (" * italic("Ferroplasma") ~ "c. and/or MAG 6)"),
      expression("Cave-only (" * italic("Ferroplasma") ~ "c. and/or MAG 6)"),
      expression("Shared: reference + " * italic("Ferroplasma") ~ "c."),
      expression("Shared: reference + MAG 6"),
      expression("Shared: reference + both cave "* italic("Ferroplasma")),
      expression("Reference-only")
    ),
    name = "Functional enrichment distribution"
  ) +
  scale_colour_identity() +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = NULL, y = "Share of enriched functions"
  ) +
  theme_minimal(base_size = 20) + 
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    
    axis.ticks.x       = element_blank(),
    axis.line          = element_line(color = "black"), 
    
    legend.position    = "right", 
    
    text               = element_text(color = "black"),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(color = "black"),
    legend.text        = element_text(color = "black"),
    legend.title       = element_text(color = "black"),
    
    plot.margin        = margin(10, 50, 10, 10) # Expanded right margin margin to prevent labels cutting off
  )

f2c
