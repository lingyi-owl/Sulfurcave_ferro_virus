# =============================================================================
# Fig. 5 — virus–microbe interactions
# Figure numbering follows the manuscript (Cave_Ferro_Virus, Nov 2026).
# Part 1 (from Fig5_virome.R): a = virus class / host genus alluvial, MAGs included.
# Part 2 (from Fig5_crispr.R): b = Ferroplasma c. MAG 8 CRISPR array, d = spacer scores.
# Panel c (Cas protein LFQ) is built in Fig4_proteome_pangenome_defense.R
#   (object ferro_crispr_lfq_p), because it reuses that script's proteomics tables.
# =============================================================================

# ############################ Part 1 ############################
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
EXTERNAL <- Sys.getenv("SULFURCAVE_EXTERNAL", unset = file.path(REPO, "external_data"))  # inputs not shipped in the repo

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)

# --- 1. DATA LOADING & PREPARATION ---
no_mags_file <- file.path(EXTERNAL, "vcat_iphop_abundance_202602/vcat_iphop_abundance_no_mags.txt")
with_mags_file <- file.path(EXTERNAL, "vcat_iphop_abundance_202602/vcat_iphop_abundance_with_mags.txt")
  
df_no_mags <- read.csv2(no_mags_file, sep = '\t')
df_with_mags <- read.csv2(with_mags_file, sep = '\t')

df <- df_with_mags   # Fig. 5a: MAGs included (Fig. S8 uses df_no_mags)

df <- df %>%
  mutate(across(starts_with("Sample"), ~ as.numeric(as.character(.x))))

# Convert from wide to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("Sample"), 
    names_to = "sample", 
    values_to = "abundance"
  ) %>%
  mutate(sample = gsub("_abundance", "", sample)) 

# --- 2. TRANSFORMATION TO LODES FORMAT ---
lodes_all <- df_long %>%
  transmute(contig_id, Virus = vcat_class, Host = iphop_genus, weight = abundance, sample) %>%
  ggalluvial::to_lodes_form(
    axes  = c("Virus", "Host"),
    key   = "axis",
    value = "stratum",
    id    = "contig_id"
  ) %>%
  mutate(
    stratum = as.character(stratum),
    # Handle missing values as "Unclassified"
    stratum = ifelse(is.na(stratum) | stratum == "", "Unclassified", stratum),
    stratum = factor(stratum)
  )

# --- 3. AESTHETICS (MAPS & COLORS) ---
name_map <- c(
  "Sample_ERR10036468" = "Biofilm 1",
  "Sample_ERR10036469" = "Biofilm 2",
  "Sample_ERR10036470" = "Lab CH4"
)

my_colors <- c(
  # --- Virus Taxa ---
  "Caudoviricetes" = "#e7298a",
  "Faserviricetes" = "#2ca02c",
  "Megaviricetes"  = "#d62728",
  "Naldaviricetes" = "#1f77b4",
  "Below 1%"       = "#8da0cb",   # virus block
  "Unclassified"   = "gray",      # virus block
  
  # --- Host Genera ---
  "Acidithiobacillus"     = "#fc8d62",   # LOCKED (group 2)
  "Below 1%"              = "#8da0cb",    # LOCKED (group 2)
  "Cuniculiplasma"        = "#e5c494",    # LOCKED (group 2)
  "Ferroplasma"           = "#ffd92f",    # LOCKED (group 2)
  "Mycobacterium"         = "#66c2a5",    # LOCKED (group 2)
  "Unclassified"          = "#b3b3b3",    # LOCKED (group 2)
  
  "Escherichia"           = "#C49C94FF",  # kept (was fine)
  "JAKAFX01"              = "#9EDAE5FF",   # kept (was fine)
  "Sulfobacillus"         = "#98DF8AFF",   # kept (was fine)
  
  "Igneacidithiobacillus" = "#4d004b",    # CHANGED (was #decbe4, clashed w/ Marinobacter)
  "Marinobacter"          = "#005a32",    # CHANGED (was #fddaec, clashed w/ Igneacidithiobacillus)
  "Lachnospira"           = "#525252",    # CHANGED (was #e5d8bd, clashed w/ Cuniculiplasma/Streptomyces/Mobilitalea)
  "Mobilitalea"           = "#02818a",    # CHANGED (was #DBDB8DFF, clashed w/ Streptomyces/Lachnospira)
  "Streptomyces"          = "#993404"     # CHANGED (was #ffffcc, clashed w/ Mobilitalea/Lachnospira)
)

# --- 4. PLOTTING ---
p <- ggplot(lodes_all,
            aes(x = axis, stratum = stratum, alluvium = contig_id, 
                y = weight * 100, fill = stratum)) +
  # Flows
  geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
  # Vertical bars
  geom_stratum(width = 0.3, color = "grey30") +
  # Scales & Colors
  scale_x_discrete(limits = c("Virus", "Host"), expand = c(.12, .12)) +
  scale_fill_manual(values = my_colors, name = "Virus / Host") +
  # Multi-panel display by Sample
  facet_wrap(~ sample, nrow = 1, labeller = labeller(sample = name_map)) +
  # Labels
  labs(
    x = "Virus class --> Host genus", y = "% mapped reads",
    # title = "Virus and predicted host abundance",
    fill = "Taxa / Hosts"
  ) +
  # Theme Customization
  theme_minimal(base_family = "Arial", base_size = 20) +
  theme(
    legend.position = "none",
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    panel.grid.major.x = element_blank() # Clean up vertical lines
  )

p

# ggsave(
#   filename = file.path(REPO, "figures", "no_mags_color_legend.svg"), 
#   plot = p, 
#   device = "svg", 
#   width = 16, 
#   height = 10, 
#   units = "in",
#   limitsize = FALSE
# )
# legend
library(ggplot2)
library(cowplot)

virus_colors <- c(
  "Caudoviricetes" = "#e7298a",
  "Faserviricetes" = "#2ca02c",
  "Megaviricetes"  = "#d62728",
  "Naldaviricetes" = "#1f77b4",
  "Below 1%"       = "#8da0cb",
  "Unclassified"   = "#b3b3b3"
)
host_colors <- c(
  "Acidithiobacillus"     = "#fc8d62",
  "Below 1%"              = "#8da0cb",
  "Cuniculiplasma"        = "#e5c494",
  "Escherichia"           = "#C49C94FF",
  "Ferroplasma"           = "#ffd92f",
  "Igneacidithiobacillus" = "#4d004b",
  "JAKAFX01"              = "#9EDAE5FF",
  "Lachnospira"           = "#525252",
  "Marinobacter"          = "#005a32",
  "Mobilitalea"           = "#02818a",
  "Mycobacterium"         = "#66c2a5",
  "Streptomyces"          = "#993404",
  "Sulfobacillus"         = "#98DF8AFF",
  "Unclassified"          = "#b3b3b3"
)

# entries that are NOT taxon names and should stay upright
non_taxa <- c("Below 1%", "Unclassified")

# build italic expression labels for taxon names, plain for the rest
make_labels <- function(cols, italicize) {
  nm <- names(cols)
  if (!italicize) return(nm)
  setNames(
    lapply(nm, function(x) if (x %in% non_taxa) x else bquote(italic(.(x)))),
    nm
  )
}

legend_only <- function(cols, title, italicize = FALSE) {
  d <- data.frame(cat = factor(names(cols), levels = names(cols)), x = 1, y = 1)
  ggplot(d, aes(x, y, fill = cat)) +
    geom_col() +
    scale_fill_manual(values = cols, name = title,
                      labels = make_labels(cols, italicize),
                      guide = guide_legend(ncol = 1, byrow = TRUE)) +
    theme_void(base_family = "Arial", base_size = 20) +
    theme(
      legend.justification = c(0, 1),
      legend.title = element_text(face = "bold", hjust = 0),
      legend.text  = element_text(),
      legend.key.size = unit(0.65, "cm"),
      legend.key.spacing.y = unit(0.35, "cm")
    )
}

leg_v <- get_legend(legend_only(virus_colors, "Virus class"))                 # virus stays upright
leg_h <- get_legend(legend_only(host_colors,  "Host genus", italicize = TRUE)) # host genera italic
combined <- plot_grid(leg_v, leg_h, ncol = 2, rel_widths = c(0.9, 1))
combined
# ggsave("legends_virus_host.png", combined, width = 7.5, height = 5, dpi = 300, bg = "white")





# ############################ Part 2 ############################
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
EXTERNAL <- Sys.getenv("SULFURCAVE_EXTERNAL", unset = file.path(REPO, "external_data"))  # inputs not shipped in the repo

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(gggenes)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# read the ferroplasma prokka annotation file
prokka_gff <- read_tsv(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/data_from_chrats/PROKKA_04262022.gff"), 
                       comment = "#",
                       col_names = FALSE)
# select the useful columns: [1]genome, [4]start, [5]stop, [7]strand, [9]annotation
prokka_annotation <- prokka_gff[,c(1,4,5,7,9)]
colnames(prokka_annotation) <- c('genome', 'start', 'stop', 'strand', 'annotation')
# extract the information of ID and prokka annotation from the annotation column
prokka_annotation$ID <- str_extract(prokka_annotation$annotation, "(?<=ID=)[^;]+")
prokka_annotation$prokka_annotation <- str_extract(prokka_annotation$annotation, "(?<=product=)[^;]+")

# read the ferroplasma eggnog annotation file
eggnog_tsv <- read_tsv(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/data_from_chrats/Ferroplasma.emapper.annotations.tsv"), 
                       comment = "##",
                       col_names = TRUE)

# select the useful columns: [1]ID, [8]eggnog_annotation
eggnog_annotation <- eggnog_tsv[,c(1,8)]
colnames(eggnog_annotation) <- c('ID', 'eggnog_annotation')

# combine the prokka and eggnog annotations
prokka_eggnog_annotation <- left_join(prokka_annotation, eggnog_annotation, by = 'ID')

# write.table(prokka_eggnog_annotation, file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/modified_coordinates/prokka_eggnog_annotation.txt"), 
#             sep = '\t', 
#             col.names = T,
#             row.names = F,
#             quote = F)
##################################################################################
# plot the whole region
##################################################################################
whole_region_df <- read.table(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/modified_coordinates/between_transposase_prokka_eggnog_annotation.txt"), sep = '\t', header = 1)

ferro_crispr_array <- ggplot(whole_region_df, aes(xmin = start, xmax = stop, y = genome, fill = Group, forward = strand == "+")) +
  geom_gene_arrow() +
  scale_fill_manual(values = c('#FDB462', '#FB8072', '#D9D9D9', '#FCCDE5', '#B3DE69')) +
  ggtitle("Ferroplasma Circular - CAS-TypeID+B") +
  labs(y="")+
  geom_text_repel(aes(x = (start+stop)/2, y = 1, label = molecule, angle=60), nudge_x = 0.15,
                  box.padding = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 1,
                  segment.angle = 10) +
  theme(legend.position="bottom") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    # panel.grid = element_blank(),       # removes both major and minor grid lines
    # panel.grid.major = element_blank(), # optional, explicitly removes major grid lines
    panel.grid.minor = element_blank()  # optional, explicitly removes minor grid lines
  )

ferro_crispr_array
##################################################################################
# plot the CRISPR spacers region with spacer origin
##################################################################################
# read the ferroplasma spacer annotation file
CRISPR_target_df <- read_tsv(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/data_from_chrats/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt"), 
                       comment = "#",
                       col_names = T)
crispr_df <- CRISPR_target_df %>% 
  select('Spacer_ID', 'Spacer_index', 'Score', 'Protospacer_description')
crispr_df <- crispr_df %>%
  separate(Spacer_index, into = c("spacer_order", "spacer_location"), sep = "\\|")
crispr_df <- crispr_df %>%
  separate(spacer_location, into = c("from", "to"), sep = "\\-")
crispr_df$from <- as.numeric(crispr_df$from)
crispr_df$to <- as.numeric(crispr_df$to)
crispr_df <- crispr_df %>%
  mutate(origin = ifelse(startsWith(Protospacer_description, "c_"), "cave virus", "viruses from \n CRISPRTarget \n databases")) 

names(crispr_df)[names(crispr_df) == 'Spacer_ID'] <- 'genome'
names(crispr_df)[names(crispr_df) == 'Score'] <- 'score'
names(crispr_df)[names(crispr_df) == 'Protospacer_description'] <- 'subgene'

crispr_df$gene <- 'CRISPR_repeats'
crispr_df$start <- 13622
crispr_df$stop <- 21729
crispr_df$strand <- '.'
crispr_df <- crispr_df[!is.na(crispr_df$origin),]

ggplot(crispr_df, aes(xmin = start, xmax = stop, y = genome, forward = strand == "+")) +
  geom_gene_arrow() +
  scale_fill_manual(
    values = c(
      "cave virus"                              = "#ffd92f",   # Cave virus
      "viruses from \n CRISPRTarget \n databases" = "#2c7fb8"  # Reference virus
    )
  ) +
  geom_subgene_arrow(
    data = crispr_df,
    aes(xsubmin = from, xsubmax = to, fill = origin),
    color = NA
  ) +
  theme_void() +
  theme(
    legend.position = "none"
  )

##################################################################################
# plot the CRISPR spacers scores
##################################################################################
color_df <- data.frame(
  score = 20:36,
  color = c("#FFC0CBFF", "#FFC0CBFF", "#FFC0CBFF", "#FFC0CBFF", "#FFC0CBFF", "#FFC0CBFF", "#FFC0CBFF",
            "#F5ABA2FF", "#E9967AFF", "#F86643FF", "#FF0000FF", "#FF0000FF", "#FF0000FF", "#FF0000FF",
            "#FF0000FF", "#FF0000FF", "#FF0000FF"),
  stringsAsFactors = FALSE
)
crispr_df <- left_join(crispr_df, color_df, by = 'score')

# CRISPR spacer score plot
ggplot(crispr_df, aes(xmin = start, xmax = stop, y = genome, forward = strand == "+")) +
  geom_gene_arrow() +
  geom_subgene_arrow(
    data = crispr_df,
    aes(xsubmin = from, xsubmax = to, fill = color)
  ) +
  scale_fill_identity() +  
  theme_void()

ggplot(crispr_df, aes(xmin = start, xmax = stop, y = genome, forward = strand == "+")) +
  geom_gene_arrow() +
  geom_subgene_arrow(
    data = crispr_df,
    aes(xsubmin = from, xsubmax = to, fill = color),
    color = NA  # remove the border
  ) +
  scale_fill_identity() +  
  theme_void()

# Extract the legend from the dummy plot
legend <- get_legend(dummy_plot)

# Display the legend
plot_grid(NULL, legend, ncol = 1, rel_heights = c(0.1, 0.9))

ggplot(crispr_df, aes(xmin = start, xmax = stop, y = genome, fill = color)) +
  geom_gene_arrow() +
  scale_fill_identity() +  # Use the colors specified in the data
  facet_wrap(~ genome, scales = "free", ncol = 1) 


##################################################################################
ferroplasma_crispr_f <- file.path(EXTERNAL, "ferroplasma_CRISPR/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt")
ferroplasma_crispr_df <- read.table(ferroplasma_crispr_f, sep = '\t', header = 1)
df <- ferroplasma_crispr_df %>%
  separate(Spacer_index, into = c("spacer_order", "spacer_location"), sep = "\\|")
df <- df %>%
  separate(spacer_location, into = c("spacer_start", "spacer_stop"), sep = "\\-")
df$spacer_start <- as.numeric(df$spacer_start)
df$spacer_stop <- as.numeric(df$spacer_stop)
df_simple <- df %>% select(Spacer_ID, spacer_start, spacer_stop, Score, Protospacer_description)
write.table(df_simple, file.path(EXTERNAL, "sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_coordinates.txt"), sep = '\t', row.names = F, col.names = T, quot = F)

ferroplasma_gff <- read_tsv(file.path(EXTERNAL, "sulfur_cave_data/complete_genomes/Ferroplasma_complete_genome/Ferroplasma_Prokka/PROKKA_04262022.gff"),
                            comment = "#",
                            col_names = FALSE)

ferroplasma_gff_simple <- ferroplasma_gff[10:30,]
ferroplasma_gff_simple2 <- ferroplasma_gff_simple[,c(1,4,5,7,9)]
colnames(ferroplasma_gff_simple2) <- c('genome', 'start', 'stop', 'strand', 'note')
ferroplasma_gff_simple2$description <- str_extract(ferroplasma_gff_simple2$note, "(?<=product=)[^;]+")

write.table(ferroplasma_gff_simple2, file.path(EXTERNAL, "sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_region_other_coordinates.txt"), sep = '\t', row.names = F, col.names = T, quot = F)


plot_df <- read.table(file.path(EXTERNAL, "sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_region_coordinates.txt"), sep = '\t', header = 1)

# CRISPR region plot
ggplot(plot_df, aes(xmin = start, xmax = stop, y = genome, fill = description, forward = strand == "+")) +
  geom_gene_arrow() +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")

plot_spacer_df <- read.table(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/ferroplasma_crispr_coordinates.txt"), sep = '\t', header = 1)

color_df <- data.frame(
  score = 20:36,
  color = c("#FFB6C1", "#FF69B4", "#FF1493", "#FF007F", "#C71585",
                "#FFA07A", "#F08080", "#FA8072", "#E9967A", "#FF7F50",
                "#FF6347", "#FF0000", "#CB4154", "#CD5C5C", "#DC143C",
                "#B22222", "#8B0000"),
  stringsAsFactors = FALSE
)
plot_spacer_color_df <- left_join(plot_spacer_df, color_df, by = 'score')

# CRISPR spacer score plot
ggplot(plot_spacer_color_df, aes(xmin = start, xmax = stop, y = genome, fill = color)) +
  geom_gene_arrow() +
  scale_fill_identity() +  # Use the colors specified in the data
  facet_wrap(~ genome, scales = "free", ncol = 1) 

# CRISPR spacer origin plot
plot_spacer_df <- plot_spacer_df %>%
  mutate(origin = if_else(grepl("c_", Protospacer_description), "cave viruses", "database"))

spacer_origin_plot <- ggplot(plot_spacer_df, aes(xmin = start, xmax = stop, y = genome, fill = origin)) +
  geom_gene_arrow() +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal()
spacer_origin_plot + theme(legend.position="right")

##############################################################
ferroplasma_crispr_original_output_df <- read.table(file.path(EXTERNAL, "ferroplasma_CRISPR_annotations/data_from_chrats/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt"), 
                                                    comment.char = "#",
                                                    sep = "\t",             # or "" if it's whitespace-separated
                                                    header = FALSE,         # or TRUE if there’s a header line after the comments
                                                    fill = TRUE,            # Fills in missing columns
                                                    quote = ""              # Avoid issues with quote characters
)

#####################################################################
# legend
#####################################################################
library(ComplexHeatmap)
library(circlize)

# Define color function
col_fun <- colorRamp2(
  c(20, 36),
  c("#FFC0CBFF", "#FF0000FF")
)

# Create horizontal color legend with black outside ticks
lgd <- Legend(
  title = NULL,
  col_fun = col_fun,
  direction = "horizontal",
  at = c(20, 25, 30, 35),
  legend_width = unit(4, "cm"),
  labels_gp = gpar(fontsize = 9, col = "black"),
  tick_length = unit(0.2, "cm"),
  tick_gp = gpar(col = "black"),
  tick_side = "top",  # for ticks *outside* the bar
  grid_height = unit(0.4, "cm")
)

# Draw the legend
draw(lgd)

#############################################################
# box plots for CRISPR spacer origin and matching scores
#############################################################
crispr_df_unique <- crispr_df %>%
  distinct(spacer_order, from, to, score, origin, .keep_all = TRUE)

crispr_df_unique$Origin <- ifelse(crispr_df_unique$origin == 'cave virus', 'Sulfur cave phage', 'Database phage')
crispr_df_unique$Origin <- factor(crispr_df_unique$Origin, levels = c('Sulfur cave phage', 'Database phage'))

violin_color <- c('Sulfur cave phage' = '#ffd92f', 'Database phage' = '#2c7fb8')
point_color <- c('Sulfur cave phage' = '#ffd92f', 'Database phage' = '#2c7fb8')

library(ggsignif)

# compute and format the p-value
pv  <- wilcox.test(score ~ Origin, data = crispr_df_unique)$p.value
lab <- paste0("p-value = ", format(pv, digits = 2))   # e.g. "p-value = 4.2e-10"

ggplot(crispr_df_unique, aes(x = Origin, y = score)) +
  geom_violin(aes(fill = Origin), trim = FALSE, color = NA, alpha = 0.6) +
  geom_beeswarm(aes(color = Origin), size = 1.5) +
  scale_fill_manual(
    values = violin_color,
    labels = c("Database phage" = "Reference", "Sulfur cave phage" = "Cave"),
    name = "Virus origin"
  ) +
  scale_color_manual(
    values = point_color,
    labels = c("Database phage" = "Reference", "Sulfur cave phage" = "Cave"),
    name = "Virus origin"
  ) +
  scale_x_discrete(labels = c(
    "Database phage"    = "Reference",
    "Sulfur cave phage" = "Cave"
  )) +
  geom_signif(
    comparisons = list(c("Database phage", "Sulfur cave phage")),
    annotations = lab,          # custom label with the prefix
    textsize    = 7,
    size        = 0.6,          # bracket line thickness (was bracket.size)
    tip_length  = 0.02          # note underscore, not dot
  ) +
  theme_minimal(base_size = 20, base_family = "Arial") +
  labs(
    x = "Virus origin",
    y = "CRISPR spacer matching score"
  ) +
  theme(
    axis.line    = element_line(color = "black", linewidth = 0.5),
    axis.ticks   = element_line(color = "black", linewidth = 0.5),
    panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.8),
    legend.key   = element_rect(fill = NA, color = NA)
  )
