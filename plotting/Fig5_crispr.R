install.packages("readr")
library(readr)
install.packages("stringr")
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
prokka_gff <- read_tsv('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/data_from_chrats/PROKKA_04262022.gff', 
                       comment = "#",
                       col_names = FALSE)
# select the useful columns: [1]genome, [4]start, [5]stop, [7]strand, [9]annotation
prokka_annotation <- prokka_gff[,c(1,4,5,7,9)]
colnames(prokka_annotation) <- c('genome', 'start', 'stop', 'strand', 'annotation')
# extract the information of ID and prokka annotation from the annotation column
prokka_annotation$ID <- str_extract(prokka_annotation$annotation, "(?<=ID=)[^;]+")
prokka_annotation$prokka_annotation <- str_extract(prokka_annotation$annotation, "(?<=product=)[^;]+")

# read the ferroplasma eggnog annotation file
eggnog_tsv <- read_tsv('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/data_from_chrats/Ferroplasma.emapper.annotations.tsv', 
                       comment = "##",
                       col_names = TRUE)

# select the useful columns: [1]ID, [8]eggnog_annotation
eggnog_annotation <- eggnog_tsv[,c(1,8)]
colnames(eggnog_annotation) <- c('ID', 'eggnog_annotation')

# combine the prokka and eggnog annotations
prokka_eggnog_annotation <- left_join(prokka_annotation, eggnog_annotation, by = 'ID')

# write.table(prokka_eggnog_annotation, '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/modified_coordinates/prokka_eggnog_annotation.txt', 
#             sep = '\t', 
#             col.names = T,
#             row.names = F,
#             quote = F)
##################################################################################
# plot the whole region
##################################################################################
whole_region_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/modified_coordinates/between_transposase_prokka_eggnog_annotation.txt', sep = '\t', header = 1)

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
CRISPR_target_df <- read_tsv('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/data_from_chrats/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt', 
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
ferroplasma_crispr_f <- '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt'
ferroplasma_crispr_df <- read.table(ferroplasma_crispr_f, sep = '\t', header = 1)
df <- ferroplasma_crispr_df %>%
  separate(Spacer_index, into = c("spacer_order", "spacer_location"), sep = "\\|")
df <- df %>%
  separate(spacer_location, into = c("spacer_start", "spacer_stop"), sep = "\\-")
df$spacer_start <- as.numeric(df$spacer_start)
df$spacer_stop <- as.numeric(df$spacer_stop)
df_simple <- df %>% select(Spacer_ID, spacer_start, spacer_stop, Score, Protospacer_description)
write.table(df_simple, '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_coordinates.txt', sep = '\t', row.names = F, col.names = T, quot = F)

ferroplasma_gff <- read_tsv('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/sulfur_cave_data/complete_genomes/Ferroplasma_complete_genome/Ferroplasma_Prokka/PROKKA_04262022.gff',
                            comment = "#",
                            col_names = FALSE)

ferroplasma_gff_simple <- ferroplasma_gff[10:30,]
ferroplasma_gff_simple2 <- ferroplasma_gff_simple[,c(1,4,5,7,9)]
colnames(ferroplasma_gff_simple2) <- c('genome', 'start', 'stop', 'strand', 'note')
ferroplasma_gff_simple2$description <- str_extract(ferroplasma_gff_simple2$note, "(?<=product=)[^;]+")

write.table(ferroplasma_gff_simple2, '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_region_other_coordinates.txt', sep = '\t', row.names = F, col.names = T, quot = F)


plot_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/sulfur_cave_data/ferroplasma_crispr/ferroplasma_crispr_region_coordinates.txt', sep = '\t', header = 1)

# CRISPR region plot
ggplot(plot_df, aes(xmin = start, xmax = stop, y = genome, fill = description, forward = strand == "+")) +
  geom_gene_arrow() +
  facet_wrap(~ genome, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")

plot_spacer_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/ferroplasma_crispr_coordinates.txt', sep = '\t', header = 1)

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
ferroplasma_crispr_original_output_df <- read.table('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/ferroplasma_CRISPR_annotations/data_from_chrats/Tue_Mar_19_02_20_20_2024_80063_CRISPRTarget_text_report.txt', 
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
