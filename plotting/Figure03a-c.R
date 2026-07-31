setwd('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/')
################################################
# fig3a
################################################
################ combine proteomics
# metaproteomics were mapped against 92103 proteins from MAGs and virous 
ProteinGroups <- readxl::read_xlsx("proteomics/proteinGroups.xlsx", sheet = 2)
ids <- readxl::read_xlsx("proteomics/proteinGroups.xlsx", sheet = 1)

data_proteomics <- as.data.frame(ProteinGroups[,1:2])
data_proteomics$LFQ_intensity_Cave_Biofilm_1 <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_1`)
data_proteomics$LFQ_intensity_Cave_Biofilm_2 <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_2`)
data_proteomics$LFQ_intensity_Culture <- as.numeric(ProteinGroups$`LFQ intensity invitro`)
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_1 <- as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_1`))
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 <- as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_2`))
data_proteomics$log10_LFQ_intensity_Culture <- as.numeric(log10(ProteinGroups$`LFQ intensity invitro`))
data_proteomics[data_proteomics == "-Inf"] <- 0

# taxonomy join
colors_pB <- c('#a6d854', "#fb8072", '#8da0cb', '#66c2a5', "grey")
info_table_DNA <- read.csv("metagenomics/info_table_DNA.csv", sep = ';')
info_table_DNA$Alternative_ID2 <- gsub("_", "", info_table_DNA$Alternative_ID2)
data_proteomics$ID <- unlist(lapply(strsplit(data_proteomics$`Protein IDs`, "_"), function(x) x[1]))
data_proteomics_tax <- merge(data_proteomics, info_table_DNA, by.x = "ID", by.y = "Alternative_ID2")
data_proteomics_tax$Genus1 <- data_proteomics_tax$Genus
data_proteomics_tax$Genus1[!data_proteomics_tax$Genus %in% names(which(table(data_proteomics_tax$Genus) > 100))] <- "Others"
data_proteomics_tax$Genus1 <- as.factor(data_proteomics_tax$Genus1)

######################################
# check the data points on each plane
######################################
x <- data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_1
y <- data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2
z <- data_proteomics_tax$log10_LFQ_intensity_Culture

sum(z == 0, na.rm = TRUE)   # x–y plane (Culture == 0)
sum(x == 0, na.rm = TRUE)   # y–z plane (Biofilm 1 == 0)
sum(y == 0, na.rm = TRUE)   # x–z plane (Biofilm 2 == 0)
sum(x != 0 & y != 0 & z != 0, na.rm = TRUE) # present in all three samples
sum(x != 0 & y != 0, na.rm = TRUE) # present in both biofilm samples

################ Panel a: static 3D scatterplot
library(scatterplot3d)

# font/style settings translated to base-R graphical parameters
# (ggplot2 theme_classic/theme_minimal calls do not apply to scatterplot3d output)
par(family = "Arial", cex.lab = 20/12, cex.axis = 20/12, font.lab = 1, font.axis = 1)

# marker shape: distinguishes proteins identified in both cave biofilm samples
data_proteomics_tax$cave_only <- rep(18, nrow(data_proteomics_tax))
data_proteomics_tax$cave_only[data_proteomics_tax$LFQ_intensity_Cave_Biofilm_1 > 0 & data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2 > 0] <- 16

colors <- colors_pB[as.numeric(data_proteomics_tax$Genus1)]
size <- c(2, 1)[factor(data_proteomics_tax$cave_only)]
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')

width  <- 944 / 300  # 3.15 inches
height <- 657 / 300  # 2.19 inches

s3d <- scatterplot3d(x = data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_1,
                     y = data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2,
                     z = data_proteomics_tax$log10_LFQ_intensity_Culture,
                     main = "", pch = data_proteomics_tax$cave_only, color = colors, box = FALSE,
                     cex.symbols = size,
                     cex.axis = 20/12, cex.lab = 20/12, cex.sub = 20/12,
                     col.axis = "black", col.lab = "black",
                     col.grid = "grey92",
                     xlab = "Biofilm 1", ylab = "Biofilm 2", zlab = "Lab CH4",
                     sub = "LFQ intensity (log10)")

addgrids3d(x = data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_1,
           y = data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2,
           z = data_proteomics_tax$log10_LFQ_intensity_Culture,
           grid = c("xy", "xz", "yz"))

dev.copy2pdf(file = "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/figures/Fig3/fig3a.pdf")

# Capture the just-drawn base plot as a grob
grid.echo()  # re-records the last base plot using grid graphics
fig3a_grob <- grid.grab()
fig3a_grob
################################################
# fig3a legend
################################################
library(ggplot2)
library(cowplot)

genus_colors <- c(
  "Ferroplasma"       = '#ffd92f',
  "Mycobacterium"     = '#66c2a5',
  "Acidithiobacillus" = '#fc8d62',
  "Cuniculiplasma"    = '#e5c494',
  "Others"            = "grey"
)

legend_plot <- ggplot(data_proteomics_tax,
                      aes(x = log10_LFQ_intensity_Cave_Biofilm_1,
                          y = log10_LFQ_intensity_Cave_Biofilm_2,
                          color = Genus1,
                          shape = factor(cave_only))) +
  geom_point(alpha = 0) +
  scale_color_manual(values = genus_colors, name = "Genus") +
  scale_shape_manual(values = c("16" = 21, "18" = 23),
                     labels = c("Both cave biofilms", "Not in both"),
                     name = "Origin") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 6)),
         shape = guide_legend(override.aes = list(alpha = 1, size = 6,
                                                  colour = "grey", fill = "white"))) +
  theme_void(base_size = 20) +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"))

legend_only <- cowplot::get_legend(legend_plot)
grid::grid.newpage()

################################################
# fig3b
################################################
# Identify proteins found in both cave biofilm replicates
data_proteomics1 <- data_proteomics[which(data_proteomics$LFQ_intensity_Cave_Biofilm_1 > 0 &
                                            data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 > 0), ]
dim(data_proteomics1)
# 2054

# Extract bin/MAG ID prefix from Protein IDs
data_proteomics1$ID <- unlist(lapply(strsplit(data_proteomics1$`Protein IDs`, "_"), function(x) x[1]))

# Count proteins per MAG/bin
sum_prot <- as.data.frame(table(data_proteomics1$ID))

# Join with taxonomy info to get Genus per MAG
sum_prot_tax <- merge(sum_prot, info_table_DNA, by.x = "Var1", by.y = "Alternative_ID2")
# sum_prot_tax <- sum_prot_tax[-2, ]

# remove the third row
sum_prot_tax <- sum_prot_tax[-3, ]
# # Rename specific MAG IDs to readable labels
# sum_prot_tax$Var <- as.character(sum_prot_tax$Var)
# sum_prot_tax$Var[1:5] <- c("M. MAG2", "M. MAG3", "Cuniculiplasma c.", "M. Methanotrophicum", "Ferroplasma c.")
# sum_prot_tax$Var <- gsub("SC", "", sum_prot_tax$Var)
sum_prot_tax$Alternative_ID3 <- factor(sum_prot_tax$Alternative_ID3, levels = as.character(sum_prot_tax$Alternative_ID3[order(sum_prot_tax$Freq, decreasing = TRUE)]))

# Recode rare genera as "Others"
sum_prot_tax$Genus1 <- sum_prot_tax$Genus
sum_prot_tax$Genus1[which(sum_prot_tax$Freq < 49)] <- "Others"
sum_prot_tax$Genus1 <- factor(sum_prot_tax$Genus1, levels = c("Ferroplasma", "Mycobacterium", "Acidithiobacillus", "Others"))

fig3b <- ggplot(sum_prot_tax, aes(x = Alternative_ID3, y = Freq, fill = Genus1)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Freq), vjust = -0.5, size = 4) +
  labs(title = "", x = "MAG", y = "Number of identified proteins\n in biofilm") +
  scale_fill_manual(
    values = genus_colors, name = "Genus",
    labels = c(
      expression(italic("Ferroplasma")),
      expression(italic("Mycobacterium")),
      expression(italic("Acidithiobacillus")),
      "Others"
    )
  ) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_set(theme_classic(base_size = 20, base_family = "Arial")) +
  theme(
    # text             = element_text(size = 24),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.5),
    # panel.grid.minor = element_line(color = "grey92", linewidth = 0.5),
    panel.border     = ggplot2::element_rect(color = "grey85", fill = NA, linewidth = 0.8),
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.5),
    axis.text        = element_text(color = "black"),
    # axis.text.x      = element_text(size = 8)
  )

fig3b

################################################
# fig3c
################################################
#combine proteomics and metagenomics 
# get abundance from figure01 mean_cov_all
mean_cov_all <- read.csv("metagenomics/DNA_mean_cov_all.csv")[, -1]

data_proteomics_formated <- data_proteomics1 %>%
  group_by(ID) %>%
  summarise(
    sum_log10_LFQ_intensity_Cave_Biofilm_1  = sum(log10_LFQ_intensity_Cave_Biofilm_1),
    sum_log10_LFQ_intensity_Cave_Biofilm_2  = sum(log10_LFQ_intensity_Cave_Biofilm_2),
    sum_log10_LFQ_intensity_Culture         = sum(log10_LFQ_intensity_Culture),
    mean_log10_LFQ_intensity_Cave_Biofilm_1 = mean(log10_LFQ_intensity_Cave_Biofilm_1),
    mean_log10_LFQ_intensity_Cave_Biofilm_2 = mean(log10_LFQ_intensity_Cave_Biofilm_2),
    mean_log10_LFQ_intensity_Culture        = mean(log10_LFQ_intensity_Culture),
    sum_LFQ_intensity_Cave_Biofilm_1        = sum(LFQ_intensity_Cave_Biofilm_1),
    sum_LFQ_intensity_Cave_Biofilm_2        = sum(LFQ_intensity_Cave_Biofilm_2),
    sum_LFQ_intensity_Culture               = sum(LFQ_intensity_Culture)  # fixed: was summing the log10 column under a non-log name
  )

data_proteomics_formated$mean_mean_log10_LFQ_intensity_Cave_Biofilm <- rowMeans(data_proteomics_formated[, 5:6])

data_comb <- cbind(
  mean_cov_all[match(c("SC_MAG_00006", "SC_MAG_00016", "SC_MAG_00008", "SC_MAG_00004"), mean_cov_all$bins), ],
  data_proteomics_formated[match(c("SCMAG00006", "KDJLIKBO", "KNPMNEEE", "SCMAG00004"), data_proteomics_formated$ID), ]
)

data_comb$DNA_mean_ab <- rowMeans(data_comb[, 5:6])
data_comb$highlight   <- c("MAG 6", "MAG 16", "MAG 8", "MAG 4")
data_comb$genus       <- factor(data_comb$genus, levels = c("Ferroplasma", "Mycobacterium", "Acidithiobacillus"))

# Same genus->color mapping used in panel a, for visual consistency across the figure
genus_colors <- c(
  "Ferroplasma"       = '#ffd92f',
  "Mycobacterium"     = '#66c2a5',
  "Acidithiobacillus" = '#fc8d62',
  "Cuniculiplasma"    = '#e5c494',
  "Others"            = "grey"
)


# Panel c: mean LFQ intensity vs mean DNA coverage, colored by genus, labeled by organism
# Log-paper grid
# ============================================================================
data_comb$DNA_mean_ab <- rowMeans(data_comb[,5:6])
data_comb$highlight   <- c("Ferroplasma MAG6", "M. methanotrophicum", "Ferroplasma c.", "Acidithiobacillus")
data_comb$genus       <- factor(data_comb$genus)

fig3c <- ggplot(data_comb, aes(DNA_mean_ab, mean_mean_log10_LFQ_intensity_Cave_Biofilm, color = genus)) +
  geom_point(size = 8) +
  scale_color_manual(values = genus_colors ) +
  geom_label_repel(aes(label = highlight), size = 6, color = "black", fill = alpha(c("white"), 0.5)) +
  labs(x = "mean DNA coverage (log10)", y = "mean LFQ intenisty (log10)") +
  theme(text = element_text(size = 14))

fig3c

library(ggplot2)
library(ggrepel)
library(scales)

library(ggplot2)
library(ggrepel)
library(scales)

fig3c <- ggplot(
  data_comb,
  aes(
    x = DNA_mean_ab,
    y = mean_LFQ_intensity,
    color = genus
  )
) +
  geom_point(size = 8) +
  geom_label_repel(
    aes(label = highlight),
    size = 6,
    color = "black",
    fill = alpha("white", 0.5)
  ) +
  scale_color_manual(values = genus_colors) +
  
  scale_x_log10(
    breaks = 10^(2:7),
    labels = label_math(10^.x)
  ) +
  scale_y_log10(
    breaks = 10^(2:12),
    labels = label_math(10^.x)
  ) +
  
  annotation_logticks(sides = "lb") +
  
  labs(
    x = "Mean DNA coverage",
    y = "Mean LFQ intensity"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(colour = "grey80"),
    panel.grid.minor = element_line(colour = "grey92"),
    legend.title = element_blank()
  )

fig3c



