# =============================================================================
# Fig. 3 — metaproteomic overview
# Figure numbering follows the manuscript (Cave_Ferro_Virus, Nov 2026).
# Part 1 (from Figure03a-c.R): a = 3D LFQ scatter, plus the protein-group bar
#   (object fig3b = manuscript 3c) and LFQ-vs-coverage (object fig3c = manuscript 3e).
# Part 2 (from FIG3d-f.R): pathway panels; objects fig3d/fig3e/fig3f are
#   manuscript 3b (Ca. M. methanotrophicum), 3d (Ferroplasma c. MAG 8), 3f (MAG 6).
# Object names follow an older panel layout; see figures/README.md.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root
setwd(REPO)

# ############################ Part 1 ############################
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

setwd(REPO)
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

dev.copy2pdf(file = file.path(REPO, "figures/main/Fig3/drafts", "fig3a.pdf"))

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





# ############################ Part 2 ############################
setwd(REPO)
# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

setwd(REPO)
library(ggplot2)

# Step 1: Load and process raw proteomics data
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

data_proteomics1 <- data_proteomics[which(data_proteomics$LFQ_intensity_Cave_Biofilm_1 > 0 & data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 > 0),]
data_proteomics1$ID <- unlist(lapply(strsplit(data_proteomics1$`Protein IDs`, "_"), function(x) x[1]))

# Step 2: Subset proteomics by organism, using each organism's unique protein-ID prefix
data_proteomics_ferro  <- data_proteomics1[grepl("KNPMNEEE",  data_proteomics1$`Majority protein IDs`),]
data_proteomics_mag6   <- data_proteomics1[grepl("SCMAG00006", data_proteomics1$`Majority protein IDs`),]
data_proteomics_M_meth <- data_proteomics1[grepl("KDJLIKBO",  data_proteomics1$`Majority protein IDs`),]

# Step 3: Define the CD-HIT cluster parser function
library(dplyr)
library(tidyr)

read.cdhit.clstr <- function(fname) {
  data.fields <- c("E.Value", "Aln", "Identity")
  read.table(fname, sep = "\t", comment.char = "", quote = "", fill = T, stringsAsFactors = F, col.names = c("Col1", "Col2")) %>%
    separate(Col1, into = c("Seq.Num", "Cluster"), sep = " ", fill = "right") %>%
    fill(Cluster) %>%
    filter(!grepl(">", Seq.Num)) %>%
    separate(Col2, into = c("Seq.Len", "Col2"), sep = "aa, >") %>%
    extract(Col2, into = c("Seq.Name", "Is.Representative", "Col2"), regex = "(.*?)[.]{3} ([*]|at) ?(.*)") %>%
    mutate(Is.Representative = Is.Representative == "*", Col2 = ifelse(Is.Representative, "100%", Col2)) %>%
    group_by(Cluster) %>%
    mutate(Representative = Seq.Name[which(Is.Representative)]) %>%
    separate_rows(Col2, sep = ",") %>%
    separate(Col2, into = data.fields, sep = "/", fill = "left", convert = T) %>%
    mutate(Identity = sub("%", "", Identity) %>% as.numeric) %>%
    group_by(Seq.Name) %>%
    mutate(level.rank = paste0(".", 1:n() - 1), level.rank = ifelse(level.rank == ".0", "", level.rank)) %>%
    pivot_wider(names_from = level.rank, values_from = data.fields, names_sep = "") %>%
    ungroup
}

# Step 4: Load pan-genome gene cluster summaries, subset to each genome
Ferroplasma_Pan_2_gene_clusters_summary.txt <- read.delim("pangenomics/Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
Ferr <- Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name == "SFerroplasmacircular"),]
Ferr <- Ferr[order(Ferr$gene_callers_id, decreasing = F),]

MAG6 <- Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name == "SMAG00006"),]
MAG6 <- MAG6[order(MAG6$gene_callers_id, decreasing = F),]

Myco_Pan_gene_clusters_summary.txt <- read.delim("pangenomics/MYCO_Pan_gene_clusters_summary.txt.gz")
Myco_Pan <- Myco_Pan_gene_clusters_summary.txt[which(Myco_Pan_gene_clusters_summary.txt$genome_name == "M_methanotrophicum"),]
Myco_Pan <- Myco_Pan[order(Myco_Pan$gene_callers_id, decreasing = F),]

# Step 5: Match CD-HIT clusters to pan-genome gene IDs, merge with proteomics — for each of the three organisms
# Ferroplasma c.
temp_Fer <- read.cdhit.clstr("proteomics/Ferroplasma_connect_100.clstr")

matchIDs <- rep(NA, nrow(Ferr))
for(i in unique(temp_Fer$Cluster)){
  temp1 <- temp_Fer[which(i == temp_Fer$Cluster),]
  if(nrow(temp1) == 2){
    idx <- which(Ferr$gene_callers_id == as.integer(temp1$Seq.Name[1]))
    matchIDs[idx] <- temp1$Seq.Name[2]
  }
}
Ferr$matchIDs <- matchIDs

data_proteomics_ferro$Majority_protein_IDs_u <- unlist(lapply(strsplit(data_proteomics_ferro$`Majority protein IDs`, ";"),
                                                              function(x) substr(x[grep("KNPMNEEE", x)], 1, 14)[1]))

Ferr_pan_only_prot <- merge(Ferr, data_proteomics_ferro, by.x = "matchIDs", by.y = "Majority_protein_IDs_u", all = TRUE)
Ferr_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm <- rowMeans(Ferr_pan_only_prot[,c(32,33)], na.rm = T)
Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm <- rowMeans(Ferr_pan_only_prot[,c(29,30)], na.rm = T)

# Ferroplasma MAG6
temp_mag6 <- read.cdhit.clstr("proteomics/MAG6_connect_100.clstr")

matchIDs <- rep(NA, nrow(MAG6))
for(i in unique(temp_mag6$Cluster)){
  temp1 <- temp_mag6[which(i == temp_mag6$Cluster),]
  if(nrow(temp1) == 2){
    idx <- which(MAG6$gene_callers_id == as.integer(temp1$Seq.Name[1]))
    matchIDs[idx] <- temp1$Seq.Name[2]
  }
}
MAG6$matchIDs <- gsub("SC_MAG_00006_", "SCMAG00006_", matchIDs)

data_proteomics_mag6$Majority_protein_IDs_u <- unlist(lapply(strsplit(data_proteomics_mag6$`Majority protein IDs`, ";"),
                                                             function(x) x[grep("SCMAG00006_", x)][1]))

MAG6_pan2_only_prot <- merge(MAG6, data_proteomics_mag6, by.x = "matchIDs", by.y = "Majority_protein_IDs_u", all = TRUE)
MAG6_pan2_only_prot$log10_LFQ_intensity_Cave_Biofilm <- rowMeans(MAG6_pan2_only_prot[,32:33], na.rm = T)
MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm <- rowMeans(MAG6_pan2_only_prot[,c(29,30)], na.rm = T)

# M. methanotrophicum
temp_Myc <- read.cdhit.clstr("proteomics/Mycobacterium_connect_100.clstr")

matchIDs <- rep(NA, nrow(Myco_Pan))
for(i in unique(temp_Myc$Cluster)){
  temp1 <- temp_Myc[which(i == temp_Myc$Cluster),]
  if(nrow(temp1) == 2){
    idx <- which(Myco_Pan$gene_callers_id == as.integer(temp1$Seq.Name[1]))
    matchIDs[idx] <- temp1$Seq.Name[2]
  }
}
Myco_Pan$matchIDs <- matchIDs

data_proteomics_M_meth$Majority_protein_IDs_u <- unlist(lapply(strsplit(data_proteomics_M_meth$`Majority protein IDs`, ";"),
                                                               function(x) substr(x[grep("KDJLIKBO", x)], 1, 14)[1]))

Myco_pan_only_prot <- merge(Myco_Pan, data_proteomics_M_meth, by.x = "matchIDs", by.y = "Majority_protein_IDs_u", all = TRUE)
Myco_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm <- rowMeans(Myco_pan_only_prot[,32:33], na.rm = T)
Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm <- rowMeans(Myco_pan_only_prot[,29:30], na.rm = T)

# Step 6: Load the KEGG pathway/module reference table, filtered
pws_mod_names <- read.csv("metagenomics/Functional_annotation/pw_info.csv")
pws_mod_names_r <- pws_mod_names[-grep("Human|Organismal Systems", pws_mod_names$V5),]
pws_mod_names_r <- pws_mod_names_r[-grep(" - fly| - worm", pws_mod_names_r$V2),]

# Panel d — M. methanotrophicum
Myco_data <- pws_mod_names_r
match_vals <- apply(pws_mod_names_r, 1, function(x)
  sum(Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), Myco_pan_only_prot$KOfam_ACC)], na.rm = T))
NR_active <- apply(pws_mod_names_r, 1, function(x)
  length(which(Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), Myco_pan_only_prot$KOfam_ACC)] > 0)))

Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm <- match_vals
Myco_data$Myco_NR_active <- NR_active
Myco_data$org <- rep("M.meth", nrow(Myco_data))
Myco_data$Myco_log10_Sum_LFQ_intensity_Cave_Biofilm <- log10(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm)
Myco_data <- Myco_data[-which(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm == 0),]
Myco_data <- Myco_data[order(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm, decreasing = T),]
Myco_data$pw_nr_active_perc <- Myco_data$Myco_NR_active / as.integer(Myco_data$V3) * 100

sub_Myco_data <- Myco_data[1:10,]
sub_Myco_data$V2 <- factor(sub_Myco_data$V2, levels = sub_Myco_data$V2)

sub_Myco_data$V2 <- gsub("Pentose phosphate pathway, non-oxidative phase, fructose 6P => ribose 5P", "Pentose phosphate pathway, \n non-oxidative phase", sub_Myco_data$V2)
sub_Myco_data <- sub_Myco_data[order(sub_Myco_data$Myco_log10_Sum_LFQ_intensity_Cave_Biofilm, decreasing = T),]

# log-spaced minor gridline positions: 1,2,...,9 × 10^n
log_minor <- as.vector(outer(1:9, 10^(0:12)))

fig3d <- ggplot(sub_Myco_data,
                aes(Myco_Sum_LFQ_intensity_Cave_Biofilm,
                    reorder(V2, -Myco_Sum_LFQ_intensity_Cave_Biofilm),
                    color = Myco_NR_active)) +
  geom_point(size = 6) +
  scale_x_log10(
    breaks       = c(1e9, 1e10),
    labels       = scales::label_log(),
    minor_breaks = log_minor
  ) +
  coord_cartesian(xlim = c(1e9, 1e11)) +
  labs(
    y = "",
    x = "",
    title = expression("C. " * italic(M) * ". methanotrophicum"),
    color = "Number of proteins"
  ) +
  scale_colour_viridis_c(option = "plasma") +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.8),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    axis.text.y = element_text(
      color = "black",
      size = 18
    ),
    axis.text.x = element_text(
      color = "black",
      size = 18
    )
  )

fig3d

# Panel e — Ferroplasma c.
Ferr_data <- pws_mod_names_r
match_vals <- apply(pws_mod_names_r, 1, function(x)
  sum(Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), Ferr_pan_only_prot$KOfam_ACC)], na.rm = T))
NR_active <- apply(pws_mod_names_r, 1, function(x)
  length(which(Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), Ferr_pan_only_prot$KOfam_ACC)] > 0)))

Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm <- match_vals
Ferr_data$Ferr_NR_active <- NR_active
Ferr_data$org <- rep("Ferroplasma c.", nrow(Ferr_data))
Ferr_data$Ferr_log10_Sum_LFQ_intensity_Cave_Biofilm <- log10(Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm)
Ferr_data <- Ferr_data[-which(Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm == 0),]
Ferr_data <- Ferr_data[order(Ferr_data$Ferr_log10_Sum_LFQ_intensity_Cave_Biofilm, decreasing = T),]
Ferr_data$pw_nr_active_perc <- Ferr_data$Ferr_NR_active / as.integer(Ferr_data$V3) * 100

sub_Ferr_data <- Ferr_data[1:10,]
sub_Ferr_data$V2 <- factor(sub_Ferr_data$V2, levels = sub_Ferr_data$V2)

fig3e <- ggplot(
  sub_Ferr_data,
  aes(Ferr_Sum_LFQ_intensity_Cave_Biofilm, V2, color = Ferr_NR_active)
) +
  geom_point(size = 6) +
  scale_x_log10(
    breaks       = c(1e9, 1e10),
    labels       = scales::label_log(),
    minor_breaks = log_minor
  ) +
  coord_cartesian(xlim = c(1e9, 1e11)) +
  labs(
    y = "",
    x = "",
    title = expression("" * italic(Ferroplasma) * " c."),
    # title = "Ferroplasma c. MAG 8",
    color = "Number of proteins"
  ) +
  scale_colour_viridis_c(option = "plasma") +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.8),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    axis.text.y = element_text(
      color = "black",
      size = 18
    ),
    axis.text.x = element_text(
      color = "black",
      size = 18
    )
  )

fig3e

# Panel f — Ferroplasma MAG6
MAG6_data <- pws_mod_names_r
match_vals <- apply(pws_mod_names_r, 1, function(x)
  sum(MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), MAG6_pan2_only_prot$KOfam_ACC)], na.rm = T))
NR_active <- apply(pws_mod_names_r, 1, function(x)
  length(which(MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5], ",")), MAG6_pan2_only_prot$KOfam_ACC)] > 0)))

MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm <- match_vals
MAG6_data$MAG6_NR_active <- NR_active
MAG6_data$org <- rep("Ferroplasma MAG6", nrow(MAG6_data))
MAG6_data$MAG6_log10_Sum_LFQ_intensity_Cave_Biofilm <- log10(MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm)
MAG6_data <- MAG6_data[-which(MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm == 0),]
MAG6_data <- MAG6_data[order(MAG6_data$MAG6_log10_Sum_LFQ_intensity_Cave_Biofilm, decreasing = T),]
MAG6_data$pw_nr_active_perc <- MAG6_data$MAG6_NR_active / as.integer(MAG6_data$V3) * 100

sub_MAG6_data <- MAG6_data[1:10,]
sub_MAG6_data$V2 <- factor(sub_MAG6_data$V2, levels = sub_MAG6_data$V2)

fig3f <- ggplot(
  sub_MAG6_data,
  aes(MAG6_Sum_LFQ_intensity_Cave_Biofilm, V2, color = MAG6_NR_active)
) +
  geom_point(size = 6) +
  scale_x_log10(
    breaks       = c(1e9, 1e10),
    labels       = scales::label_log(),
    minor_breaks = log_minor
  ) +
  coord_cartesian(xlim = c(1e9, 1e11)) +
  labs(
    y = "",
    x = "",
    title = expression("" * italic(Ferroplasma) * " MAG 6"),
    color = "Number of proteins"
  ) +
  scale_colour_viridis_c(option = "plasma") +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.8),
    axis.line    = element_line(color = "black", linewidth = 0.8),
    legend.position = "none",
    axis.text.y = element_text(
      color = "black",
      size = 18
    ),
    axis.text.x = element_text(
      color = "black",
      size = 18
    )
  )

fig3d
fig3e
fig3f

library(cowplot)
library(ggplot2)
library(viridis)

# Build fig3d as you already have it (legend included, so don't set legend.position = "none")
fig3d_forlegend <- ggplot(sub_Myco_data, aes(Myco_log10_Sum_LFQ_intensity_Cave_Biofilm, V2, color = Myco_NR_active)) +
  geom_point(size = 6) +
  labs(
    y = "Top 10 pathways",
    x = "",
    title = expression("C. " * italic(M) * ". methanotrophicum MAG 16"),
    color = "Number of proteins"
  ) +
  scale_colour_viridis_c(option = "plasma") +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_text(color = "black", size = 14),
    axis.text.x = element_text(color = "black", size = 14)
  )

fig3d_forlegend
# Extract just the legend grob
fig3d_legend <- get_legend(fig3d_forlegend)

# Draw it alone on a blank canvas
legend_only <- ggdraw(fig3d_legend)
legend_only
