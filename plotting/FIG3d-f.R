setwd('/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/')
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
