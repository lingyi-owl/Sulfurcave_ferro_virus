#Figure2
#pan-genome ferroplasma
library(ggplot2)

Ferroplasma <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
Myco<- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/MYCO_Pan_gene_clusters_summary.txt.gz")
CAVE <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/CAVE_R_Pan_gene_clusters_summary.txt.gz")
hist(Ferroplasma$num_genomes_gene_cluster_has_hits/13*100)
hist(Myco$num_genomes_gene_cluster_has_hits/69*100)
hist(CAVE$num_genomes_gene_cluster_has_hits/21*100)

#################################################################################
##unique statistics
#################################################################################

CAVE_f<-CAVE[CAVE$genome_name == "SFerroplasmacircularcontig",]
table(CAVE_f$num_genomes_gene_cluster_has_hits)
annotations<- CAVE_f$KOfam_ACC
# Calculate percentage
total_count <- length(annotations)
annotated_count <- sum(annotations != "")
percentage_annotated <- (annotated_count / total_count) * 100
percentage_annotated
# 45.15%
# 1   2   3   4   5   6   7   8  10  11  12  13  14  15  18  19  20  21 
# 479 893 293 414   8   9   4   6   2   2   1   2   1   4   1   1   2   2 
annotations<- CAVE_f$KOfam_ACC[CAVE_f$num_genomes_gene_cluster_has_hits==1]
# Calculate percentage
total_count <- length(annotations)
annotated_count <- sum(annotations != "")
percentage_annotated <- (annotated_count / total_count) * 100
percentage_annotated
# 16.07
CAVE_fMAG6<-CAVE[CAVE$genome_name == "SMAG00006",]
table(CAVE_fMAG6$num_genomes_gene_cluster_has_hits)
annotations<- CAVE_fMAG6$KOfam_ACC
# Calculate percentage
total_count <- length(annotations)
annotated_count <- sum(annotations != "")
percentage_annotated <- (annotated_count / total_count) * 100
percentage_annotated
# 40.13
# 1   2   3   4   5   6   7   8  10  11  12  13  15  18  19  20  21 
# 775 870 280 456  12  26   5   2   3   1   1   3   3   1   1   2   6 
annotations<- CAVE_fMAG6$KOfam_ACC[CAVE_fMAG6$num_genomes_gene_cluster_has_hits==1]
# Calculate percentage
total_count <- length(annotations)
annotated_count <- sum(annotations != "")
percentage_annotated <- (annotated_count / total_count) * 100
percentage_annotated
#22.06

#################################################################################
#
#################################################################################

x1<-data.frame(gene_cluster_freq= Ferroplasma$num_genomes_gene_cluster_has_hits,gene_cluster_freq_perc=Ferroplasma$num_genomes_gene_cluster_has_hits/13*100,Pan_genome=rep("Ferroplasma",nrow(Ferroplasma)))
x2<-data.frame(gene_cluster_freq= Myco$num_genomes_gene_cluster_has_hits,gene_cluster_freq_perc=Myco$num_genomes_gene_cluster_has_hits/69*100,Pan_genome=rep("Mycobacterium",nrow(Myco)))
x3<-data.frame(gene_cluster_freq= CAVE$num_genomes_gene_cluster_has_hits,gene_cluster_freq_perc=CAVE$num_genomes_gene_cluster_has_hits/21*100,Pan_genome=rep("Cave community",nrow(CAVE)))
xx<-rbind(x1,x2,x3)
library(ggridges)
ggplot(xx, aes(x = gene_cluster_freq_perc, y = Pan_genome, fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1)

colors_r<-c('#fc8d62','#8da0cb','#66c2a5')
ggplot(xx, aes(x = gene_cluster_freq_perc, y = Pan_genome, fill = Pan_genome)) +
  geom_density_ridges(alpha = .6) +labs(x="Gene cluster frequencies (%)",y="Pan-genomes",fill= "")+
  theme(text = element_text(size = 24)) +scale_fill_manual(values=colors_r) +theme(legend.position = "none")
ggplot(xx, aes(x = gene_cluster_freq_perc, fill = Pan_genome)) +
  geom_histogram(binwidth = 3, alpha = 0.7) +labs(x="Gene cluster frequencies (%)",y="Pan-genomes",fill= "")+
  theme(text = element_text(size = 24))+ facet_wrap(~ Pan_genome, ncol = 1) +scale_fill_manual(values=colors_r) +theme(legend.position = "none")

ferroplasma_table<- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/ferroplasma_table.txt")
ferroplasma_table$groups<-rep("NA",nrow(ferroplasma_table))
Group1<-c("SMAG00006","S25910037","S259100318","S259100328")
Group2<-c("S259100341","SFerroplasmacircular" , "S33314612","S749695")
Group3<-c("S259100335","S259100337","S259100333","S25910039")
ferroplasma_table$groups[ferroplasma_table$layer%in%Group1]<-"Group1"
ferroplasma_table$groups[ferroplasma_table$layer%in%Group2]<-"Group2"
ferroplasma_table$groups[ferroplasma_table$layer%in%Group3]<-"Group3"
ferroplasma_tablef<-ferroplasma_table[!ferroplasma_table$groups=="NA",]

####### perform anivo functional enrichment
Ferro_cave_un_enriched.KEGG_Module <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/Ferro_cave_un_enriched-KEGG_Module.txt")
Ferro_cave_un_enriched.KEGG_Module$category <- "KEGG_Module"
colnames(Ferro_cave_un_enriched.KEGG_Module)[1] <- "Function"
Ferro_cave_un_enriched.KOfam <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/Ferro_cave_un_enriched-KOfam.txt")
Ferro_cave_un_enriched.KOfam$category <- "KOfam"
colnames(Ferro_cave_un_enriched.KOfam)[1] <- "Function"
data_enrichment <- rbind(Ferro_cave_un_enriched.KEGG_Module[,c(1:7,9,14)],Ferro_cave_un_enriched.KOfam[,c(1:7,10,14)])

# Create scatterplot
data_enrichment$Significance <- ifelse(data_enrichment$adjusted_q_value < 0.1 & data_enrichment$p_rest ==0, "Yes", "No")

names_to_label <- c(
  "phosphoenolpyruvate---glycerone phosphotransferase", 
  "CRISPR", 
  "glycerol uptake facilitator protein"
)


data_enrichment$Function_highlight<-data_enrichment$Function
data_enrichment$Function_highlight[-unlist(lapply(names_to_label, function(x) grep(x,data_enrichment$Function)))] <- NA
data_enrichment$Function_highlight<- unlist(lapply(strsplit(data_enrichment$Function_highlight,"-"), function(x) x[1]))

data_enrichment$pvalue<- - log10(data_enrichment$adjusted_q_value)

colors_s<-c('#f1a340','#8da0cb')
f2b<-ggplot(data_enrichment, aes(x = enrichment_score, y = pvalue, color = Significance, shape = category)) +
  geom_jitter(size = 5, width = 0.3, height = 0.3,alpha=0.7) +
  theme_minimal() +
  labs(
    title = "",
    x = "Enrichment score",
    y = "-Log10(pvalue)",
    color = "Significance",
    shape = "Category"
  ) + theme(text = element_text(size = 24)) + scale_color_manual(values = colors_s)
f2b

KOs_ferro_in_cave <- data_enrichment$accession[which(data_enrichment$Significance=="Yes"& data_enrichment$category=="KOfam")]
# +
#   geom_text_repel(
#     data = data_enrichment[unlist(lapply(names_to_label, function(x) grep(x,data_enrichment$Function))), ],
#     aes(label = unlist(lapply(strsplit(Function,"-"), function(x) x[1]))),
#     size = 8,
#     box.padding = 0.5,
#     point.padding = 0.3, max.overlaps=50
#   ) 
data_table_sig<-data_enrichment[data_enrichment$Significance=="Yes",]
data_table_sig$enrichment_score <- round(data_table_sig$enrichment_score,digits = 3)
data_table_sig$unadjusted_p_value <- round(data_table_sig$unadjusted_p_value,digits = 3)
data_table_sig$adjusted_q_value <- round(data_table_sig$adjusted_q_value,digits = 3)
#write.csv(data_table_sig,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/cave/figure2b_sig_KOs.csv")
#cave function plots:
##################################################################
# read all eggnog and make a matrix of 1 and 0 for genomes and KOs
##################################################################

KEGG_KO_list<-list()
OG_list<-list()
i<-1
for(j in unique(CAVE$genome_name)){
  KEGG_KO_list[[i]]<- unique(CAVE$KOfam_ACC[which(CAVE$genome_name==j)])
  OG_list[[i]]<- unique(CAVE$gene_cluster_id[which(CAVE$genome_name==j)])
  i<-i+1
}
names(KEGG_KO_list)<-unique(CAVE$genome_name)
names(OG_list)<-unique(CAVE$genome_name)

KEGG_u_KO<-unique(CAVE$KOfam_ACC)
KEGG_u_KO<-KEGG_u_KO[-1]
KO_matrix<-c()
for (i in 1:length(KEGG_KO_list)){
  t<-lapply(KEGG_u_KO, function(x)  ifelse(any(grepl(x,KEGG_KO_list[[i]])), 1, 0))
  KO_matrix<-rbind(KO_matrix, unlist(t))
}
colnames(KO_matrix)<-KEGG_u_KO
rownames(KO_matrix)<-unique(CAVE$genome_name)
KO_matrix<-as.matrix(KO_matrix)

OG_u<-unique(CAVE$gene_cluster_id)
OG_matrix<-c()
for (i in 1:length(OG_list)){
  t<-lapply(OG_u, function(x)  ifelse(any(grepl(x,OG_list[[i]])), 1, 0))
  OG_matrix<-rbind(OG_matrix, unlist(t))
}
colnames(OG_matrix)<-OG_u
rownames(OG_matrix)<-unique(CAVE$genome_name)
OG_matrix<-as.matrix(OG_matrix)
#write.csv(OG_matrix,"/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/OG_matrix_cave.csv")
OG_matrix <- read.csv("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/OG_matrix_cave.csv", row.names = "X")
#write.csv(matrix,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/KO_matrix_cave.csv")
KO_matrix <- read.csv("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/KO_matrix_cave.csv", row.names = "X")

########
library(UpSetR)


# Define a vector of colors for the main bars
bar_colors <- c("#66c2a5", "#66c2a5", "#66c2a5", "#8da0cb", "#66c2a5", "#8da0cb", "#8da0cb", "#66c2a5", "#66c2a5", "#66c2a5")
rownames(KO_matrix)[c(2,8,19,20,21)] <- c("F. c.","F. MAG6","M. MAG2", "M.MAG3" , "C. M. methanotrophicum")
# Create the UpSet plot with differently colored bars
data_upset <- upset(
  as.data.frame(t(KO_matrix)),
  main.bar.color = bar_colors,
  sets.bar.color = "black",
  order.by = "freq",
  nintersects = NA,
  nsets = 21,
  sets =  c("F. c.","F. MAG6","M. MAG2", "M.MAG3" , "C. M. methanotrophicum"),
  text.scale = 2
)
data_upset
library(ggplotify)
f2c <- as.ggplot(data_upset)
#ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/figure2c.pdf",f2c, width =16, height = 5)
###########
#run figure 1 to get mean_cov_all
##########
matrix<-KO_matrix[-which(rownames(KO_matrix)=="SMAG00003"),]

# Load necessary libraries
library(ggplot2)
library(factoextra)

# Perform PCA
# Function to remove columns with all zeros or all ones
remove_all_zeros_ones_columns <- function(matrix) {
  non_constant_columns1 <- apply(matrix, 2, function(col) all(col ==1))
  return(matrix[, !non_constant_columns1, drop = FALSE])
}

# Remove columns with all zeros or all ones
filtered_matrix <- remove_all_zeros_ones_columns(matrix)


pca_result1 <- prcomp(filtered_matrix)

# Visualize PCA using a biplot
fviz_pca_ind(pca_result1)

# Print the summary of the PCA
summary(pca_result1)

# ####################
# 
# library(ComplexUpset)
# upset(CAVE, genome_name, name='genome_name', width_ratio=0.1)
# # Transform the data for ggupset
# data_long <- CAVE %>%
#   group_by(gene_cluster_id) %>%
#   summarize(genomes = list(genome_name))
# 
# # Create the ggupset plot
# ggplot(data_long, aes(x = genomes)) +
#   geom_bar() +
#   scale_x_upset() +
#   labs(
#     title = "Gene Cluster Associations with Genomes",
#     x = "Genomes",
#     y = "Count"
#   ) +
#   theme_minimal()

############# calculate pathway completeness

# Load required library
# Load required packages for data manipulation
library(dplyr)  # for data manipulation (filter, group_by, summarise, etc.)
library(tidyr)  # for splitting and reshaping data

# -------------------------------
# STEP 1: Prepare unique KEGG ortholog data
# -------------------------------

# Start with your main dataset 'CAVE'
# separate_rows() splits KEGG_Module_ACC if multiple modules are in one cell, separated by '|'
# select() keeps only the relevant columns
# distinct() removes duplicate rows (same KO-module-genome combination)
unique_kos <- CAVE %>%
  separate_rows(KEGG_Module_ACC, sep = "\\|") %>%  # Split multiple modules in one cell
  select(KOfam_ACC, KEGG_Module_ACC, genome_name)  # Keep relevant columns
distinct()  # Remove exact duplicates

# -------------------------------
# STEP 2: Count unique KOfam per KEGG module for each genome
# -------------------------------

result <- unique_kos %>%
  group_by(genome_name, KEGG_Module_ACC) %>%  # Group by genome and module
  summarise(unique_KOs = n_distinct(KOfam_ACC),  # Count unique KOfam IDs
            .groups = "drop")  # Drop grouping after summarise

# Remove rows where KEGG_Module_ACC is empty
result <- result[-which(result$KEGG_Module_ACC == ""),]

# -------------------------------
# STEP 3: Calculate stats per KEGG module across genomes
# -------------------------------

result_SD <- result %>%
  group_by(KEGG_Module_ACC) %>%  # Group by KEGG module
  summarise(
    # Perform Kruskal-Wallis test across genomes if more than 1 genome
    p_value = if(length(unique(genome_name)) > 1) {
      kruskal.test(unique_KOs ~ genome_name)$p.value
    } else {
      NA  # If only one genome, test cannot be performed
    },
    SD_unique_KOs = sd(unique_KOs, na.rm = TRUE),  # Standard deviation
    avg_unique_KOs = mean(unique_KOs, na.rm = TRUE),  # Overall mean
    # Average for Ferroplasma genomes
    avg_unique_KOs_ferroplasma = mean(unique_KOs[genome_name %in% c("SFerroplasmacircularcontig","SMAG00006")], na.rm = TRUE),
    # Average for all other genomes
    avg_unique_KOs_other = mean(unique_KOs[!genome_name %in% c("SFerroplasmacircularcontig","SMAG00006")], na.rm = TRUE),
    # Logical: is this module more abundant in Ferroplasma than others?
    is_ferroplasma = avg_unique_KOs_ferroplasma > avg_unique_KOs_other,
    .groups = "drop"
  )

# Replace NaN values (from empty averages) with 0
result_SD$avg_unique_KOs_ferroplasma[is.nan(result_SD$avg_unique_KOs_ferroplasma)] <- 0
result_SD$avg_unique_KOs_other[is.nan(result_SD$avg_unique_KOs_other)] <- 0

# Replace NA logicals with "NA" string to avoid plotting issues
result_SD$is_ferroplasma[is.na(result_SD$is_ferroplasma)] <- "NA"

# -------------------------------
# STEP 4: Retrieve KEGG module information
# -------------------------------
library(KEGGREST)

pw_info <- list()  # Initialize empty list to store KEGG info
for(i in 1:length(result_SD$KEGG_Module_ACC)){
  tryCatch({
    # Retrieve module info from KEGG database
    pw_info[[i]] <- keggGet(result_SD$KEGG_Module_ACC[i])
  }, error=function(e){
    message('An Error Occurred')  # Handle errors gracefully
    print(e)
    return(NA)
  })
}

# -------------------------------
# STEP 5: Extract module names, class, pathways, and orthology counts
# -------------------------------
pws_names <- c()
for(i in 1:length(pw_info)){
  temp <- c(
    pw_info[[i]][[1]]$ENTRY,  # KEGG Module ID
    pw_info[[i]][[1]]$NAME,   # Module name
    pw_info[[i]][[1]]$CLASS,  # Functional class
    unlist(paste(unlist(names(pw_info[[i]][[1]]$PATHWAY)), collapse = ",")),  # Pathway IDs
    unlist(paste(unlist(pw_info[[i]][[1]]$PATHWAY), collapse = ",")),        # Pathway names
    length(unique(pw_info[[i]][[1]]$ORTHOLOGY))  # Number of orthologs in module
  )
  pws_names <- rbind(pws_names, temp)
}

# Convert to data.frame and assign column names
pws_names <- as.data.frame(pws_names)
colnames(pws_names) <- c("ENTRY","NAME","CLASS","PATHWAY-ID","PATHWAY","NR-ORTHOLOGY")

# -------------------------------
# STEP 6: Calculate ratio and merge KEGG info
# -------------------------------
result_SD$ratio <- result_SD$avg_unique_KOs_ferroplasma / result_SD$avg_unique_KOs_other

# Merge with KEGG module metadata
result_SD_pw <- merge(result_SD, pws_names, by.x = "KEGG_Module_ACC", by.y = "ENTRY", all.x = TRUE)

# Simplify module class for plotting
result_SD_pw$class <- unlist(lapply(strsplit(result_SD_pw$CLASS, ";"), function(x) x[2]))
result_SD_pw$class[is.na(result_SD_pw$class)] <- "Energy metabolism"
result_SD_pw$SD_unique_KOs[is.na(result_SD_pw$SD_unique_KOs)] <- 0

# Optionally, save or read CSV (commented out or read from previous export)
# write.csv(result_SD_pw,"/path/to/result_SD_pw.csv")
result_SD_pw <- read.csv("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/result_SD_pw.csv")

# -------------------------------
# STEP 7: Plot data with ggplot2 + beeswarm + labels
# -------------------------------
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)

colors_kegg <- c('#7570b3','#66c2a5')  # Define custom colors
result_SD_pw$is_ferroplasma <- as.factor(result_SD_pw$is_ferroplasma)  # Convert to factor for coloring

f2c <- ggplot(result_SD_pw, aes(x = SD_unique_KOs, y = class)) +
  geom_beeswarm(aes(color = is_ferroplasma, size = avg_unique_KOs)) +  # Plot points with beeswarm layout
  labs(
    size = "Average KOs",
    x = "Std Dev KOs",
    y = "KEGG Module Functional Class",
    color = expression(italic("Ferroplasma"))
  ) +
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = colors_kegg)
# +
#   geom_label_repel(
#     data = subset(result_SD_pw, SD_unique_KOs > 3 & is_ferroplasma == TRUE),  # Only label high-variance Ferroplasma modules
#     aes(label = KEGG_Module_ACC),
#     max.overlaps = Inf  # Ensure all labels are plotted
#   )
#################################
# new plot with new legen 2026-05-18
# 1. Load the data
result_SD_pw <- read.csv("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/pangenomics/result_SD_pw.csv")
# Step 1: convert logical to factor directly using logical values
result_SD_pw$is_ferroplasma <- factor(result_SD_pw$is_ferroplasma, 
                                      levels = c(TRUE, FALSE),
                                      labels = c("Ferroplasma", "Not Ferroplasma"))

# Step 2: plot — NA values stay as NA and are handled by na.value
f2c <- ggplot(result_SD_pw, aes(x = SD_unique_KOs, y = class)) +
  geom_beeswarm(aes(color = is_ferroplasma, size = avg_unique_KOs)) +
  labs(
    size = "Average KOs",
    x = "Std Dev KOs",
    y = "KEGG Module Functional Class",
    color = "Taxonomy"
  ) +
  theme(text = element_text(size = 20)) +
  scale_color_manual(values = c("Ferroplasma"     = "#66c2a5",
                                "Not Ferroplasma" = "#7570b3"),
                     na.value = "grey50")


f2c

library(patchwork)
patchwork2 <- (f2a  + f2b) / (f2a  + f2b)
patchwork2 + plot_annotation(tag_levels = 'a')
ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/figure2.pdf", width =24, height = 10)


##########


# #########################
# #alternative
# ########################
# 
# library(KEGGREST)
# kegg_KOs<-colnames(matrix)[grep("ko",colnames(matrix))]
# #remove Global and overview maps
# to_remove<-c("ko01100","ko01110","ko01120","ko01200","ko01210","ko01212","ko01230","ko01232","ko01250","ko01240","ko01220")
# kegg_KOs<-kegg_KOs[!kegg_KOs%in%to_remove]
# pw_info<-list()
# for(i in 1:length(kegg_KOs)){
#   tryCatch({
#     pw_info[[i]]<- keggGet(kegg_KOs[i])
#     if(is.null(pw_info[[i]][[1]]$ORTHOLOGY)){
#       module<-list()
#       for (j in 1:length(pw_info[[i]][[1]]$MODULE)){
#         module[[j]]<-keggGet(names(pw_info[[i]][[1]]$MODULE)[j])
#       }
#       KOs<-lapply(module, function(x) x[[1]]$ORTHOLOGY)
#       names(KOs)<-lapply(module, function(x) x[[1]]$ENTRY)
#       pw_info[[i]][[1]]$ORTHOLOGY<-unlist(KOs)
#     }
#   }, error=function(e){
#     message('An Error Occurred')
#     print(e)
#     return(NA)
#   })
# }
# #save(pw_info, file = "kegg_KOs_list.RData")
# pw_info<-load( file = "kegg_KOs_list.RData")
# pw_info_c <- pw_info[sapply(pw_info, function(x) !is.null(x))]
# 
# pws_names<-c()
# for(i in 1:length(pw_info_c)){
#   temp<-c(
#     unlist(paste(unlist(names(pw_info_c[[i]][[1]]$PATHWAY)),collapse = ",")),
#     unlist(paste(unlist(pw_info_c[[i]][[1]]$PATHWAY),collapse = ",")),
#     length(unique(pw_info_c[[i]][[1]]$ORTHOLOGY))
#   )
#   #print(temp)
#   pws_names<-rbind(pws_names,temp)
# }
# pws_names<-as.data.frame(pws_names)
# colnames(pws_names)<-c("ID","PATHWAY","Number-KOs")
# pws_names$`Number-KOs`<-as.integer(pws_names$`Number-KOs`)
# rm<-which(pws_names$`Number-KOs`==0)
# pws_names<-pws_names[-rm,]
# pw_info_c<-pw_info_c[-rm]
# pws_names$ID
# 
# #count kos in pathways
# pw_matrix<-c()
# for (i in 1:length(pw_info_c)){
#   KOs<-paste("ko:",gsub("^M......","",names(pw_info_c[[i]][[1]]$ORTHOLOGY)),sep = "")
#   t<-unlist(lapply(KEGG_u$KO, function(x) sum(x%in%KOs)))
#   pw_matrix<-cbind(pw_matrix, t)
# }
# rownames(pw_matrix)<-gsub(".eggnog.emapper.annotations","",flist)
# colnames(pw_matrix)<-pws_names$ID
# 
# pw_matrix_comp<-sweep(pw_matrix,2,pws_names$`Number-KOs`,FUN="/")*100
# write.csv(pw_matrix_comp,"/home/chrats/Desktop/Projects/yeast_bacteria_interactions/chunxu/ANNOATIONS_F/BAC/PW_comp_matrix.csv")
# #calculate pathway completeness
# 
# sd_all<-apply(pw_matrix_comp, 2, function(x) sd(x))
# pw_matrix_compf<-pw_matrix_comp[,which(sd_all>quantile(sd_all)[4])] #sd > 75% quantile of all sd
# pw_matrix_compf<-pw_matrix_comp[,colnames(pw_matrix_comp)%in%names(sort(sd_all,decreasing = T)[1:20])] #top 20 in sd
# library(reshape2)
# data<-melt(pw_matrix_compf)
# 
# 
# ggplot(data,aes(Var2,Var1))+geom_point(aes(size=value))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 
# 
# # Hierarchical clustering of rows and columns
# row_order <- hclust(dist(pw_matrix_compf))$order
# col_order <- hclust(dist(t(pw_matrix_compf)))$order
# 
# # Create a ggplot heatmap with clustered rows and columns
# ggheatmap <- ggplot(data, aes(x = Var2, y = Var1, fill = value)) +
#   geom_tile() +
#   geom_text(aes(label = round(value, 2)), vjust = 1) +  # Display labels
#   scale_fill_gradient(low = "blue", high = "red") +  # Adjust color gradient as needed
#   labs(title = "Clustered Heatmap Example", x = "Columns", y = "Rows") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
#   scale_x_discrete(limits = colnames(pw_matrix_compf)[col_order]) +  # Apply column order
#   scale_y_discrete(limits = rownames(pw_matrix_compf)[row_order])    # Apply row order
# 
# # Print the ggplot heatmap
# print(ggheatmap)
# exploration of pan-genomics

library(tidyverse)

Ferroplasma <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/Ferroplasma/Ferroplasma_Pan_2/SUMMARY_few_classes/Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
Myco<- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/mycobacterium/MYCO-SUMMARY/MYCO_Pan_gene_clusters_summary.txt.gz")
CAVE <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/cave/CAVE_R_Pan/CAVE_R-SUMMARY/CAVE_R_Pan_gene_clusters_summary.txt.gz")

CAVE_KOs<-unique(CAVE$KOfam_ACC)
unique_CAVE_KEGG <- CAVE %>%
  select(KOfam_ACC,KOfam, KEGG_Module_ACC,  KEGG_Module  ) %>%
  distinct(KOfam_ACC,KOfam, KEGG_Module_ACC,  KEGG_Module)
# unique_CAVE_COG <- CAVE %>%
#   select( COG20_FUNCTION_ACC, COG20_FUNCTION, COG20_CATEGORY_ACC,COG20_CATEGORY   ) %>%
#   distinct()
unique_CAVE_KEGG_acc<-unique(unique_CAVE_KEGG$KOfam_ACC)
unique_CAVE_KEGG_acc<-unique_CAVE_KEGG_acc[-which(unique_CAVE_KEGG_acc=="")]
CAVE_KOs_prev<-c()
for(i in unique_CAVE_KEGG_acc){
  prev<-length(unique(CAVE$genome_name[which(CAVE$KOfam_ACC==i)]))
  CAVE_KOs_prev<-c(CAVE_KOs_prev,prev)
}
CAVE_KOs_prev<-CAVE_KOs_prev/21*100
cave_ko_prev<-cbind(unique_CAVE_KEGG_acc,CAVE_KOs_prev)

unique_Myco_KOs<-unique(Myco$KOfam_ACC)
unique_Myco_KOs<-unique_Myco_KOs[-which(unique_Myco_KOs=="")]
Myco_KOs_prev<-c()
for(i in unique_Myco_KOs){
  prev<-length(unique(Myco$genome_name[which(Myco$KOfam_ACC==i)]))
  Myco_KOs_prev<-c(Myco_KOs_prev,prev)
}
Myco_KOs_prev<-Myco_KOs_prev/69*100
Myco_ko_prev<-cbind(unique_Myco_KOs,Myco_KOs_prev)

cave_myco<-merge(cave_ko_prev,Myco_ko_prev,by.x="unique_CAVE_KEGG_acc",by.y="unique_Myco_KOs", all=T)

unique_Ferroplasma_KOs<-unique(Ferroplasma$KOfam_ACC)
unique_Ferroplasma_KOs<-unique_Ferroplasma_KOs[-which(unique_Ferroplasma_KOs=="")]

Ferr_KOs_prev<-c()
for(i in unique_Ferroplasma_KOs){
  prev<-length(unique(Ferroplasma$genome_name[which(Ferroplasma$KOfam_ACC==i)]))
  Ferr_KOs_prev<-c(Ferr_KOs_prev,prev)
}
Ferr_KOs_prev<-Ferr_KOs_prev/13*100
Ferr_KOs_prev<-cbind(unique_Ferroplasma_KOs,Ferr_KOs_prev)
cave_myco_ferro<-merge(cave_myco,Ferr_KOs_prev,by.x="unique_CAVE_KEGG_acc",by.y="unique_Ferroplasma_KOs", all=T)

cave_myco_ferro$CAVE_KOs_prev<-as.numeric(cave_myco_ferro$CAVE_KOs_prev)
cave_myco_ferro$Myco_KOs_prev<-as.numeric(cave_myco_ferro$Myco_KOs_prev)
cave_myco_ferro$Ferr_KOs_prev<-as.numeric(cave_myco_ferro$Ferr_KOs_prev)
cave_myco_ferrof <- cave_myco_ferro %>%
  mutate(across(everything(), ~replace_na(., 0)))


ggplot(cave_myco_ferrof, aes(x=CAVE_KOs_prev,y=Myco_KOs_prev))+geom_jitter()


KOs_ferro_in_cave
KO_Ferro<-CAVE$KOfam_ACC[grep("Amino acid transport and metabolism",CAVE$COG20_CATEGORY[CAVE$genome_name=="SFerroplasmacircularcontig"]) ]
KO_FMAG6<-CAVE$KOfam_ACC[grep("Amino acid transport and metabolism",CAVE$COG20_CATEGORY[CAVE$genome_name=="SMAG00006"]) ]
KO_FERRO_AA<- unique(c(KO_Ferro,KO_FMAG6))[-1]
KO_AA_ALL<-CAVE$KOfam_ACC[grep("Amino acid transport and metabolism",CAVE$COG20_CATEGORY) ]
KO_AA_ALL_clean<- unique(KO_AA_ALL)[-1]

library(KEGGREST)
AA_map<-keggGet("map01230")
modules_AA<-AA_map[[1]]$MODULE
AA_info<-list()
for(i in 1:length(names(modules_AA))){
  tryCatch({
    AA_info[[i]]<- keggGet(names(modules_AA)[i])
  }, error=function(e){
    message('An Error Occurred')
    print(e)
    return(NA)
  })
}
AA_modules_names<-c()
KOsall <-unique(unlist(strsplit(unlist(lapply(AA_info, function(x) names(x[[1]]$ORTHOLOGY))),",|[+]")))
KOsall_names <- unique(unlist(strsplit(unlist(lapply(AA_info, function(x) x[[1]]$ORTHOLOGY)),",|[+]")))

for(i in 1:length(KOsall)){
  my_list <- lapply(AA_info,function(x) grep(KOsall[i], names(x[[1]]$ORTHOLOGY)))
  
  non_empty_indices <- which(sapply(my_list, function(x) is.numeric(x) && length(x) > 0))
  print(non_empty_indices)
  if( length(non_empty_indices)<2){
    temp<-c(
      KOsall[i],
      AA_info[[non_empty_indices]][[1]]$ENTRY,
      AA_info[[non_empty_indices]][[1]]$NAME,
      AA_info[[non_empty_indices]][[1]]$CLASS,
      unlist(paste(unlist(names(AA_info[[non_empty_indices]][[1]]$PATHWAY)),collapse = ",")),
      unlist(paste(unlist(AA_info[[non_empty_indices]][[1]]$PATHWAY),collapse = ",")),
      length(unique(AA_info[[non_empty_indices]][[1]]$ORTHOLOGY))
    )
    #print(temp)
    AA_modules_names<-rbind(AA_modules_names,temp)
  }else{
    for (j in 1:length(non_empty_indices)){
      temp<-c(
        KOsall[i],
        AA_info[[non_empty_indices[j]]][[1]]$ENTRY,
        AA_info[[non_empty_indices[j]]][[1]]$NAME,
        AA_info[[non_empty_indices[j]]][[1]]$CLASS,
        unlist(paste(unlist(names(AA_info[[non_empty_indices[j]]][[1]]$PATHWAY)),collapse = ",")),
        unlist(paste(unlist(AA_info[[non_empty_indices[j]]][[1]]$PATHWAY),collapse = ",")),
        length(unique(AA_info[[non_empty_indices[j]]][[1]]$ORTHOLOGY))
      )
      #print(temp)
      AA_modules_names<-rbind(AA_modules_names,temp)
    }
   
  }

}
AA_modules_names<-as.data.frame(AA_modules_names)
colnames(AA_modules_names)<-c("KO_NAME","ENTRY","NAME","CLASS","PATHWAY-ID","PATHWAY","NR-ORTHOLOGY")
modules_AA_sel <- unique(AA_modules_names$ENTRY)
aa_data_modules <- c()
for(i in 1:length(modules_AA_sel)){
  kos<-AA_modules_names$KO_NAME[which(AA_modules_names$ENTRY==modules_AA_sel[i])]
  data <- cave_myco_ferrof
  data$AA_KEGG<- data$unique_CAVE_KEGG_acc%in%kos
  
  data_melted<-reshape2::melt(data)
  t <- unique(AA_modules_names$NAME[which(AA_modules_names$ENTRY==modules_AA_sel[i])])
  temp<-data_melted[data_melted$AA_KEGG == T,]
  if(nrow(temp)>0){
    temp$AA_module <- t
    aa_data_modules <- rbind(aa_data_modules,temp)
  }else{
    print(paste("module not found: ",    print(i), sep = "") )
  }


}
#ggplot(aa_data_modules, aes(x = value, y = aa_data_modules, group = variable)) + geom_point()




cave_myco_ferrof$unique_Ferro_cave <- cave_myco_ferrof$unique_CAVE_KEGG_acc%in%KOs_ferro_in_cave
cave_myco_ferrof$AA_Ferro<- cave_myco_ferrof$unique_CAVE_KEGG_acc%in%KO_FERRO_AA
cave_myco_ferrof$AA<- cave_myco_ferrof$unique_CAVE_KEGG_acc%in%KO_AA_ALL_clean
cave_myco_ferrof$AA_KEGG<- cave_myco_ferrof$unique_CAVE_KEGG_acc%in%KOsall

ggplot(cave_myco_ferrof, aes(x=CAVE_KOs_prev,y=Ferr_KOs_prev, color=AA))+geom_jitter()

library(beeswarm)
ggplot(cave_myco_ferrof, aes(x=CAVE_KOs_prev,y=Myco_KOs_prev, color=AA))+geom_jitter()


cave_myco_ferrof_melted<-reshape2::melt(cave_myco_ferrof)
ggplot(cave_myco_ferrof_melted[cave_myco_ferrof_melted$unique_Ferro_cave == T,], 
       aes(x = variable, y = value, group = unique_CAVE_KEGG_acc)) +
  geom_point(aes(color = variable), 
             position = position_jitter(width = 0.2, height = 0), 
             alpha = 0.6, size = 6) +  # Added jitter and transparency
  geom_line(color = "grey", linetype = "dashed", alpha = 0.2) +  # Dashed, transparent line
  theme_minimal() +  # Clean background
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Tilt x-axis labels
        text = element_text(size = 24), 
        legend.position = "none") +  # Remove legend
  scale_x_discrete(labels = c("Cave", expression(italic("Mycobacterium")), expression(italic("Ferroplasma")))) +  # Update x-axis labels
  labs(title = "", x = "Pangenomes", y = "Prevalence (%)") + # Title and labels
scale_color_manual(values = c('#fc8d62','#66c2a5','#8da0cb'))

  # K00373 K00371 K00374 K00370 K21908 K00038 high in mycobacterium
#check on which myco are present 
Myco_data_enriched_cave_Ferro<-Myco[Myco$KOfam_ACC%in%KOs_ferro_in_cave,]

mycoKOs_unique<-unique(Myco_data_enriched_cave_Ferro$KOfam_ACC)
for (i in mycoKOs_unique){
  tempdata_m<-Myco_data_enriched_cave_Ferro[which(Myco_data_enriched_cave_Ferro$KOfam_ACC==i),]
  tempdata_m_cave_index<-tempdata_m$genome_name%in%c("M_methanotrophicum","M_Bin2","M_Bin3")
  print(table(tempdata_m_cave_index))
}


KOs_cave_functions<-cave_myco_ferrof$unique_CAVE_KEGG_acc[cave_myco_ferrof$CAVE_KOs_prev >30 & cave_myco_ferrof$Myco_KOs_prev<10 & cave_myco_ferrof$Ferr_KOs_prev< 10]
KOs_cavehigh_lowmycoFerroplasma<- CAVE[CAVE$KOfam_ACC%in%KOs_cave_functions,]
CAVE[CAVE$KOfam_ACC%in%c("K17713", "K00754" , "K00148" , "K05878", "K05879","K05881","K06919","K08221"),]



ggplot(cave_myco_ferrof_melted[cave_myco_ferrof_melted$AA_KEGG==T,], aes(x=variable,y=value))+geom_boxplot()

data<-cave_myco_ferrof
rownames(data)<-data[,1]
library(apcluster)
## compute similarity matrix and run affinity propagation 
## (p defaults to median of similarity)
apres <- apcluster(negDistMat(r=2), data[,-1], details=TRUE)

## show details of clustering results
show(apres)

## plot clustering result
plot(apres, x)



library(plotly)
plot_ly(cave_myco_ferrof, x = ~CAVE_KOs_prev, y = ~Myco_KOs_prev, z = ~Ferr_KOs_prev) %>%
  add_markers(size =0.5) %>%
  layout(scene = list(
    xaxis = list(title = 'Cave'),
    yaxis = list(title = 'Myco'),
    zaxis = list(title = 'Ferro')
  ))

##################cog
############################



unique_CAVE_COG <- CAVE %>%
  select( COG20_FUNCTION_ACC, COG20_FUNCTION, COG20_CATEGORY_ACC,COG20_CATEGORY   ) %>%
  distinct()

unique_CAVE_COG<-unique_CAVE_COG[-which(unique_CAVE_COG$COG20_FUNCTION_ACC==""),]
CAVE_KOs_prev<-c()
for(i in unique_CAVE_COG$COG20_FUNCTION_ACC){
  prev<-length(unique(CAVE$genome_name[which(CAVE$COG20_FUNCTION_ACC==i)]))
  CAVE_KOs_prev<-c(CAVE_KOs_prev,prev)
}
CAVE_KOs_prev<-CAVE_KOs_prev/21*100
cave_cog_prev<-cbind(unique_CAVE_COG,CAVE_KOs_prev)

unique_Myco_COG <- Myco %>%
  select( COG20_FUNCTION_ACC, COG20_FUNCTION, COG20_CATEGORY_ACC,COG20_CATEGORY   ) %>%
  distinct()
unique_Myco_COG<-unique_Myco_COG[-which(unique_Myco_COG$COG20_FUNCTION_ACC==""),]
Myco_KOs_prev<-c()
for(i in unique_Myco_COG$COG20_FUNCTION_ACC){
  prev<-length(unique(Myco$genome_name[which(Myco$COG20_FUNCTION_ACC==i)]))
  Myco_KOs_prev<-c(Myco_KOs_prev,prev)
}
Myco_KOs_prev<-Myco_KOs_prev/69*100
Myco_cog_prev<-cbind(unique_Myco_COG,Myco_KOs_prev)

cave_myco_cog<-merge(cave_cog_prev,Myco_cog_prev,by.x="COG20_FUNCTION_ACC",by.y="COG20_FUNCTION_ACC", all=T)

unique_Ferroplasma_cog <- Ferroplasma %>%
  select( COG20_FUNCTION_ACC, COG20_FUNCTION, COG20_CATEGORY_ACC,COG20_CATEGORY   )%>%
  distinct()
unique_Ferroplasma_cog<-unique_Ferroplasma_cog[-which(unique_Ferroplasma_cog$COG20_FUNCTION_ACC==""),]
Ferr_KOs_prev<-c()
for(i in unique_Ferroplasma_cog$COG20_FUNCTION_ACC){
  prev<-length(unique(Ferroplasma$genome_name[which(Ferroplasma$COG20_FUNCTION_ACC==i)]))
  Ferr_KOs_prev<-c(Ferr_KOs_prev,prev)
}
Ferr_KOs_prev<-Ferr_KOs_prev/13*100
Ferr_cog_prev<-cbind(unique_Ferroplasma_cog,Ferr_KOs_prev)
cave_myco_ferro_cog<-merge(cave_myco_cog,Ferr_cog_prev,by.x="COG20_FUNCTION_ACC",by.y="COG20_FUNCTION_ACC", all=T)

cave_myco_ferro_cog_f<-cave_myco_ferro_cog %>%
  mutate(COG20_FUNCTION = coalesce(COG20_FUNCTION.x, COG20_FUNCTION.y, COG20_FUNCTION))%>%
  mutate(COG20_CATEGORY_ACC = coalesce(COG20_CATEGORY_ACC.x, COG20_CATEGORY_ACC.y, COG20_CATEGORY_ACC))%>%
  mutate(COG20_CATEGORY = coalesce(COG20_CATEGORY.x, COG20_CATEGORY.y, COG20_CATEGORY))%>%
  distinct(COG20_FUNCTION_ACC,COG20_FUNCTION,COG20_CATEGORY_ACC,COG20_CATEGORY,CAVE_KOs_prev,Myco_KOs_prev,Ferr_KOs_prev, .keep_all = FALSE) %>%
  select(COG20_FUNCTION_ACC,COG20_FUNCTION,COG20_CATEGORY_ACC,COG20_CATEGORY,CAVE_KOs_prev,Myco_KOs_prev,Ferr_KOs_prev)
cave_myco_ferro_cog_f <- cave_myco_ferro_cog_f %>%
  mutate(across(everything(), ~replace_na(., 0)))

plot_ly(cave_myco_ferro_cog_f, x = ~CAVE_KOs_prev, y = ~Myco_KOs_prev, z = ~Ferr_KOs_prev) %>%
  add_markers() %>%
  layout(scene = list(
    xaxis = list(title = 'Cave'),
    yaxis = list(title = 'Myco'),
    zaxis = list(title = 'Ferro')
  ))
