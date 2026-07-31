setwd("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git")

# =====================================================================
# Reproducible extract: construction of data frame `Myco_pan_only_prot`
# Source: plotting/Figure03_chrats.R (Sulfurcave_ferro_virus_git)
# Run from the repository root so the relative paths resolve.
#
# Dependency chain (original script line numbers in comments):
#   proteinGroups.xlsx -> data_proteomics -> data_proteomics1
#                                          -> data_proteomics_M_meth
#   MYCO_Pan_..._summary -> Myco_Pan (+ matchIDs via CD-HIT clusters)
#   Myco_Pan  +  data_proteomics_M_meth  --merge-->  Myco_pan_only_prot
# =====================================================================

library(dplyr)
library(tidyr)   # separate / fill / extract / pivot_wider used by read.cdhit.clstr

## --- 1. Proteomics table + LFQ intensities  (orig L5-15) -------------
ProteinGroups <- readxl::read_xlsx("proteomics/proteinGroups.xlsx", sheet = 2)

data_proteomics <- as.data.frame(ProteinGroups[, 1:2])
data_proteomics$LFQ_intensity_Cave_Biofilm_1        <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_1`)
data_proteomics$LFQ_intensity_Cave_Biofilm_2        <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_2`)
data_proteomics$LFQ_intensity_Culture               <- as.numeric(ProteinGroups$`LFQ intensity invitro`)
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_1  <- as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_1`))
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2  <- as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_2`))
data_proteomics$log10_LFQ_intensity_Culture         <- as.numeric(log10(ProteinGroups$`LFQ intensity invitro`))
data_proteomics[data_proteomics == "-Inf"] <- 0

## --- 2. Keep proteins detected in both cave-biofilm replicates (orig L83)
data_proteomics1 <- data_proteomics[which(data_proteomics$LFQ_intensity_Cave_Biofilm_1 > 0 &
                                            data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 > 0), ]

## --- 3. Subset to M. methanotrophicum proteins (tag KDJLIKBO) (orig L194)
data_proteomics_M_meth <- data_proteomics1[grepl("KDJLIKBO", data_proteomics1$`Majority protein IDs`), ]

## --- 4. CD-HIT .clstr reader  (orig L201-219) -----------------------
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

## --- 5. Mycobacterium pan-genome gene-cluster table  (orig L231-233)
Myco_Pan_gene_clusters_summary.txt <- read.delim("pangenomics/MYCO_Pan_gene_clusters_summary.txt.gz")
Myco_Pan <- Myco_Pan_gene_clusters_summary.txt[which(Myco_Pan_gene_clusters_summary.txt$genome_name == "M_methanotrophicum"), ]
Myco_Pan <- Myco_Pan[order(Myco_Pan$gene_callers_id, decreasing = F), ]

## --- 6. Map pan-genome gene_callers_id -> proteomics IDs via CD-HIT (orig L369-378)
temp_Myc <- read.cdhit.clstr("proteomics/Mycobacterium_connect_100.clstr")

matchIDs <- rep(NA, nrow(Myco_Pan))
for (i in unique(temp_Myc$Cluster)) {
  temp1 <- temp_Myc[which(i == temp_Myc$Cluster), ]
  if (nrow(temp1) == 2) {
    idx <- which(Myco_Pan$gene_callers_id == as.integer(temp1$Seq.Name[1]))
    matchIDs[idx] <- temp1$Seq.Name[2]
  }
}
Myco_Pan$matchIDs <- matchIDs

## --- 7. Extract the KDJLIKBO majority-protein ID key (orig L383) -----
data_proteomics_M_meth$Majority_protein_IDs_u <-
  unlist(lapply(strsplit(data_proteomics_M_meth$`Majority protein IDs`, ";"),
                function(x) substr(x[grep("KDJLIKBO", x)], 1, 14)[1]))

## --- 8. Merge pan-genome + proteomics  (orig L386) ------------------
Myco_pan_only_prot <- merge(Myco_Pan, data_proteomics_M_meth,
                            by.x = "matchIDs", by.y = "Majority_protein_IDs_u", all == T)

## --- 9. Derived columns  (orig L390, L391, L396) -------------------
Myco_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm <- rowMeans(Myco_pan_only_prot[, 32:33], na.rm = T)
Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm       <- rowMeans(Myco_pan_only_prot[, 29:30], na.rm = T)
Myco_pan_only_prot$gene_cluster_perc                <- Myco_pan_only_prot$num_genomes_gene_cluster_has_hits / 69 * 100

######################################################################
myco_defense_finder_systems <- read.delim2("metagenomics/defencefinder/Myco/defense_finder_systems.tsv")
pall<-unlist(strsplit(myco_defense_finder_systems$protein_in_syst,","))
Myco_pan_only_prot[Myco_pan_only_prot$matchIDs%in%pall,]

Ferr_c_defense_finder_systems <- read.delim2("metagenomics/defencefinder/Ferro/defense_finder_systems.tsv")
pall<-unlist(strsplit(Ferr_c_defense_finder_systems$protein_in_syst,","))
Ferr_pan_only_prot[Ferr_pan_only_prot$matchIDs%in%pall,]

MAG6_c_defense_finder_systems <- read.delim2("metagenomics/defencefinder/MAG6/defense_finder_systems.tsv")
pall<-unlist(strsplit(MAG6_c_defense_finder_systems$protein_in_syst,","))
MAG6_pan2_only_prot[gsub("SCMAG00006_","",MAG6_pan2_only_prot$matchIDs)%in%pall,]

library(ggplot2)
library(ggsignif)
library(reshape2)

def_table<-c()
for(i in 1:nrow(myco_defense_finder_systems)){
  p<-unlist(strsplit(myco_defense_finder_systems$protein_in_syst[i],","))
  if(any(Myco_pan_only_prot$matchIDs%in%p==T)){
    def_temp<-Myco_pan_only_prot[ Myco_pan_only_prot$matchIDs%in%p,]
    def_temp$system<-myco_defense_finder_systems$subtype[i]
  }else{
    def_temp<-  as.data.frame(matrix(NA, nrow = length(p), ncol = ncol(Myco_pan_only_prot)) )
    def_temp$system<-myco_defense_finder_systems$subtype[i]
  }
  
  def_table<-rbind(def_table,def_temp)
}
Myco_ribo<-Myco_pan_only_prot[grepl("ribosomal",Myco_pan_only_prot$KOfam),]
Myco_ribo$system<-"Ribosomal"
colnames(def_table)<-colnames(Myco_ribo)
data_myco_def<-rbind(def_table,Myco_ribo)

df_melted <- melt(data_myco_def, measure.vars = c("log10_LFQ_intensity_Cave_Biofilm_1", "log10_LFQ_intensity_Cave_Biofilm_2"))
df_melted$Sample<-gsub("log10_LFQ_intensity_Cave_Biofilm_","Sample ",df_melted$variable)
#['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c']
myco<-ggplot(df_melted, aes(x = system, y = value)) +
  geom_jitter(aes( shape =Sample), size=6, width = 0.2, alpha = 0.5, color="#e31a1c") +ylim(6.3,10.5)+
  geom_boxplot(aes(fill = system), alpha = 0.6, outlier.shape = NA) +
  theme_minimal() +
  labs(title = "M. methanotrophicum",
       x = "Functional category",
       y = "LFQ intensity (log10) Cave Biofilm",
       color = "Functional category",
       fill = "Functional category")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(text = element_text(size = 18))

