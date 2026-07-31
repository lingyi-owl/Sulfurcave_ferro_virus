MAG6_pan2_only_prot_crispr<-MAG6_pan2_only_prot[grep("CRISPR",MAG6_pan2_only_prot$KOfam.x),]

x1<-gsub("SC_MAG_00006_","",MAG6_pan2_only_prot_crispr$Majority_protein_IDs_u)


table(substr(x1,1,14))

# c_000000015663 c_000000049168 c_000000101225 c_000000170665 c_000000273288 
# 1              7              4              4              8 
SC_MAG_00006.gene_calls <- read.delim2("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/metagenomics/bin_by_bin/SC_MAG_00006/SC_MAG_00006-gene_calls.txt")

SC_MAG_00006.gene_calls_crispr<-SC_MAG_00006.gene_calls[SC_MAG_00006.gene_calls$contig%in%substr(x1,1,14),]
library(ggplot2)
library(gggenes)
SC_MAG_00006.gene_calls_crispr$groups<-"Unknown"
SC_MAG_00006.gene_calls_crispr$groups[which(SC_MAG_00006.gene_calls_crispr$COG20_FUNCTION!="")]<-"Others" 
SC_MAG_00006.gene_calls_crispr$groups[grepl("CRISPR",SC_MAG_00006.gene_calls_crispr$COG20_FUNCTION)]<-"CRISPR" 
SC_MAG_00006.gene_calls_crispr$CAS_types<-plyr::mapvalues(SC_MAG_00006.gene_calls_crispr$contig,unique(SC_MAG_00006.gene_calls_crispr$contig),c("~CAS-TypeID","CAS-TypeIB","CAS-TypeID+B","CAS-TypeIU","CAS-TypeIIIA"))
#write.csv(SC_MAG_00006.gene_calls_crispr, file = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/SC_MAG_00006_all_info.csv" )
SC_MAG_00006.gene_calls_crispr<-read.csv("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/SC_MAG_00006_all_info.csv")
SC_MAG_00006.gene_calls_crispr$orientation<-ifelse(SC_MAG_00006.gene_calls_crispr$direction=="f",1,0)

mag6_crispr<-ggplot(SC_MAG_00006.gene_calls_crispr, aes(xmin = start, xmax = stop, y = contig,fill=groups,forward = orientation)) +
  geom_gene_arrow() +
  facet_wrap(~ CAS_types, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3")+ggtitle("Ferroplasma MAG6") +
  geom_text_repel(aes(x = stop - ((stop-start)/2), y = 1, label = cog_simple, angle=45), nudge_x = .15,
                  box.padding = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 1,
                  segment.angle = 10)+theme_void()
mag6_crispr
library(Biostrings)
seq<-DNAStringSet(SC_MAG_00006.gene_calls_crispr$dna_sequence)
names(seq)<-paste(SC_MAG_00006.gene_calls_crispr$contig,SC_MAG_00006.gene_calls_crispr$gene_callers_id,sep = "_")
writeXStringSet(seq,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/Ferro_mag6_crispr.ffn")
###########
# ferroplassam circular
#########
Ferr_pan_only_prot_crispr<-Ferr_pan_only_prot[grep("CRISPR",Ferr_pan_only_prot$KOfam.x),]
x1<-gsub("SC_MAG_00006_","",Ferr_pan_only_prot_crispr$Majority_protein_IDs_u)

table(substr(x1,1,14))

PROKKA_04262022 <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/GENOMES/Ferroplasma_complete_genome/Ferroplasma_Prokka/PROKKA_04262022.gff", header=FALSE, comment.char="#")
PROKKA_04262022_crispr<-PROKKA_04262022[1:27,]
Ferr_pan_only_prot$COG20_CATEGORY.x[match(gsub("ID=","",unlist(lapply(strsplit(PROKKA_04262022_crispr$V9, ";"), function(x) x[1]))),Ferr_pan_only_prot$Majority_protein_IDs_u)]
PROKKA_04262022_crispr$cog<-Ferr_pan_only_prot$COG20_FUNCTION.x[match(gsub("ID=","",unlist(lapply(strsplit(PROKKA_04262022_crispr$V9, ";"), function(x) x[1]))),Ferr_pan_only_prot$Majority_protein_IDs_u)]
PROKKA_04262022_crispr$groups<-NA


# PROKKA_04262022_crispr$groups[grepl("CRISPR",PROKKA_04262022_crispr$cog)]<-PROKKA_04262022_crispr$cog[grepl("CRISPR",PROKKA_04262022_crispr$cog)]
PROKKA_04262022_crispr$groups[grepl("CRISPR",PROKKA_04262022_crispr$cog)]<-"CRISPR"
PROKKA_04262022_crispr$groups[grepl("Transposase",PROKKA_04262022_crispr$cog)]<-"Transposase"
PROKKA_04262022_crispr$groups[grepl("tRNA",PROKKA_04262022_crispr$cog)]<-"tRNA"
PROKKA_04262022_crispr$groups[which(PROKKA_04262022_crispr$V3=="repeat_region")]<-"repeat_region"
PROKKA_04262022_crispr$groups[which(PROKKA_04262022_crispr$cog=="")]<-"Unknown"
PROKKA_04262022_crispr$groups[is.na(PROKKA_04262022_crispr$groups)]<-"Others"
PROKKA_04262022_crispr$groups[16]<-"CRISPR" #BASED ON DATICATED TOOL
library(ggrepel)
#write.csv(PROKKA_04262022_crispr, file = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/FC_crispr_all_info.csv")
PROKKA_04262022_crispr<-read.csv("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/FC_crispr_all_info.csv")
PROKKA_04262022_crispr$orientation<-ifelse(PROKKA_04262022_crispr$V7=="+",1,0)
Fc_crispr<-ggplot(PROKKA_04262022_crispr, aes(xmin = V4, xmax = V5, y = V1,fill=groups,forward = orientation)) +
  geom_gene_arrow() +
  scale_fill_brewer(palette = "Set3")+ggtitle("Ferroplasma Circular - CAS-TypeID+B") +labs(y="")+
  geom_text_repel(aes(x = V5 - ((V5-V4)/2), y = 1, label = cog_simple, angle=45), nudge_x = .15,
                  box.padding = 0.5,
                  nudge_y = 0.2,
                  segment.curvature = -0.1,
                  segment.ncp = 1,
                  segment.angle = 10)+theme_void()
Fc_crispr
library(patchwork)
mag6_crispr / Fc_crispr



test<-readDNAStringSet("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/GENOMES/Ferroplasma_complete_genome/Ferroplasma_Prokka/PROKKA_04262022.fna")
xx<-subseq(test, start=1, end=28000) 
writeXStringSet(xx,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/CRISPR/Ferroc_crispr_conting.fasta")
