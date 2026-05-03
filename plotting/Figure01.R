# SAMPLES_MERGED.COVs <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/contings-coverages/SAMPLES_MERGED-COVs.txt")
# unlist(lapply(strsplit(CAT$genus,":"), function(x) x[1]))
# CAT <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/contings-coverages/CAT_run_20230309.contig2classification.names.txt")
# colnames(CAT)[1]<- "contig"  
# CAT$genus_clean<-unlist(lapply(strsplit(CAT$genus,":"), function(x) x[1]))
# CAT$class_clean<-unlist(lapply(strsplit(CAT$class,":"), function(x) x[1]))
# CAT$superkingdom_clean<-unlist(lapply(strsplit(CAT$superkingdom,":"), function(x) x[1]))
# CAT$phylum<-unlist(lapply(strsplit(CAT$phylum,":"), function(x) x[1]))
# 
# contigs<-Biostrings::readDNAStringSet("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/contings-coverages/SAMPLES_MERGED-CONTIGS.fa")
# dataframe_contigs_width<-data.frame(names(contigs))
# colnames(dataframe_contigs_width)[1]<-"contig"
# dataframe_contigs_width$width<-contigs@ranges@width
# SAMPLES_MERGED.COVs_cat1<-merge(SAMPLES_MERGED.COVs,dataframe_contigs_width,by="contig", all=T)
# SAMPLES_MERGED.COVs_cat2<-merge(SAMPLES_MERGED.COVs_cat1,CAT,by="contig", all=T)
# library(tidyverse)
# class<-SAMPLES_MERGED.COVs_cat2 %>% group_by(class_clean) %>% summarise(Sample_ERR10036468_m = sum(Sample_ERR10036468)
#                                                                        ,Sample_ERR10036469_m = sum(Sample_ERR10036469),
#                                                                        Sample_ERR10036470_m = sum(Sample_ERR10036470))
# superkingdom_clean<-SAMPLES_MERGED.COVs_cat2 %>% group_by(superkingdom_clean) %>% summarise(Sample_ERR10036468_m = sum(Sample_ERR10036468/width)
#                                                                        ,Sample_ERR10036469_m = sum(Sample_ERR10036469/width),
#                                                                        Sample_ERR10036470_m = sum(Sample_ERR10036470/width))

library(rlang)
library(ggplot2)
setwd("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/KAIJU_reads")
kaiju_reads<-list.files(pattern = "kaiju_ERR")

kaiju_phylum<-c()
for(i in kaiju_reads[grep("phylum",kaiju_reads)]){
  temp<-read.delim(i)
  kaiju_phylum<-rbind(kaiju_phylum,temp)
}
kaiju_phylum_fac<-factor(kaiju_phylum$file)
levels(kaiju_phylum_fac)[levels(kaiju_phylum_fac) == "ERR10036468.out"] <- "Cave biofilm 1"
levels(kaiju_phylum_fac)[levels(kaiju_phylum_fac) == "ERR10036469.out"] <- "Cave biofilm 2"
levels(kaiju_phylum_fac)[levels(kaiju_phylum_fac) == "ERR10036470.out"] <- "Lab: CH4"
kaiju_phylum$Sample_name<-as.character(kaiju_phylum_fac)
temp<-kaiju_phylum$taxon_name
temp[which(kaiju_phylum$percent<1)]<-"Below 1%"
kaiju_phylum$taxon_final<-temp
kaiju_phylum$taxon_final[kaiju_phylum$taxon_final=="cannot be assigned to a (non-viral) phylum"]<-"unclassified"
kaiju_phylum$taxon_final[kaiju_phylum$taxon_final=="Candidatus Thermoplasmatota"]<-"Thermoplasmatota"
colors_p<-c('#66c2a5','#fc8d62','#e78ac3','#a6d854','#8da0cb','#ffd92f')
kaiju_phylum_ggplot<-ggplot(kaiju_phylum, aes(fill=taxon_final, y=percent, x=Sample_name)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=colors_p)+
  labs(x="",y="Relative abundance (%)",fill= "Taxonomy (phylum)",subtitle = "", title="Phylum")+
  theme(text = element_text(size = 24)) + theme(axis.text.x = element_text(angle = 45,vjust=0.5))

kaiju_phylum_ggplot
kaiju_genus<-c()
for(i in kaiju_reads[grep("genus",kaiju_reads)]){
  temp<-read.delim(i)
  kaiju_genus<-rbind(kaiju_genus,temp)
}
kaiju_genus_fac<-factor(kaiju_genus$file)
levels(kaiju_genus_fac)[levels(kaiju_genus_fac) == "ERR10036468.out"] <- "Cave biofilm 1"
levels(kaiju_genus_fac)[levels(kaiju_genus_fac) == "ERR10036469.out"] <- "Cave biofilm 2"
levels(kaiju_genus_fac)[levels(kaiju_genus_fac) == "ERR10036470.out"] <- "Lab: CH4"
kaiju_genus$Sample_name<-as.character(kaiju_genus_fac)
temp<-kaiju_genus$taxon_name
temp[which(kaiju_genus$percent<1)]<-"Below 1%"
kaiju_genus$taxon_final<-temp
kaiju_genus$taxon_final[kaiju_genus$taxon_final=="cannot be assigned to a (non-viral) genus"]<-"unclassified"
colors_g<-c('#377eb8','#a6d854','#fc8d62','#8da0cb','#66c2a5','#ffd92f')
kaiju_genus_ggplot<-ggplot(kaiju_genus, aes(fill=taxon_final, y=percent, x=Sample_name)) + 
  geom_bar(position="stack", stat="identity")+scale_fill_manual(values=colors_g)+
  labs(x="",y="Relative abundance (%)",fill= "Taxonomy (genus)",subtitle = "", title="Genus")+
  theme(text = element_text(size = 24))  + theme(axis.text.x = element_text(angle = 45,vjust=0.5))
kaiju_genus_ggplot

# #3
# library(vegan)
# test<-transformation(kaiju_genus)
# ##################
# 
# # Calculate row means
# row_means <- rowMeans(test,na.rm = T)
# 
# # Filter rows where the average count is below 10
# filtered_count_matrix <- test[row_means >= 10, ]
# filtered_count_matrix[is.na(filtered_count_matrix)] <- 0
# 
# richness <- specnumber(t(filtered_count_matrix))
# shannon <- diversity(t(filtered_count_matrix), index = "shannon")


#####################

####
mean_coverage <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/SAMPLES_SUMMARY_SC_MAGS/bins_across_samples/mean_coverage.txt")
mean_coverage_log10<-apply(mean_coverage[,-1], 2, log10)
colnames(mean_coverage_log10)<-paste("log10_",colnames(mean_coverage_log10),sep = "")
mean_coverage<-cbind(mean_coverage,mean_coverage_log10)

abundance <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/SAMPLES_SUMMARY_SC_MAGS/bins_across_samples/abundance.txt")
relab<-apply(abundance[,-1], 2, function(x) x/sum(x))
colnames(relab)<-paste("relab_",colnames(relab),sep = "")
abundance<- cbind(abundance,relab)
bins_percent_recruitment <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/SAMPLES_SUMMARY_SC_MAGS/bins_across_samples/bins_percent_recruitment.txt")
detection <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/SAMPLES_SUMMARY_SC_MAGS/bins_across_samples/detection.txt")

bin_sum <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/SAMPLES_SUMMARY_SC_MAGS/bins_summary.txt")
library(ggrepel)
mean_cov_all<-merge(mean_coverage,bin_sum,by="bins" )
theme_set(theme_bw())
mean_cov_all$MAG_quality<-ifelse(grepl("_MAG_",mean_cov_all$bins)==T,"High","Low")

mean_cov_all$phylum<-factor(mean_cov_all$t_phylum)
mean_cov_all$phylum[mean_cov_all$phylum==""]<-"unclassified"
mean_cov_all$genus<-ifelse(mean_cov_all$log10_Sample_ERR10036468>2|mean_cov_all$log10_Sample_ERR10036469>2,mean_cov_all$t_genus,"" )
#write.csv(mean_cov_all, "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/DNA_mean_cov_all.csv")

colors_pB<-c('#66c2a5','#e5c494','#e78ac3',"#e31a1c",'#a6d854','#8da0cb','#ffd92f')
bins1<-ggplot(mean_cov_all,aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469,  color=log10_Sample_ERR10036470,shape=MAG_quality))+geom_point(size=5)+
  scale_color_viridis_c()+scale_fill_manual(values=colors_pB)+labs(x="Cave biofilm 1",y="Cave biofilm 2",caption= "Axes and color scale represent the mean DNA coverage (log10).",color= "Lab: CH4")+
  theme(text = element_text(size = 12)) +
   guides(
    fill = guide_legend(
      override.aes = aes(label = "")
    )
  ) + geom_segment(aes(x = 0, y = 0, xend = 4, yend = 4), linetype="dashed", color = "red")
bins2<-ggplot(mean_cov_all,aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469,  color=log10_Sample_ERR10036470,shape=MAG_quality))+geom_point(size=5)+
  scale_color_viridis_c()+scale_fill_manual(values=colors_pB)+labs(x="Cave biofilm 1",y="Cave biofilm 2",caption= "Axes and color scale represent the mean DNA coverage (log10).",color= "Lab: CH4")+
  theme(text = element_text(size = 12)) +
  geom_label_repel(box.padding = 0.5, max.overlaps = Inf,aes( label = genus, fill=phylum),color="black",size=3) + guides(
    fill = guide_legend(
      override.aes = aes(label = "")
    )
  )+ geom_segment(aes(x = 0, y = 0, xend = 4, yend = 4), linetype="dashed", color = "red")
mean_cov_all$highlight<-""
mean_cov_all$highlight[match(c("SC_MAG_00008","SC_MAG_00016","SC_MAG_00006"),mean_cov_all$bins)]<-c("Ferroplasma MAG8","C. M. methanotrophicum","Ferroplasma MAG6")
bins3<-ggplot(mean_cov_all,aes(log10_Sample_ERR10036468, log10_Sample_ERR10036469))+geom_point(aes(color=log10_Sample_ERR10036470,shape=MAG_quality), size=8)+
  scale_color_viridis_c()+labs(x="Cave biofilm 1 - mean DNA coverage (log10)",y="Cave biofilm 2 - mean DNA coverage (log10)",caption= "",color= "Lab: CH4", shape="MAG quality")+
  theme(text = element_text(size = 24)) +
  geom_label_repel(box.padding = 0.5,aes(label = highlight),color="black",size=10, fill = alpha(c("white"),0.5),ylim=c( NA,1),arrow = arrow(
    length = unit(0.02, "npc"), type = "closed", ends = "first"
  )) + guides(
    fill = guide_legend(
      override.aes = aes(label = "")
    )
  )+ geom_segment(aes(x = 0, y = 0, xend = 4, yend = 4), linetype="dashed", color = "red") 
bins3
library(patchwork)
patchwork <- ((kaiju_phylum_ggplot + kaiju_genus_ggplot)+ plot_layout( guides = "collect"))  / bins1
patchwork + plot_annotation(tag_levels = 'a')

patchwork2 <- ((kaiju_phylum_ggplot + kaiju_genus_ggplot)+ plot_layout( guides = "collect"))  / bins2
patchwork2 + plot_annotation(tag_levels = 'a')

patchwork <- ((kaiju_phylum_ggplot + kaiju_genus_ggplot)+ plot_layout( guides = "collect"))  + bins3
patchwork + plot_annotation(tag_levels = 'a')
ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/figure1.pdf", width =24, height = 10)
