#Figure 3 combine metaproteomics



#load pan genome Ferr and myco

Ferroplasma_Pan_2_gene_clusters_summary.txt <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/Ferroplasma/Ferroplasma_Pan_2/SUMMARY_few_classes/Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
Ferr<-Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name=="SFerroplasmacircular"),]
Ferr<-Ferr[order(Ferr$gene_callers_id,decreasing = F),]
MAG6<-Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name=="SMAG00006"),] 
MAG6<-MAG6[order(MAG6$gene_callers_id,decreasing = F),]

Myco_Pan_gene_clusters_summary.txt <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/mycobacterium/MYCO-SUMMARY/MYCO_Pan_gene_clusters_summary.txt.gz")
Myco_Pan<-Myco_Pan_gene_clusters_summary.txt[which(Myco_Pan_gene_clusters_summary.txt$genome_name=="M_methanotrophicum"),]
Myco_Pan<-Myco_Pan[order(Myco_Pan$gene_callers_id,decreasing = F),]

#####eggnog
sulfur_cave_ref_proteins_eggnog.emapper <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/eggnog/sulfur_cave_ref_proteins_eggnog.emapper.annotations", header=FALSE, comment.char="#")

library(KEGGREST)
kegg_pws<-unique(unlist(strsplit(sulfur_cave_ref_proteins_eggnog.emapper$V13,",")))
kegg_pws<-kegg_pws[grep("ko",kegg_pws)]
pwinfo<-list()
for(i in 1:length(kegg_pws)){
  tryCatch({
    pwinfo[[i]]<-keggGet(kegg_pws[i])
  }, error=function(e){})
}

nr_KOs<-unlist(lapply(pwinfo, function(x) length(x[[1]]$ORTHOLOGY)))

pws_names<-c()
for(i in 1:length(pwinfo)){
  temp<-c(unlist(pwinfo[[i]][[1]]$ENTRY),unlist(pwinfo[[i]][[1]]$NAME[1]),
          unlist(length(pwinfo[[i]][[1]]$ORTHOLOGY)),
          unlist(paste(unlist(names(pwinfo[[i]][[1]]$ORTHOLOGY)),collapse = ",")),
          unlist(pwinfo[[i]][[1]]$CLASS)
          )
  #print(temp)
  pws_names<-rbind(pws_names,temp)
}
pws_names<-as.data.frame(pws_names)

AA_mod<-names(pwinfo[[18]][[1]]$MODULE)
modinfo<-list()
for(i in 1:length(AA_mod)){
  tryCatch({
    modinfo[[i]]<-keggGet(AA_mod[i])
  }, error=function(e){})
}
AAmod_names<-c()
for(i in 1:length(modinfo)){
  temp<-c(unlist(modinfo[[i]][[1]]$ENTRY),unlist(modinfo[[i]][[1]]$NAME[1]),
          unlist(length(modinfo[[i]][[1]]$ORTHOLOGY)),
          unlist(paste(unlist(names(modinfo[[i]][[1]]$ORTHOLOGY)),collapse = ",")),
          unlist(modinfo[[i]][[1]]$CLASS) )
  #print(temp)
  AAmod_names<-rbind(AAmod_names,temp)
}
AAmod_names<-as.data.frame(AAmod_names)
pws_mod_names<-rbind(pws_names,AAmod_names)
pws_mod_names<-pws_mod_names[-which(pws_mod_names$V3==0),]
#write.csv(as.data.frame(pws_mod_names),"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/eggnog/pw_info.csv")
pws_mod_names <- read.csv("~/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/eggnog/pw_info.csv")


################ combine proteomics

ProteinGroups<-readxl::read_xlsx("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/proteinGroups.xlsx",sheet = 2)
ids<-readxl::read_xlsx("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/proteinGroups.xlsx",sheet = 1)
data_proteomics<-as.data.frame(ProteinGroups[,1:2])
data_proteomics$LFQ_intensity_Cave_Biofilm_1 <-as.numeric(ProteinGroups$`LFQ intensity cave_wcl_1`)
data_proteomics$LFQ_intensity_Cave_Biofilm_2 <-as.numeric(ProteinGroups$`LFQ intensity cave_wcl_2`)
data_proteomics$LFQ_intensity_Culture <-as.numeric(ProteinGroups$`LFQ intensity invitro`)
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_1 <-as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_1`))
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 <-as.numeric(log10(ProteinGroups$`LFQ intensity cave_wcl_2`))
data_proteomics$log10_LFQ_intensity_Culture <-as.numeric(log10(ProteinGroups$`LFQ intensity invitro`))
data_proteomics[data_proteomics=="-Inf"]<-0

data_proteomics_ferro<-data_proteomics[grepl("KNPMNEEE",data_proteomics$`Majority protein IDs`),]
data_proteomics_mag6<-data_proteomics[grepl("SCMAG00006",data_proteomics$`Majority protein IDs`),]
data_proteomics_M_meth<-data_proteomics[grepl("KDJLIKBO",data_proteomics$`Majority protein IDs`),]


################ sum lfg proteomics per pw
pws_mod_names_r<-pws_mod_names[-grep("Human|Organismal Systems",pws_mod_names$V5),]
pws_mod_names_r<-pws_mod_names_r[-grep(" - fly| - worm",pws_mod_names_r$V2),]

Ferr_data<-pws_mod_names_r
# Ferr_pan_only_prot_reduced<-Ferr_pan_only_prot[apply(Ferr_pan_only_prot[,32:33], 1, function(x) ( all(x>0))),]
# Ferr_pan_only_prot_reduced$gene_cluster_perc<-Ferr_pan_only_prot_reduced$num_genomes_gene_cluster_has_hits/13*100
data_temp_proteins<-Ferr_pan_only_prot # from other file
match<-apply(pws_mod_names_r, 1, function(x) sum(data_temp_proteins$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")) ,data_temp_proteins$KOfam_ACC) ] ,na.rm = T))
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(data_temp_proteins$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")) ,data_temp_proteins$KOfam_ACC)]>0)))
Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm<-match
Ferr_data$Ferr_NR_active<-NR_active
Ferr_data$org<-rep("Ferroplasma c.",nrow(Ferr_data))

Ferr_data$Ferr_log10_Sum_LFQ_intensity_Cave_Biofilm<-log10(Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm)
Ferr_data<- Ferr_data[-which(Ferr_data$Ferr_Sum_LFQ_intensity_Cave_Biofilm==0),]
Ferr_data<- Ferr_data[order(Ferr_data$Ferr_log10_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),]
Ferr_data$pw_nr_active_perc<-Ferr_data$Ferr_NR_active/as.integer(Ferr_data$V3)*100

sub_Ferr_data<-Ferr_data[1:10,]
sub_Ferr_data$V2<-factor(sub_Ferr_data$V2,levels = sub_Ferr_data$V2)
gg1<-ggplot(sub_Ferr_data,aes(Ferr_log10_Sum_LFQ_intensity_Cave_Biofilm,V2,color=Ferr_NR_active))+geom_point(size=3)+
  labs(y="",x="",title = "Ferroplasma c.",color="N. proteins")+
  scale_colour_continuous(type = "viridis")
gg1

###########################
MAG6_data<-pws_mod_names_r
data_temp_proteins<-MAG6_pan2_only_prot
match<-apply(pws_mod_names_r, 1, function(x) sum(data_temp_proteins$LFQ_intensity_Cave_Biofilm[match( unlist(strsplit(x[5],",")),data_temp_proteins$KOfam_ACC ) ] ,na.rm = T))
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(data_temp_proteins$LFQ_intensity_Cave_Biofilm[match( unlist(strsplit(x[5],",")),data_temp_proteins$KOfam_ACC ) ]>0)))
MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm<-match
MAG6_data$MAG6_NR_active<-NR_active
MAG6_data$org<-rep("Ferroplasma MAG6",nrow(MAG6_data))

MAG6_data$MAG6_log10_Sum_LFQ_intensity_Cave_Biofilm<-log10(MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm)
MAG6_data<- MAG6_data[-which(MAG6_data$MAG6_Sum_LFQ_intensity_Cave_Biofilm==0),]
MAG6_data<- MAG6_data[order(MAG6_data$MAG6_log10_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),]
MAG6_data$pw_nr_active_perc<-MAG6_data$MAG6_NR_active/as.integer(MAG6_data$V3)*100

sub_MAG6_data<-MAG6_data[1:10,]
sub_MAG6_data$V2<-factor(sub_MAG6_data$V2,levels = sub_MAG6_data$V2)
gg2<-ggplot(sub_MAG6_data,aes(MAG6_log10_Sum_LFQ_intensity_Cave_Biofilm,V2,color=MAG6_NR_active))+geom_point(size=3)+
  labs(y="",x="",title = "Ferroplasma MAG6",color="N. proteins")+
  scale_colour_continuous(type = "viridis")
gg2


#combine ferroplasma 
colnames(Ferr_data)<-c( "X","Pathway","Pathway_name"
                       , "Pathway_total_kos","KOs","BRITE"
                       , "Sum_LFQ_intensity_Cave_Biofilm",   "NR_active","org"                                    
                       ,"log10_Sum_LFQ_intensity_Cave_Biofilm" ,"pw_nr_active_perc" )
colnames(MAG6_data)<-colnames(Ferr_data)
data_ferro_c_MAG6<-rbind(Ferr_data, MAG6_data)
selection<-unique(c(MAG6_data[1:10,3],Ferr_data[1:10,3]))

data_ferro_c_MAG6_sel<-data_ferro_c_MAG6[data_ferro_c_MAG6$Pathway_name%in%selection,]
data_ferro_c_MAG6_sel$Pathway_name[order(data_ferro_c_MAG6_sel$Sum_LFQ_intensity_Cave_Biofilm)]
data_ferro_c_MAG6_sel$Pathway_name<-factor(data_ferro_c_MAG6_sel$Pathway_name,levels = unique(data_ferro_c_MAG6_sel$Pathway_name[order(data_ferro_c_MAG6_sel$Sum_LFQ_intensity_Cave_Biofilm)]))
gg_ferr2<-ggplot(data_ferro_c_MAG6_sel,aes(log10_Sum_LFQ_intensity_Cave_Biofilm,Pathway_name,color=NR_active, shape=org))+geom_point(size=10)+
  labs(y="Pathways (top10)",x="Sum protein LFQ intensity (log10) - Cave biofilm",title = "",color="N. proteins",shape="Organism")+
  scale_colour_continuous(type = "viridis")+theme(text = element_text(size = 24))
gg_ferr2
gg_poster<-gb+theme(legend.direction="horizontal",legend.position="top") + prot_dna+theme(legend.position="none") 
gg_poster 

gg_ferr2

  gg3+theme(legend.direction="horizontal",legend.position="top")+gg1+theme(legend.direction="horizontal",legend.position="top")+gg2+theme(legend.direction="horizontal",legend.position="top")
wrap_elements(gg) +
  labs(tag = "Sum protein LFQ intensity (log10) - Cave biofilm") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )
##########################
Myco_data<-pws_mod_names_r
match<-apply(pws_mod_names_r, 1, function(x) sum(Myco_pan_only_prot_reduced$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[4],",")),Myco_pan_only_prot_reduced$KOfam_ACC ) ] ,na.rm = T))
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(Myco_pan_only_prot_reduced$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[4],",")),Myco_pan_only_prot_reduced$KOfam_ACC )  ]>0)))
Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm<-match
Myco_data$Myco_NR_active<-NR_active
Myco_data$org<-rep("M.meth",nrow(Myco_data))

Myco_data$Myco_log10_Sum_LFQ_intensity_Cave_Biofilm<-log10(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm)
Myco_data<- Myco_data[-which(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm==0),]
Myco_data<- Myco_data[order(Myco_data$Myco_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),]
Myco_data$pw_nr_active_perc<-Myco_data$Myco_NR_active/as.integer(Myco_data$V3)*100

sub_Myco_data<-Myco_data[1:10,]
sub_Myco_data$V2<-factor(sub_Myco_data$V2,levels = sub_Myco_data$V2)
gg3<-ggplot(sub_Myco_data,aes(Myco_log10_Sum_LFQ_intensity_Cave_Biofilm,V2,color=Myco_NR_active))+geom_point(size=3)+
  labs(y="Pathways (top10)",x="",title = "M.methanotrophicum",color="N. proteins")+
  scale_colour_continuous(type = "viridis")
gg3

gg<-gg3+theme(legend.direction="horizontal",legend.position="top")+gg1+theme(legend.direction="horizontal",legend.position="top")+gg2+theme(legend.direction="horizontal",legend.position="top")
wrap_elements(gg) +
  labs(tag = "Sum protein LFQ intensity (log10) - Cave biofilm") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )
library(reshape2)
top_pw<-list(ferr=pws_mod_names_r[order(pws_mod_names_r$Ferr_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),2][1:10],
             mag6=pws_mod_names_r[order(pws_mod_names_r$MAG6_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),2][1:10],
             Myco=pws_mod_names_r[order(pws_mod_names_r$Myco_Sum_LFQ_intensity_Cave_Biofilm,decreasing = T),2][1:10])

unique(unlist(top_pw))

library(gplots)
v.table <- venn(top_pw)
v.table

which(pws_mod_names$Pathway=="ko00720")# 
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[71,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[71,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")


which(pws_mod_names$Pathway=="ko00710")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[79,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[79,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")


sum(MAG6_pan2_only_prot_reduced$LFQ_intensity_Cave_Biofilm[match( unlist(strsplit(pws_mod_names[70,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  )] ,na.rm = T)

#############
#Idividula
###########
#KNPMNEEE ferro
#mag6

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
temp_Fer<-read.cdhit.clstr("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/connect_proteomics/Ferroplasma_connect_100.clstr")

##########Ferr_pan2
matchIDs<-rep(NA,nrow(Ferr))
for(i in unique(temp_Fer$Cluster)){
  temp1<-temp_Fer[which(i==temp_Fer$Cluster),]
  if(nrow(temp1)==2){
    idx<-which(Ferr$gene_callers_id==as.integer(temp1$Seq.Name[1]))
    matchIDs[idx]<-temp1$Seq.Name[2]
  }
}
Ferr$matchIDs<-matchIDs
strsplit(data_proteomics_ferro$`Majority protein IDs`,";")
data_proteomics_ferro$Majority_protein_IDs_u<-unlist(lapply(strsplit(data_proteomics_ferro$`Majority protein IDs`,";"), function(x) substr(x[grep("KNPMNEEE",x)], 1, 14)[1]))
data_proteomics_ferro$Majority_protein_IDs_all<-unlist(lapply(strsplit(data_proteomics_ferro$`Majority protein IDs`,";"), function(x) paste(substr(x[grep("KNPMNEEE",x)], 1, 14),collapse = ";")))

Ferr_pan_only_prot<-merge(Ferr,data_proteomics_ferro,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)

ggplot(Ferr_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=bin_name))+geom_point()
Ferr_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm<- rowMeans(Ferr_pan_only_prot[,33:34],na.rm = T)
Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm<- rowMeans(Ferr_pan_only_prot[,30:31],na.rm = T)
Ferr_pan_only_prot_reduced<-Ferr_pan_only_prot[apply(Ferr_pan_only_prot[,33:34], 1, function(x) ( all(x>0))),]
library(ggbeeswarm)
ggplot(Ferr_pan_only_prot_reduced,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
Ferr_pan_only_prot_reduced$gene_cluster_perc<-Ferr_pan_only_prot_reduced$num_genomes_gene_cluster_has_hits/13*100

ggplot(Ferr_pan_only_prot_reduced,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 10)) + scale_y_continuous(limits = c(6, 10))
sor<-c("K16952")
Sulf<-c("K17218","K01739")
succ_succ_coac<-c("K01902","K01903")
oxoglutare<-c("K00174","K00175","K00176")
isocitrate<-c("K00031")
citrate<-c("K01681","K01682")
SMMO<-c("K16157","K16158","K16159","K16160","K16161","K16162")
CRISP<-c("K07725","K15342","K19091","K07464","K19120","K19121","K19122","K09951","K07012","K19116","K19075","K19115","K09002","K19090","K19136","K07016","K19114","K19138","K19140","K19135","K19132","K19131","K19139")
unique(Ferroplasma_Pan_2_gene_clusters_summary.txt$KOfam_ACC[grep("CRISP",Ferroplasma_Pan_2_gene_clusters_summary.txt$KOfam)])
rTCA<-c(succ_succ_coac,oxoglutare,isocitrate,citrate)
library(ggrepel)
pp1<-ggplot(Ferr_pan_only_prot_reduced,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
   geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "x/y axes shows the protein LFQ intensity (log10)",title="Ferroplasma c.")
# - Cave biofilm
pp1
##########MAG6
temp_mag6<-read.cdhit.clstr("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/connect_proteomics/MAG6_connect_100.clstr")
matchIDs<-rep(NA,nrow(MAG6))
for(i in unique(temp_mag6$Cluster)){
  temp1<-temp_mag6[which(i==temp_mag6$Cluster),]
  if(nrow(temp1)==2){
    idx<-which(MAG6$gene_callers_id==as.integer(temp1$Seq.Name[1]))
    matchIDs[idx]<-temp1$Seq.Name[2]
  }
}
MAG6$matchIDs<-matchIDs
strsplit(data_proteomics_mag6$`Majority protein IDs`,";")
temp_list<-lapply(strsplit(data_proteomics_mag6$`Majority protein IDs`,";"), function(x) x[grep("SC_MAG_00006",x)])
data_proteomics_mag6$Majority_protein_IDs_u<-unlist(lapply(lapply(temp_list, function(x) unlist(lapply(strsplit(x,"[:]"), function(y) y[1]))),function(z) z[1]))
#data_proteomics_mag6$Majority_protein_IDs_all<-unlist(lapply(temp_list, function(x) strsplit(x,"[:]")))

MAG6_pan2_only_prot<-merge(MAG6,data_proteomics_mag6,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)
colors<-c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')
ggplot(MAG6_pan2_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=bin_name))+geom_point()+scale_color_manual(values=colors)

MAG6_pan2_only_prot$log10_LFQ_intensity_Cave_Biofilm<- rowMeans(MAG6_pan2_only_prot[,33:34],na.rm = T)
MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm<- rowMeans(MAG6_pan2_only_prot[,30:31],na.rm = T)
MAG6_pan2_only_prot_reduced<-MAG6_pan2_only_prot[apply(MAG6_pan2_only_prot[,33:34], 1, function(x) ( all(x>0))),]
library(ggbeeswarm)
ggplot(MAG6_pan2_only_prot_reduced,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
MAG6_pan2_only_prot_reduced$gene_cluster_perc<-MAG6_pan2_only_prot_reduced$num_genomes_gene_cluster_has_hits/13*100
pp2<-ggplot(MAG6_pan2_only_prot_reduced,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
   geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
  #  geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "",title = "Ferroplasma MAG6")
# - Cave biofilm
pp2



MAG6_pan2_only_prot[MAG6_pan2_only_prot$KOfam_ACC.x%in%all,]
View(MAG6_pan2_only_prot[order(MAG6_pan2_only_prot$log10_LFQ_intensity_Cave_Biofilm,decreasing = T),])

##########Myco
temp_Myc<-read.cdhit.clstr("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/connect_proteomics/Mycobacterium_connect_100.clstr")

matchIDs<-rep(NA,nrow(Myco_Pan))
for(i in unique(temp_Myc$Cluster)){
  temp1<-temp_Myc[which(i==temp_Myc$Cluster),]
  if(nrow(temp1)==2){
    idx<-which(Myco_Pan$gene_callers_id==as.integer(temp1$Seq.Name[1]))
    matchIDs[idx]<-temp1$Seq.Name[2]
  }
}
Myco_Pan$matchIDs<-matchIDs
strsplit(data_proteomics_M_meth$`Majority protein IDs`,";")
data_proteomics_M_meth$Majority_protein_IDs_u<-unlist(lapply(strsplit(data_proteomics_M_meth$`Majority protein IDs`,";"), function(x) substr(x[grep("KDJLIKBO",x)], 1, 14)[1]))
data_proteomics_M_meth$Majority_protein_IDs_all<-unlist(lapply(strsplit(data_proteomics_M_meth$`Majority protein IDs`,";"), function(x) paste(substr(x[grep("KDJLIKBO",x)], 1, 14),collapse = ";")))

Myco_pan_only_prot<-merge(Myco_Pan,data_proteomics_M_meth,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)

ggplot(Myco_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2))+geom_point()
Myco_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm<- rowMeans(Myco_pan_only_prot[,33:34],na.rm = T)
Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm<- rowMeans(Myco_pan_only_prot[,30:31],na.rm = T)
Myco_pan_only_prot_reduced<-Myco_pan_only_prot[apply(Myco_pan_only_prot[,33:34], 1, function(x) ( all(x>0))),]
ggplot(Myco_pan_only_prot_reduced,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
View(Myco_pan_only_prot_reduced[order(Myco_pan_only_prot_reduced$log10_LFQ_intensity_Cave_Biofilm,decreasing = T),])
Myco_pan_only_prot_reduced$gene_cluster_perc<-Myco_pan_only_prot_reduced$num_genomes_gene_cluster_has_hits/69*100
SMMO<-c("K16157","K16158","K16159","K16160","K16161","K16162")
pp3<-ggplot(Myco_pan_only_prot_reduced,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
#   geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
#   geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
#    geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
#    geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "",title = "M.methanotrophicum")
# - Cave biofilm
pp3

pp3+pp1+ theme(legend.position = "none")+pp2+ theme(legend.position = "none")+ plot_layout(guides = "collect")


# 02010 ABC transporters [PATH:ko02010]
# 02060 Phosphotransferase system (PTS) [PATH:ko02060]
# 03070 Bacterial secretion system [PATH:ko03070]

i<-which(pws_mod_names$Pathway=="ko02010")# 
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")


i<-which(pws_mod_names$Pathway=="ko02060")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")

i<-which(pws_mod_names$Pathway=="ko03070")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")


i<-which(pws_mod_names$Pathway=="ko00500")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")



i<-which(pws_mod_names$Pathway=="ko00680")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")

i<-which(pws_mod_names$Pathway=="ko00920")# photosynth
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(MAG6_pan2_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( unlist(strsplit(pws_mod_names[i,4],",")),Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")


aa_kos<-unique(unlist(strsplit(paste(AAmod_names$V4,collapse = ","),",|[+]")))
paste(unique(MAG6_pan2_only_prot_reduced[match(aa_kos,MAG6_pan2_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
paste(unique(Myco_pan_only_prot_reduced[match(aa_kos,Myco_pan_only_prot_reduced$KOfam_ACC  ),24]),collapse = "+")
paste(unique(Ferr_pan_only_prot_reduced[match( aa_kos,Ferr_pan_only_prot_reduced$KOfam_ACC  ),16]),collapse = "+")
