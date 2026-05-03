################ combine proteomics
# metaproteomics were mapped against 92103 proteins from MAGs and virous 
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

dim(data_proteomics)
#[1] 3276    9

all_prot_id <- unique(c(unlist(strsplit(data_proteomics$`Majority protein IDs`,";")),unlist(strsplit(data_proteomics$`Protein IDs`,";"))))
length(all_prot_id)
# [1] 5352

length(grep(";",data_proteomics$`Majority protein IDs`))
#[1] 630

View(data_proteomics_tax)
colors_pB<-c('#a6d854',"#fb8072",'#8da0cb','#66c2a5',"grey")

info_table_DNA <- read.csv("~/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/info_table_DNA.csv")
info_table_DNA$Alternative_ID2<-gsub("_","",info_table_DNA$Alternative_ID2)
data_proteomics$ID<-unlist(lapply(strsplit(data_proteomics$`Protein IDs`,"_"),function(x) x[1]))
data_proteomics_tax<-merge(data_proteomics,info_table_DNA,by.x = "ID",by.y = "Alternative_ID2")
data_proteomics_tax$Genus1<-data_proteomics_tax$Genus

data_proteomics_tax$Genus1[!data_proteomics_tax$Genus%in%names(which(table(data_proteomics_tax$Genus)>100))]<-"Others"
data_proteomics_tax$Genus1<-as.factor(data_proteomics_tax$Genus1)
library(plotly)
axx <- list(
  title = "Cave Biofilm 1"
)

axy <- list(
  title = "Cave Biofilm 2"
)

axz <- list(
  title = "Culture"
)
data_proteomics_tax$cave_only<-rep(18,nrow(data_proteomics_tax))
data_proteomics_tax$cave_only[data_proteomics_tax$LFQ_intensity_Cave_Biofilm_1>0 & data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2>0]<-16

g3d<-plot_ly(data_proteomics_tax, x = ~log10_LFQ_intensity_Cave_Biofilm_1, 
             y = ~log10_LFQ_intensity_Cave_Biofilm_2, z = ~log10_LFQ_intensity_Culture, color = ~Genus1, colors = colors_pB,
             type = "scatter3d", mode = "markers",marker = list(size = 5))%>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))%>%
  layout(title = 'log10 LFQ intensity', legend = list(title=list(text='<b> Genus </b>')))

df1 <- data_proteomics_tax[data_proteomics_tax$cave_only == "square", ]
df2 <- data_proteomics_tax[data_proteomics_tax$cave_only != "square", ]

f<-plot_ly() %>%
  add_trace(x = ~df2$log10_LFQ_intensity_Cave_Biofilm_1, 
            y = ~df2$log10_LFQ_intensity_Cave_Biofilm_2, z = ~df2$log10_LFQ_intensity_Culture, color = ~df2$Genus1, colors = colors_pB,
            type = "scatter3d", mode = "markers",marker = list(symbol = "circle",size = 5)) %>%
  add_trace(x = ~df1$log10_LFQ_intensity_Cave_Biofilm_1, 
            y = ~df1$log10_LFQ_intensity_Cave_Biofilm_2, z = ~df1$log10_LFQ_intensity_Culture, color = ~df1$Genus1, colors = colors_pB,
            type = "scatter3d", mode = "markers",marker = list(symbol = "square",size = 3))%>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))%>%
  layout(title = 'log10 LFQ intensity', legend = list(title=list(text='<b> Genus </b>')))
f

# size indicate the protein identification on both biofilm samples  
library(scatterplot3d)
colors <- colors_pB[as.numeric(data_proteomics_tax$Genus1)]
size <- c(2,1)[factor(data_proteomics_tax$cave_only)]
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
scatterplot3d(x=data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_1, y=data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2, z=data_proteomics_tax$log10_LFQ_intensity_Culture,
                 main="",pch = data_proteomics_tax$cave_only,color=colors, box=F,cex.symbols=size,
              xlab="Cave Biofilm 1", ylab="Cave Biofilm 3", zlab="Culture",sub="log10 LFQ intensity")
addgrids3d(x=data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_1, y=data_proteomics_tax$log10_LFQ_intensity_Cave_Biofilm_2, z=data_proteomics_tax$log10_LFQ_intensity_Culture, grid = c("xy", "xz", "yz"))
#plot(data_proteomics$log10_LFQ_intensity_Cave_Biofilm_1,data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2)
dim(data_proteomics)

data_proteomics1<-data_proteomics[which(data_proteomics$LFQ_intensity_Cave_Biofilm_1>0 & data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2>0),]
dim(data_proteomics1)
# 2054
plot(data_proteomics1$log10_LFQ_intensity_Cave_Biofilm_1,data_proteomics1$log10_LFQ_intensity_Cave_Biofilm_2)
data_proteomics1$ID<-unlist(lapply(strsplit(data_proteomics1$`Protein IDs`,"_"),function(x) x[1]))


sum_prot<-as.data.frame(table(data_proteomics1$ID))
sum_prot_tax<-merge(sum_prot,info_table_DNA,by.x = "Var1",by.y = "Alternative_ID2")
sum_prot_tax<-sum_prot_tax[-2,]

sum_prot_tax$Var<-as.character(sum_prot_tax$Var)
sum_prot_tax$Var[1:5]<-c("M. MAG2", "M. MAG3","Cuniculiplasma c.", "M. Methanotrophicum", "Ferroplasma c.")
sum_prot_tax$Var<-gsub("SC","",sum_prot_tax$Var)
sum_prot_tax$Var <- factor(sum_prot_tax$Var, levels = as.character(sum_prot_tax$Var[order(sum_prot_tax$Freq,decreasing = T)]))

colors_p<-c('#a6d854','#8da0cb','#66c2a5',"grey")
sum_prot_tax$Genus1<-sum_prot_tax$Genus
sum_prot_tax$Genus1[which(sum_prot_tax$Freq<49)]<-"Others"

gb<-ggplot(sum_prot_tax, aes(x = Var, y = Freq,fill=Genus1)) +
  geom_bar(stat = "identity") +
  labs(title = "",x="MAGs",y="Identified Proteins in Biofilm") +scale_fill_manual(values=colors_p)+
  theme_minimal()+scale_x_discrete(guide = guide_axis(angle = 45)) + guides(fill=guide_legend(title="Genus"))+theme(text = element_text(size = 14),axis.text.x = element_text(size = 8))#+theme(axis.text.x = element_blank()) 
gb

################# combine proteomics and metagenomics 
# get abundance from figure01 mean_cov_all
mean_cov_all <- read.csv("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/DNA_mean_cov_all.csv")[,-1]
data_proteomics_formated<-data_proteomics1 %>% group_by(ID) %>% summarise(sum_log10_LFQ_intensity_Cave_Biofilm_1 = sum(log10_LFQ_intensity_Cave_Biofilm_1)
                                                                       ,sum_log10_LFQ_intensity_Cave_Biofilm_2 = sum(log10_LFQ_intensity_Cave_Biofilm_2),
                                                                       sum_log10_LFQ_intensity_Culture = sum(log10_LFQ_intensity_Culture),
                                                                       mean_log10_LFQ_intensity_Cave_Biofilm_1 = mean(log10_LFQ_intensity_Cave_Biofilm_1)
                                                                       ,mean_log10_LFQ_intensity_Cave_Biofilm_2 = mean(log10_LFQ_intensity_Cave_Biofilm_2),
                                                                       mean_log10_LFQ_intensity_Culture = mean(log10_LFQ_intensity_Culture),
                                                                       sum_LFQ_intensity_Cave_Biofilm_1 = sum(LFQ_intensity_Cave_Biofilm_1)
                                                                       ,sum_LFQ_intensity_Cave_Biofilm_2 = sum(LFQ_intensity_Cave_Biofilm_2),
                                                                       sum_LFQ_intensity_Culture = sum(log10_LFQ_intensity_Culture))

# data_proteomics_formated$mean_sum_log10_LFQ_intensity_Cave_Biofilm<-rowMeans(data_proteomics_formated[,2:3])
data_proteomics_formated$mean_mean_log10_LFQ_intensity_Cave_Biofilm<-rowMeans(data_proteomics_formated[,5:6])
# data_proteomics_formated$mean_mean_LFQ_intensity_Cave_Biofilm<-rowMeans(data_proteomics_formated[,5:6])


data_comb<-cbind(mean_cov_all[match(c("SC_MAG_00006","SC_MAG_00016","SC_MAG_00008","SC_MAG_00004"),mean_cov_all$bins),],
data_proteomics_formated[match(c("SCMAG00006","KDJLIKBO","KNPMNEEE","SCMAG00004"),data_proteomics_formated$ID),])
colors_p1<-c('#a6d854','#8da0cb','#66c2a5')
data_comb$DNA_mean_ab<-rowMeans(data_comb[,5:6])
data_comb$highlight<-c("Ferroplasma MAG6",    "M. methanotrophicum", "Ferroplasma c." ,   "Acidithiobacillus" )
data_comb$genus<-factor(data_comb$genus)
prot_dna<-ggplot(data_comb, aes(DNA_mean_ab,mean_mean_log10_LFQ_intensity_Cave_Biofilm,color=genus))+
  geom_point(size=8)+scale_color_manual(values=colors_p1)+
  geom_label_repel(aes(label=highlight), size=6, color="black",fill = alpha(c("white"),0.5))+labs(x="mean DNA coverage (log10)",y="mean LFQ intenisty (log10)")+theme(text = element_text(size = 14))
prot_dna  







##################
#################
#methane, sulfur,photosynth, Co2 fixation
library(KEGGREST)
kegg_KOs<-c("ko00680","ko00920","ko00710","ko00720")

pw_info<-list()
for(i in 1:length(kegg_KOs)){
  tryCatch({
    pw_info[[i]]<-keggGet(kegg_KOs[i])
  }, error=function(e){})

    module<-list()
    for (j in 1:length(pw_info[[i]][[1]]$MODULE)){
      module[[j]]<-keggGet(names(pw_info[[i]][[1]]$MODULE)[j])
    }
    KOs<-lapply(module, function(x) x[[1]]$ORTHOLOGY)
    names(KOs)<-lapply(module, function(x) x[[1]]$ENTRY)
    pw_info[[i]][[1]]$ORTHOLOGY_MODULE<-unlist(KOs)
  
  
}

pws_names<-c()
for(i in 1:length(pw_info)){
  temp<-c(
    unlist(paste(unlist(names(pw_info[[i]][[1]]$PATHWAY)),collapse = ",")),
    unlist(paste(unlist(pw_info[[i]][[1]]$PATHWAY),collapse = ",")),
    length(unique(pw_info[[i]][[1]]$ORTHOLOGY))
  )
  #print(temp)
  pws_names<-rbind(pws_names,temp)
}
pws_names<-as.data.frame(pws_names)
colnames(pws_names)<-c("ID","PATHWAY","Number-KOs")

list_mod<-list()
names_mod<-c()
h<-1
for( i in 1:length(pw_info)){
  mod<-unique(unlist(lapply(strsplit(names(pw_info[[i]][[1]]$ORTHOLOGY_MODULE),"[.]"), function(x) x[1])))
  for(j in 1:length(mod)){
    names_temp<-names(pw_info[[i]][[1]]$ORTHOLOGY_MODULE)[grepl( mod[j],names(pw_info[[i]][[1]]$ORTHOLOGY_MODULE))]
    names_temp1<-gsub(paste(mod[j],".",sep = ""),"",names_temp)
    names_ko<-unlist(strsplit(names_temp1,",|[+]|-"))
    list_mod[[h]]<-unlist(lapply(strsplit(names_ko," "), function(x) x[1]))
    h<-h+1
  }
  names_mod<-c(names_mod,paste(mod,pw_info[[i]][[1]]$NAME,sep = "."))
}
names(list_mod)<-names_mod



data_proteomics_ferro<-data_proteomics1[grepl("KNPMNEEE",data_proteomics1$`Majority protein IDs`),]
data_proteomics_mag6<-data_proteomics1[grepl("SCMAG00006",data_proteomics1$`Majority protein IDs`),]
data_proteomics_M_meth<-data_proteomics1[grepl("KDJLIKBO",data_proteomics1$`Majority protein IDs`),]



########## connect pan-genomics with 
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
Ferroplasma_Pan_2_gene_clusters_summary.txt <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/Ferroplasma/Ferroplasma_Pan_2/SUMMARY_few_classes/Ferroplasma_Pan_2_gene_clusters_summary.txt.gz")
Ferr<-Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name=="SFerroplasmacircular"),]
Ferr<-Ferr[order(Ferr$gene_callers_id,decreasing = F),]

MAG6<-Ferroplasma_Pan_2_gene_clusters_summary.txt[which(Ferroplasma_Pan_2_gene_clusters_summary.txt$genome_name=="SMAG00006"),] 
MAG6<-MAG6[order(MAG6$gene_callers_id,decreasing = F),]

Myco_Pan_gene_clusters_summary.txt <- read.delim("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/Genomes_selection_clean/mycobacterium/MYCO-SUMMARY/MYCO_Pan_gene_clusters_summary.txt.gz")
Myco_Pan<-Myco_Pan_gene_clusters_summary.txt[which(Myco_Pan_gene_clusters_summary.txt$genome_name=="M_methanotrophicum"),]
Myco_Pan<-Myco_Pan[order(Myco_Pan$gene_callers_id,decreasing = F),]


#####################################

matchIDs<-rep(NA,nrow(Ferr))
for(i in unique(temp_Fer$Cluster)){
  temp1<-temp_Fer[which(i==temp_Fer$Cluster),]
  if(nrow(temp1)==2){
    idx<-which(Ferr$gene_callers_id==as.integer(temp1$Seq.Name[1]))
    matchIDs[idx]<-temp1$Seq.Name[2]
  }
}
Ferr$matchIDs<-matchIDs

#strsplit(data_proteomics_ferro$`Majority protein IDs`,";")
data_proteomics_ferro$Majority_protein_IDs_u<-unlist(lapply(strsplit(data_proteomics_ferro$`Majority protein IDs`,";"), function(x) substr(x[grep("KNPMNEEE",x)], 1, 14)[1]))
#data_proteomics_ferro$Majority_protein_IDs_all<-unlist(lapply(strsplit(data_proteomics_ferro$`Majority protein IDs`,";"), function(x) paste(substr(x[grep("KNPMNEEE",x)], 1, 14),collapse = ";")))

Ferr_pan_only_prot<-merge(Ferr,data_proteomics_ferro,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)
Ferr_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm<-rowMeans(Ferr_pan_only_prot[,c(32,33)],na.rm = T)
Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm<-rowMeans(Ferr_pan_only_prot[,c(29,30)],na.rm = T)

ggplot(Ferr_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=bin_name))+geom_point()

#Ferr_pan_only_prot_reduced<-Ferr_pan_only_prot[apply(Ferr_pan_only_prot[,33:34], 1, function(x) ( all(x>0))),]
library(ggbeeswarm)
ggplot(Ferr_pan_only_prot,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
Ferr_pan_only_prot$gene_cluster_perc<-Ferr_pan_only_prot$num_genomes_gene_cluster_has_hits/13*100

#write.csv(Ferr_pan_only_prot,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/Combined_data/Ferr_pan_only_prot.csv")
methane_metabolism<-names(pw_info[[1]][[1]]$ORTHOLOGY)
sulfor_metabolism<-names(pw_info[[2]][[1]]$ORTHOLOGY)
co2_fix_photo<-names(pw_info[[3]][[1]]$ORTHOLOGY)
co2_fix_pro<-names(pw_info[[4]][[1]]$ORTHOLOGY)
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
pp1<-ggplot(Ferr_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
  # geom_label_repel(data=subset(Ferr_pan_only_prot_reduced, Ferr_pan_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
  geom_label_repel(data=subset(Ferr_pan_only_prot, Ferr_pan_only_prot$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "x/y axes shows the protein LFQ intensity (log10)",title="Ferroplasma c.")
# - Cave biofilm
pp1
data_temp<-Ferr_pan_only_prot
ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc))+geom_point(alpha=0.8)+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#ff7f00",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='sor')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='Sulf')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%succ_succ_coac), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='succ_succ_coac')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%oxoglutare), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='oxoglutare')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%isocitrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='isocitrate')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%citrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='citrate')+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#33a02c",label='SMMO')+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%methane_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#b2df8a")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%sulfor_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#e31a1c")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_photo), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#1f78b4")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_pro), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#cab2d6")


##################


temp_mag6<-read.cdhit.clstr("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/PANGENOMES/Thermoplasmatales_order/connect_proteomics/MAG6_connect_100.clstr")
matchIDs<-rep(NA,nrow(MAG6))
for(i in unique(temp_mag6$Cluster)){
  temp1<-temp_mag6[which(i==temp_mag6$Cluster),]
  if(nrow(temp1)==2){
    idx<-which(MAG6$gene_callers_id==as.integer(temp1$Seq.Name[1]))
    matchIDs[idx]<-temp1$Seq.Name[2]
  }
}
MAG6$matchIDs<-gsub("SC_MAG_00006_","SCMAG00006_",matchIDs)
strsplit(data_proteomics_mag6$`Majority protein IDs`,";")

#data_proteomics_mag6$Majority_protein_IDs_u<-unlist(lapply(strsplit(data_proteomics_mag6$`Majority protein IDs`,";"), function(x) gsub("SCMAG00006_","",x[grep("SCMAG00006_",x)])[1]))
data_proteomics_mag6$Majority_protein_IDs_u<-unlist(lapply(strsplit(data_proteomics_mag6$`Majority protein IDs`,";"), function(x) x[grep("SCMAG00006_",x)][1]))

# temp_list<-lapply(strsplit(data_proteomics_mag6$`Majority protein IDs`,";"), function(x) x[grep("SC_MAG_00006",x)])
# data_proteomics_mag6$Majority_protein_IDs_u<-unlist(lapply(lapply(temp_list, function(x) unlist(lapply(strsplit(x,"[:]"), function(y) y[1]))),function(z) z[1]))
#data_proteomics_mag6$Majority_protein_IDs_all<-unlist(lapply(temp_list, function(x) strsplit(x,"[:]")))

MAG6_pan2_only_prot<-merge(MAG6,data_proteomics_mag6,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)

colors<-c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f')
ggplot(MAG6_pan2_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=bin_name))+geom_point()+scale_color_manual(values=colors)

MAG6_pan2_only_prot$log10_LFQ_intensity_Cave_Biofilm<- rowMeans(MAG6_pan2_only_prot[,32:33],na.rm = T)
MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm<-rowMeans(MAG6_pan2_only_prot[,c(29,30)],na.rm = T)
# MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm<- rowMeans(MAG6_pan2_only_prot[,30:31],na.rm = T)
#MAG6_pan2_only_prot_reduced<-MAG6_pan2_only_prot[apply(MAG6_pan2_only_prot[,33:34], 1, function(x) ( all(x>0))),]
library(ggbeeswarm)
#write.csv(MAG6_pan2_only_prot,"/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/Combined_data/MAG6_pan2_only_prot.csv")

ggplot(MAG6_pan2_only_prot,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
MAG6_pan2_only_prot$gene_cluster_perc<-MAG6_pan2_only_prot$num_genomes_gene_cluster_has_hits/13*100
pp2<-ggplot(MAG6_pan2_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  # geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
  geom_label_repel(data=subset(MAG6_pan2_only_prot, MAG6_pan2_only_prot$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
  #  geom_label_repel(data=subset(MAG6_pan2_only_prot_reduced, MAG6_pan2_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "",title = "Ferroplasma MAG6")
# - Cave biofilm
pp2
data_temp<-MAG6_pan2_only_prot
ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc))+geom_point(alpha=0.8)+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#ff7f00",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='sor')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='Sulf')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%succ_succ_coac), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='succ_succ_coac')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%oxoglutare), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='oxoglutare')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%isocitrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='isocitrate')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%citrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='citrate')+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#33a02c",label='SMMO')+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%methane_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#b2df8a")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%sulfor_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#e31a1c")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_photo), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#1f78b4")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_pro), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#cab2d6")




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
#data_proteomics_M_meth$Majority_protein_IDs_all<-unlist(lapply(strsplit(data_proteomics_M_meth$`Majority protein IDs`,";"), function(x) paste(substr(x[grep("KDJLIKBO",x)], 1, 14),collapse = ";")))

Myco_pan_only_prot<-merge(Myco_Pan,data_proteomics_M_meth,by.x = "matchIDs",by.y="Majority_protein_IDs_u",all==T)

ggplot(Myco_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2))+geom_point()

Myco_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm<- rowMeans(Myco_pan_only_prot[,32:33],na.rm = T)
Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm<- rowMeans(Myco_pan_only_prot[,29:30],na.rm = T)

#Myco_pan_only_prot_reduced<-Myco_pan_only_prot[apply(Myco_pan_only_prot[,33:34], 1, function(x) ( all(x>0))),]
ggplot(Myco_pan_only_prot,aes(num_genomes_gene_cluster_has_hits,log10_LFQ_intensity_Cave_Biofilm))+geom_beeswarm()
View(Myco_pan_only_prot[order(Myco_pan_only_prot$log10_LFQ_intensity_Cave_Biofilm,decreasing = T),])
Myco_pan_only_prot$gene_cluster_perc<-Myco_pan_only_prot$num_genomes_gene_cluster_has_hits/69*100
SMMO<-c("K16157","K16158","K16159","K16160","K16161","K16162")
pp3<-ggplot(Myco_pan_only_prot,aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2,color=gene_cluster_perc))+geom_point(alpha=0.8)+
  scale_x_continuous(limits = c(6, 11)) + scale_y_continuous(limits = c(6, 11))+
  #   geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#fdc086",label='Sulfur')+
  #   geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#ff7f00",label='sor')+
  geom_label_repel(data=subset(Myco_pan_only_prot, Myco_pan_only_prot$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#984ea3",label='CRISPR')+
  #    geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%rTCA), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#e41a1c",label='rTCA')+
  #    geom_label_repel(data=subset(Myco_pan_only_prot_reduced, Myco_pan_only_prot_reduced$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm_1,log10_LFQ_intensity_Cave_Biofilm_2), size=3, color="#7fc97f",label='sMMO')+
  labs(x="Cave biofilm 1",y="Cave biofilm 2",color= "GC (%)",caption = "",title = "M. methanotrophicum")
# - Cave biofilm
pp3



data_temp<-Myco_pan_only_prot
ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc))+geom_point(alpha=0.8)+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#ff7f00",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%sor), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='sor')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%Sulf), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='Sulf')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%succ_succ_coac), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='succ_succ_coac')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%oxoglutare), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='oxoglutare')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%isocitrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='isocitrate')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%citrate), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#984ea3",label='citrate')+
  geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#33a02c",label='sMMO')+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%methane_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#b2df8a")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%sulfor_metabolism), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#e31a1c")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_photo), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#1f78b4")+
  geom_point(data=subset(data_temp, data_temp$KOfam_ACC%in%co2_fix_pro), aes(log10_LFQ_intensity_Cave_Biofilm,gene_cluster_perc), size=3, color="#cab2d6")

pp3+pp1+ theme(legend.position = "none")+pp2+ theme(legend.position = "none")+ plot_layout(guides = "collect")

###################################
###################################

# pan-genome and proteom highlights  FIGURE4

###################################
###################################

num_top<-5
data_temp<-Myco_pan_only_prot

#ribosomal
length(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])

#go per module
selected_mod<-rep(NA,nrow(data_temp))
for (i in 1:length(list_mod)){
  selected_mod[data_temp$KOfam_ACC%in%list_mod[[i]]]<-names(list_mod)[i]
}
data_temp$test_mod<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[1]))
data_temp$test_color<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[2]))

temp_list_pw<-list(methane_metabolism=names(pw_info[[1]][[1]]$ORTHOLOGY), sulfor_metabolism=names(pw_info[[2]][[1]]$ORTHOLOGY),co2_fix_photo=names(pw_info[[3]][[1]]$ORTHOLOGY),
                   co2_fix_pro=names(pw_info[[4]][[1]]$ORTHOLOGY), CRISPR=CRISP)
selected_pw<-rep("Other",nrow(data_temp))
for (i in 1:length(temp_list_pw)){
  selected_pw[data_temp$KOfam_ACC%in%temp_list_pw[[i]]]<-names(temp_list_pw)[i]
}
data_temp$selected_pw<-selected_pw
#myco
colors_pw<-c('#e7298a','#7570b3','#1b9e77',"lightgrey",'#d95f02')
# shapes<-c(15,17,19,18,21,8)
shapes<-c(15,17,18,21,8)
data_temp$gene_cluster_perc_f<-factor(round(data_temp$gene_cluster_perc,digits = 2))
library(arules)
step<-(max(data_temp$num_genomes_gene_cluster_has_hits)-2)/3
step<-ceiling(step)
temp_group<-data_temp$num_genomes_gene_cluster_has_hits
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>1 & data_temp$num_genomes_gene_cluster_has_hits <= step)]<-"2-23"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>step & data_temp$num_genomes_gene_cluster_has_hits <= 2*step)]<-"24-46"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>2*step & data_temp$num_genomes_gene_cluster_has_hits < 3*step)]<-"47-68"

  
data_temp$num_genomes_gene_cluster_has_hits_f<-factor(temp_group)
#data_temp$num_genomes_gene_cluster_has_hits_f<-discretize( data_temp$num_genomes_gene_cluster_has_hits,breaks = 5,method = "interval")


data_temp$KOfam_clean<-gsub(" [[].*","",data_temp$KOfam)
data_temp$KOfam_clean<-gsub("-associated protein|","",data_temp$KOfam_clean)
data_temp$KOfam_clean<-gsub("methane monooxygenase component","sMMO",data_temp$KOfam_clean)
data_temp$KOfam_clean[data_temp$KOfam_clean==""]<-"unknown"
pp1<-ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,color=selected_pw,shape=selected_pw))+ geom_beeswarm(size=6,alpha=0.7)+scale_color_manual(values=colors_pw)+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP),max.overlaps=Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="purple",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="#33a02c",label='sMMO')+
  coord_cartesian(clip = "off") +xlim(6,11)+
#  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$KOfam%in%data_temp$KOfam[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:5]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=5, color="black")+
  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$matchIDs%in%data_temp$matchIDs[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:num_top]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=6, color="black")+
  labs(x="",y="Presence of gene in pan-genome",color="Selected pathways",shape="Selected pathways",title = "M. methanotrophicum")+
  theme(text = element_text(size = 18))+scale_shape_manual(values=shapes)+
  theme(legend.position = "bottom")
pp1 

data_temp<-Ferr_pan_only_prot

length(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])

#go per module
selected_mod<-rep(NA,nrow(data_temp))
for (i in 1:length(list_mod)){
  selected_mod[data_temp$KOfam_ACC%in%list_mod[[i]]]<-names(list_mod)[i]
}
data_temp$test_mod<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[1]))
data_temp$test_color<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[2]))

temp_list_pw<-list(methane_metabolism=names(pw_info[[1]][[1]]$ORTHOLOGY), sulfor_metabolism=names(pw_info[[2]][[1]]$ORTHOLOGY),co2_fix_photo=names(pw_info[[3]][[1]]$ORTHOLOGY),
                   co2_fix_pro=names(pw_info[[4]][[1]]$ORTHOLOGY), CRISPR=CRISP)
selected_pw<-rep("Other",nrow(data_temp))
for (i in 1:length(temp_list_pw)){
  selected_pw[data_temp$KOfam_ACC%in%temp_list_pw[[i]]]<-names(temp_list_pw)[i]
}
data_temp$selected_pw<-selected_pw

#ferr
colors_pw<-c('#e7298a','#7570b3',"#8dd3c7",'#1b9e77',"lightgrey",'#d95f02')
shapes<-c(15,17,19,18,21,8)

data_temp$gene_cluster_perc_f<-factor(round(data_temp$gene_cluster_perc,digits = 2))
library(arules)
step<-(max(data_temp$num_genomes_gene_cluster_has_hits)-2)/3
step<-ceiling(step)
temp_group<-data_temp$num_genomes_gene_cluster_has_hits
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>1 & data_temp$num_genomes_gene_cluster_has_hits <= step)]<-"2-4"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>step & data_temp$num_genomes_gene_cluster_has_hits <= 2*step)]<-"5-8"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>2*step & data_temp$num_genomes_gene_cluster_has_hits <= 3*step)]<-"9-12"


data_temp$num_genomes_gene_cluster_has_hits_f<-factor(temp_group,levels = c("1"  ,  "2-4" , "5-8" , "9-12", "13" ))

#data_temp$num_genomes_gene_cluster_has_hits_f<-discretize( data_temp$num_genomes_gene_cluster_has_hits,breaks = 20)
data_temp$KOfam_clean<-gsub(" [[].*","",data_temp$KOfam)
data_temp$KOfam_clean<-gsub("-associated protein|","",data_temp$KOfam_clean)
data_temp$KOfam_clean<-gsub("methane monooxygenase component","sMMO",data_temp$KOfam_clean)
data_temp$KOfam_clean[data_temp$KOfam_clean==""]<-"unknown"
pp2<-ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,color=selected_pw,shape=selected_pw))+ geom_beeswarm(size=6,alpha=0.7)+scale_color_manual(values=colors_pw)+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP),max.overlaps=Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="purple",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="#33a02c",label='sMMO')+
  coord_cartesian(clip = "off") +xlim(6,11)+
  #  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$KOfam%in%data_temp$KOfam[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:5]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=5, color="black")+
  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$matchIDs%in%data_temp$matchIDs[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:num_top]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=6, color="black")+
  labs(x="mean LFQ intensity (log10) Cave Biofilm",y="",color="Selected pathways",shape="Selected pathways",title = "Ferroplasma c.")+
  theme(text = element_text(size = 18))+scale_shape_manual(values=shapes)+
  theme(legend.position = "bottom")
pp2 



data_temp<-MAG6_pan2_only_prot

length(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Cave_Biofilm[grep("ribosom",data_temp$KOfam)])
mean(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])
sd(data_temp$log10_LFQ_intensity_Culture[grep("ribosom",data_temp$KOfam)])

#go per module
selected_mod<-rep(NA,nrow(data_temp))
for (i in 1:length(list_mod)){
  selected_mod[data_temp$KOfam_ACC%in%list_mod[[i]]]<-names(list_mod)[i]
}
data_temp$test_mod<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[1]))
data_temp$test_color<-unlist(lapply(strsplit(selected_mod,"[.]"), function(x) x[2]))

temp_list_pw<-list(methane_metabolism=names(pw_info[[1]][[1]]$ORTHOLOGY), sulfor_metabolism=names(pw_info[[2]][[1]]$ORTHOLOGY),co2_fix_photo=names(pw_info[[3]][[1]]$ORTHOLOGY),
                   co2_fix_pro=names(pw_info[[4]][[1]]$ORTHOLOGY), CRISPR=CRISP)
selected_pw<-rep("Other",nrow(data_temp))
for (i in 1:length(temp_list_pw)){
  selected_pw[data_temp$KOfam_ACC%in%temp_list_pw[[i]]]<-names(temp_list_pw)[i]
}
data_temp$selected_pw<-selected_pw

#mag6
colors_pw<-c('#e7298a','#7570b3',"#8dd3c7",'#1b9e77',"lightgrey",'#d95f02')
shapes<-c(15,17,19,18,21,8)
data_temp$gene_cluster_perc_f<-factor(round(data_temp$gene_cluster_perc,digits = 2))
library(arules)
step<-(max(data_temp$num_genomes_gene_cluster_has_hits)-2)/3
step<-ceiling(step)
temp_group<-data_temp$num_genomes_gene_cluster_has_hits
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>1 & data_temp$num_genomes_gene_cluster_has_hits <= step)]<-"2-4"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>step & data_temp$num_genomes_gene_cluster_has_hits <= 2*step)]<-"5-8"
temp_group[which(data_temp$num_genomes_gene_cluster_has_hits>2*step & data_temp$num_genomes_gene_cluster_has_hits <= 3*step)]<-"9-12"


data_temp$num_genomes_gene_cluster_has_hits_f<-factor(temp_group,levels = c("1"  ,  "2-4" , "5-8" , "9-12", "13" ))
#data_temp$num_genomes_gene_cluster_has_hits_f<-discretize( data_temp$num_genomes_gene_cluster_has_hits,breaks = 20)
data_temp$KOfam_clean<-gsub(" [[].*","",data_temp$KOfam)
data_temp$KOfam_clean<-gsub("-associated protein|","",data_temp$KOfam_clean)
data_temp$KOfam_clean<-gsub("methane monooxygenase component","sMMO",data_temp$KOfam_clean)
data_temp$KOfam_clean[data_temp$KOfam_clean==""]<-"unknown"
pp3<-ggplot(data_temp,aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,color=selected_pw,shape=selected_pw))+ geom_beeswarm(size=6,alpha=0.7)+scale_color_manual(values=colors_pw)+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%CRISP),max.overlaps=Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="purple",label='CRISPR')+
  # geom_label_repel(data=subset(data_temp, data_temp$KOfam_ACC%in%SMMO), aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f), size=3, color="#33a02c",label='sMMO')+
  coord_cartesian(clip = "off") +xlim(6,11)+
  #  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$KOfam%in%data_temp$KOfam[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:5]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=5, color="black")+
  geom_label_repel(xlim = c(-Inf, Inf), ylim = c(-Inf, Inf), fill = alpha(c("white"),0.7),data=subset(data_temp, data_temp$matchIDs%in%data_temp$matchIDs[order(data_temp$log10_LFQ_intensity_Cave_Biofilm,decreasing = T)[1:num_top]]),box.padding = 0.5, max.overlaps = Inf, aes(log10_LFQ_intensity_Cave_Biofilm,num_genomes_gene_cluster_has_hits_f,label=KOfam_clean), size=6, color="black")+
  labs(x="",y="",color="Selected pathways",shape="Selected pathways",title = "Ferroplasma MAG6")+
  theme(text = element_text(size = 18))+scale_shape_manual(values=shapes)+
  theme(legend.position = "bottom")
pp3 


pp_final<-pp1+theme(legend.position="none")+pp2+pp3+theme(legend.position="none")+
  plot_annotation(tag_levels = 'a')

pp_final
ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/Figure4_abc.pdf",pp_final,width = 20,height = 8)
###################################
#focus on defence 
##################
myco_defense_finder_systems <- read.delim2("~/Desktop/Projects/Mycobacterium_sulfur_cave/Defencefinder/myco_defense_finder_systems.tsv")
pall<-unlist(strsplit(myco_defense_finder_systems$protein_in_syst,","))
Myco_pan_only_prot[Myco_pan_only_prot$matchIDs%in%pall,]

Ferr_c_defense_finder_systems <- read.delim2("~/Desktop/Projects/Mycobacterium_sulfur_cave/Defencefinder/Ferr_c_defense_finder_systems.tsv")
pall<-unlist(strsplit(Ferr_c_defense_finder_systems$protein_in_syst,","))
Ferr_pan_only_prot[Ferr_pan_only_prot$matchIDs%in%pall,]

MAG6_c_defense_finder_systems <- read.delim2("~/Desktop/Projects/Mycobacterium_sulfur_cave/Defencefinder/MAG6_defense_finder_systems.tsv")
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


###############Fero
def_table<-c()
for(i in 1:nrow(Ferr_c_defense_finder_systems)){
  p<-unlist(strsplit(Ferr_c_defense_finder_systems$protein_in_syst[i],","))
  if(any(Ferr_pan_only_prot$matchIDs%in%p==T)){
    def_temp<-Ferr_pan_only_prot[ Ferr_pan_only_prot$matchIDs%in%p,]
    def_temp$system<-Ferr_c_defense_finder_systems$subtype[i]
  }else{
    def_temp<-  matrix(NA, nrow = length(p), ncol = ncol(Ferr_pan_only_prot)) 
    def_temp$system<-Ferr_c_defense_finder_systems$subtype[i]
  }
 
  def_table<-rbind(def_table,def_temp)
}
ferr_ribo<-Ferr_pan_only_prot[grepl("ribosomal",Ferr_pan_only_prot$KOfam),]
ferr_ribo$system<-"Ribosomal"
data_ferro_def<-rbind(def_table,ferr_ribo)

df_melted <- melt(data_ferro_def, measure.vars = c("log10_LFQ_intensity_Cave_Biofilm_1", "log10_LFQ_intensity_Cave_Biofilm_2"))
df_melted$Sample<-gsub("log10_LFQ_intensity_Cave_Biofilm_","Sample ",df_melted$variable)
#c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c')

df_melted$system<-gsub("_"," ",df_melted$system)
df_melted$system<-gsub("Class","C",df_melted$system)
df_melted$system<-gsub("Subtype","S",df_melted$system)

Ferroplasma_ribo_crisp<-ggplot(df_melted, aes(x = system, y = value)) +
  geom_jitter(aes(color = system, shape =Sample), size=6, width = 0.2, alpha = 0.5) +
  geom_boxplot(aes(fill = system), alpha = 0.6, outlier.shape = NA) +
  scale_color_manual(values = c( "#33a02c","#e31a1c", "#fb9a99")) +
  scale_fill_manual(values = c( "#33a02c","#e31a1c", "#fb9a99")) +
  geom_signif(comparisons = list(c("CAS C1-S-I-D", "Ribosomal")),
              map_signif_level = TRUE, textsize = 4) +ylim(6.3,10.5)+
  theme_minimal() +
  labs(title = "Ferroplasma c.",
       x = "",
       y = "",
       color = "Functional category",
       fill = "Functional category")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(text = element_text(size = 18))
###############################


###############################
data<-MAG6_pan2_only_prot
data$matchIDs<-gsub("SCMAG00006_","",data$matchIDs)
def_table<-c()
for(i in 1:nrow(MAG6_c_defense_finder_systems)){
  p<-unlist(strsplit(MAG6_c_defense_finder_systems$protein_in_syst[i],","))
  if(any(data$matchIDs%in%p==T)){
    def_temp<-data[ data$matchIDs%in%p,]
    def_temp$system<-MAG6_c_defense_finder_systems$subtype[i]
  }else{
    def_temp<-  as.data.frame(matrix(NA, nrow = length(p), ncol = ncol(data)) )
    if(nrow(def_temp)==1){
      def_temp[1,ncol(def_temp)+1]<-MAG6_c_defense_finder_systems$subtype[i]
    }else{
      def_temp$system<-MAG6_c_defense_finder_systems$subtype[i]
    }
 
  }
  colnames(def_temp)<-c(colnames(data), "system")
  def_table<-rbind(def_table,def_temp)
}
ferr_ribo<-data[grepl("ribosomal",data$KOfam),]
ferr_ribo$system<-"Ribosomal"
colnames(def_table)<-colnames(Myco_ribo)
data_ferro_def<-rbind(def_table,ferr_ribo)

df_melted <- melt(data_ferro_def, measure.vars = c("log10_LFQ_intensity_Cave_Biofilm_1", "log10_LFQ_intensity_Cave_Biofilm_2"))
df_melted$Sample<-gsub("log10_LFQ_intensity_Cave_Biofilm_","Sample ",df_melted$variable)
#c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',"#ff7f00")

df_melted$system<-gsub("_"," ",df_melted$system)
df_melted$system<-gsub("Class","C",df_melted$system)
df_melted$system<-gsub("Subtype","S",df_melted$system)

MAG6_ribo_crisp<-ggplot(df_melted, aes(x = system, y = value)) +
  geom_jitter(aes(color = system, shape =Sample), size=6, width = 0.2, alpha = 0.5) +
  geom_boxplot(aes(fill = system), alpha = 0.6, outlier.shape = NA) +
  scale_color_manual(values = c('#a6cee3','#33a02c','#1f78b4','#b2df8a','#fb9a99',"#ff7f00",'#e31a1c',"#ffff33")) +
  scale_fill_manual(values = c('#a6cee3','#33a02c','#1f78b4','#b2df8a',"#ff7f00",'#e31a1c','#fb9a99',"#ffff33")) +
  geom_signif(comparisons = list(c("CAS C1-S-I-D", "Ribosomal"),c("CAS C1-S-I-B", "Ribosomal"),c("CAS C1-S-I-G", "Ribosomal"),c("CAS C1-S-III-A", "Ribosomal")),
              map_signif_level = TRUE, textsize = 4,y_position = c(9.7,10,9.4,9.1)) +ylim(6.3,10.5)+
  theme_minimal() +
  labs(title = "Ferroplasma MAG6",
       x = "",
       y = "",
       color = "Functional category",
       fill = "Functional category")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+ guides(color = "none")+theme(text = element_text(size = 18))
MAG6_ribo_crisp

def_final<-myco+theme(legend.position="none")+Ferroplasma_ribo_crisp+theme(legend.position="none")+MAG6_ribo_crisp+theme(legend.position="none")+
  plot_annotation(tag_levels = 'a')

def_final

figure4 <- (pp1+theme(legend.position="none")+pp2+pp3+theme(legend.position="none"))/(myco+theme(legend.position="none")+Ferroplasma_ribo_crisp+theme(legend.position="none")+MAG6_ribo_crisp+theme(legend.position="none"))+
  plot_annotation(tag_levels = 'a')
ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/Figure4.pdf",figure4,width = 18,height = 12)
# def_final<-Ferroplasma_ribo_crisp+MAG6_ribo_crisp+ plot_layout(guides = "collect")
# 
# def_final
################################# visulize pathway enrichment 

################ sum lfg proteomics per pw
pws_mod_names <- read.csv("~/Desktop/Projects/Mycobacterium_sulfur_cave/METAPROTEOMICS/eggnog/pw_info.csv")
pws_mod_names_r<-pws_mod_names[-grep("Human|Organismal Systems",pws_mod_names$V5),]
pws_mod_names_r<-pws_mod_names_r[-grep(" - fly| - worm",pws_mod_names_r$V2),]

Ferr_data<-pws_mod_names_r
match<-apply(pws_mod_names_r, 1, function(x) sum(Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")) ,Ferr_pan_only_prot$KOfam_ACC) ] ,na.rm = T))
# match[is.nan(match)]<-0
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(Ferr_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")) ,Ferr_pan_only_prot$KOfam_ACC)]>0)))

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
  scale_colour_continuous(type = "viridis")+theme(text = element_text(size = 14))
gg1

###########################
MAG6_data<-pws_mod_names_r
match<-apply(pws_mod_names_r, 1, function(x) sum(MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm[match( unlist(strsplit(x[5],",")),MAG6_pan2_only_prot$KOfam_ACC ) ] ,na.rm = T))
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(MAG6_pan2_only_prot$LFQ_intensity_Cave_Biofilm[match( unlist(strsplit(x[5],",")),MAG6_pan2_only_prot$KOfam_ACC ) ]>0)))
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
  scale_colour_continuous(type = "viridis")+theme(text = element_text(size = 14))
gg2

####################


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
  scale_colour_continuous(type = "viridis")+theme(text = element_text(size = 14))
gg_ferr2


##########################
Myco_data<-pws_mod_names_r
match<-apply(pws_mod_names_r, 1, function(x) sum(Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")),Myco_pan_only_prot$KOfam_ACC ) ] ,na.rm = T))
NR_active<-apply(pws_mod_names_r, 1, function(x) length(which(Myco_pan_only_prot$LFQ_intensity_Cave_Biofilm[match(unlist(strsplit(x[5],",")),Myco_pan_only_prot$KOfam_ACC )  ]>0)))
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
  labs(y="Pathways (top10)",x="",title = "M. methanotrophicum",color="N. proteins")+
  scale_colour_continuous(type = "viridis")+theme(text = element_text(size = 14))
gg3

gg<-gg3+theme(legend.direction="horizontal",legend.position="top")+gg1+theme(legend.direction="horizontal",legend.position="top")+gg2+theme(legend.direction="horizontal",legend.position="top")
gg<-wrap_elements(gg) +
  labs(tag = "Sum protein LFQ intensity (log10) - Cave biofilm") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom"
  )
gg_f<-(gb+gb+prot_dna)/(gg3+theme(legend.direction="horizontal",legend.position="top")+gg1+theme(legend.direction="horizontal",legend.position="top")+gg2+theme(legend.direction="horizontal",legend.position="top") +
                          plot_annotation(tag_levels = 'a')+
                          labs(tag = "Sum protein LFQ intensity (log10) - Cave biofilm") +
                          theme(
                            plot.tag = element_text(size = rel(1)),
                            plot.tag.position = "bottom"
                          ))

ggsave("/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/FINAL/Figure3.pdf",gg_f,width = 20,height = 9)

################# 
################# investigating pathwys C02 fixation
################# 

co2fix_KO<-names(pw_info[[4]][[1]]$ORTHOLOGY)

t<-MAG6_pan2_only_prot$KOfam_ACC[MAG6_pan2_only_prot$KOfam_ACC%in%co2fix_KO]
paste(unique(t),collapse = "+")


f<-Ferr_pan_only_prot$KOfam_ACC[Ferr_pan_only_prot$KOfam_ACC%in%co2fix_KO]
paste(unique(f),collapse = "+")

m<-Myco_pan_only_prot$KOfam_ACC[Myco_pan_only_prot$KOfam_ACC%in%co2fix_KO]
paste(unique(m),collapse = "+")
Myco_pan_only_prot
#photo
co2fix_p_KO<-names(pw_info[[3]][[1]]$ORTHOLOGY)

t<-MAG6_pan2_only_prot$KOfam_ACC[MAG6_pan2_only_prot$KOfam_ACC%in%co2fix_p_KO]
paste(unique(t),collapse = "+")


f<-Ferr_pan_only_prot$KOfam_ACC[Ferr_pan_only_prot$KOfam_ACC%in%co2fix_p_KO]
paste(unique(f),collapse = "+")

#sulf
sulf_KO<-names(pw_info[[2]][[1]]$ORTHOLOGY)

t<-MAG6_pan2_only_prot$KOfam_ACC[MAG6_pan2_only_prot$KOfam_ACC%in%sulf_KO]
paste(unique(t),collapse = "+")


f<-Ferr_pan_only_prot$KOfam_ACC[Ferr_pan_only_prot$KOfam_ACC%in%sulf_KO]
paste(unique(f),collapse = "+")

#methane
#sulf
meth_KO<-names(pw_info[[1]][[1]]$ORTHOLOGY)
meth_KO<-unlist(strsplit(meth_KO,"[+]"))
t<-MAG6_pan2_only_prot$KOfam_ACC[MAG6_pan2_only_prot$KOfam_ACC%in%meth_KO]
paste(unique(t),collapse = "+")


f<-Ferr_pan_only_prot$KOfam_ACC[Ferr_pan_only_prot$KOfam_ACC%in%sulf_KO]
paste(unique(f),collapse = "+")
