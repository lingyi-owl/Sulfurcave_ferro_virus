ProteinGroups<-readxl::read_xlsx("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/proteomics/proteinGroups.xlsx",sheet = 2)
ids<-readxl::read_xlsx("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/proteomics/proteinGroups.xlsx",sheet = 1)
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




data_proteomics$ID<- gsub("_","",data_proteomics$`Majority protein IDs`)
sulfur_cave_ref_proteins_eggnog.emapper <- read.delim("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/metagenomics/Functional_annotation/sulfur_cave_ref_proteins_eggnog.emapper.annotations", header=FALSE, comment.char="#")
sulfur_cave_ref_proteins_eggnog.emapper$tax_id <- sapply(
  strsplit(sulfur_cave_ref_proteins_eggnog.emapper$V1, "_"),
  function(x) {
    if (length(x) >= 3 && x[1] == "SC" && x[2] %in% c("MAG", "Bin")) {
      paste0(x[1:3], collapse = "")  # SC_MAG_00014 → SCMAG00014
    } else if (length(x) >= 2 && x[1] == "Crude") {
      paste0(x[1:2], collapse = "")  # Crude_X → CrudeX
    } else {
      x[1]
    }
  }
)
sulfur_cave_ref_proteins_eggnog.emapper$ID <- gsub("_","",sulfur_cave_ref_proteins_eggnog.emapper$V1)

library(dplyr)
library(tidyr)

data_proteomics_prot <- data_proteomics %>%
  separate_rows(ID, sep = ";") %>%
  right_join(sulfur_cave_ref_proteins_eggnog.emapper, by = "ID")

ferr_mag6KO<- data_proteomics_prot$V12[data_proteomics_prot$tax_id == "SCMAG00006" & data_proteomics_prot$LFQ_intensity_Cave_Biofilm_1 >0]
# reverse TCA
ferroplasma_rtca_markers <- data.frame(
  gene = c(
    "porA", "porB", "porC", "porD",
    "korA/oorA", "korB/oorB", "korC/oorC", "korD/oorD",
    "aclA", "aclB",
    "ccs", "ccl",
    "frdA", "frdB", "frdC", "frdD",
    "icd", "acnA/acnB", "fum", "mdh", "sucC", "sucD",
    "ppsA", "ppdk"
  ),
  KO = c(
    "K00169", "K00170", "K00171", "K00172",
    "K00174", "K00175", "K00176", "K00177",
    "K15230", "K15231",
    "K15232", "K15233",
    "K00244", "K00245", "K00246", "K00247",
    "K00031", "K01681", "K01676", "K00024", "K01902", "K01903",
    "K01006", "K01007"
  ),
  enzyme = c(
    "Pyruvate:ferredoxin oxidoreductase subunit A",
    "Pyruvate:ferredoxin oxidoreductase subunit B",
    "Pyruvate:ferredoxin oxidoreductase subunit C",
    "Pyruvate:ferredoxin oxidoreductase subunit D",
    "2-oxoglutarate:ferredoxin oxidoreductase subunit A",
    "2-oxoglutarate:ferredoxin oxidoreductase subunit B",
    "2-oxoglutarate:ferredoxin oxidoreductase subunit C",
    "2-oxoglutarate:ferredoxin oxidoreductase subunit D",
    "ATP-citrate lyase subunit A",
    "ATP-citrate lyase subunit B",
    "Citryl-CoA synthetase",
    "Citryl-CoA lyase",
    "Fumarate reductase flavoprotein subunit",
    "Fumarate reductase iron-sulfur subunit",
    "Fumarate reductase membrane subunit C",
    "Fumarate reductase membrane subunit D",
    "Isocitrate dehydrogenase",
    "Aconitate hydratase / aconitase",
    "Fumarate hydratase",
    "Malate dehydrogenase",
    "Succinyl-CoA synthetase beta subunit",
    "Succinyl-CoA synthetase alpha subunit",
    "Phosphoenolpyruvate synthase",
    "Pyruvate phosphate dikinase"
  ),
  pathway_block = c(
    rep("rTCA core ", 8),
    rep("Citrate cleavage", 2),
    rep("Alternative citrate cleavage", 2),
    rep("Reductive succinate branch", 4),
    rep("rTCA backbone", 6),
    rep("Gluconeogenesis", 2)
  ),
  importance = c(
    "core", "core", "core", "core",
    "core", "core", "core", "core",
    "diagnostic", "diagnostic",
    "diagnostic", "diagnostic",
    "supportive", "supportive", "supportive", "supportive",
    "supportive", "supportive", "supportive", "supportive", "supportive", "supportive",
    "important", "important"
  ),
  ferroplasma_interpretation = c(
    "Not detected; weakens a full canonical rTCA claim",
    "Not detected; weakens a full canonical rTCA claim",
    "Not detected; weakens a full canonical rTCA claim",
    "Not detected; weakens a full canonical rTCA claim",
    "Detected; supports ferredoxin-dependent reductive carbon metabolism",
    "Detected; supports ferredoxin-dependent reductive carbon metabolism",
    "Not detected in current proteome list",
    "Not detected in current proteome list",
    "Not detected; canonical ATP-citrate lyase route not supported",
    "Not detected; canonical ATP-citrate lyase route not supported",
    "Not detected; alternative archaeal citrate cleavage route not supported",
    "Not detected; alternative archaeal citrate cleavage route not supported",
    "Not detected in current proteome list",
    "Not detected in current proteome list",
    "Not detected in current proteome list",
    "Not detected in current proteome list",
    "Detected; supports active TCA/rTCA-associated carbon metabolism",
    "Detected; supports active TCA/rTCA-associated carbon metabolism",
    "Detected; supports active TCA/rTCA-associated carbon metabolism",
    "Detected; supports active TCA/rTCA-associated carbon metabolism",
    "Detected; supports succinyl-CoA interconversion in central metabolism",
    "Detected; supports succinyl-CoA interconversion in central metabolism",
    "Detected; strong support for gluconeogenic carbon flow and biomass formation",
    "Detected; strong support for gluconeogenic carbon flow and biomass formation"
  ),
  comment = c(
    "Missing POR means pyruvate-forming reductive carboxylation is not yet demonstrated",
    "Missing POR means pyruvate-forming reductive carboxylation is not yet demonstrated",
    "Missing POR means pyruvate-forming reductive carboxylation is not yet demonstrated",
    "Missing POR means pyruvate-forming reductive carboxylation is not yet demonstrated",
    "Strongest detected rTCA-like signal in this dataset",
    "Strongest detected rTCA-like signal in this dataset",
    "Could still be present genomically but not detected proteomically",
    "Could still be present genomically but not detected proteomically",
    "Without citrate cleavage, full rTCA remains incomplete",
    "Without citrate cleavage, full rTCA remains incomplete",
    "Without citrate cleavage, full rTCA remains incomplete",
    "Without citrate cleavage, full rTCA remains incomplete",
    "Would support reductive branch if present",
    "Would support reductive branch if present",
    "Would support reductive branch if present",
    "Would support reductive branch if present",
    "Central carbon metabolism clearly active",
    "Central carbon metabolism clearly active",
    "Central carbon metabolism clearly active",
    "Central carbon metabolism clearly active",
    "Supports TCA backbone and reductive/oxidative flexibility",
    "Supports TCA backbone and reductive/oxidative flexibility",
    "Important because it links central carbon metabolism to biosynthesis",
    "Important because it links central carbon metabolism to biosynthesis"
  ),
  stringsAsFactors = FALSE
)

ferroplasma_rtca_markers



colors_pB<-c('#a6d854',"#fb8072",'#8da0cb','#66c2a5',"grey")
info_table_DNA <- read.csv("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/metagenomics/info_table_DNA.csv")
info_table_DNA$Alternative_ID2<-gsub("_","",info_table_DNA$Alternative_ID2)
data_proteomics_tax_prot<-merge(data_proteomics_prot,info_table_DNA,by.x = "tax_id",by.y = "Alternative_ID2")
library(dplyr)
library(tidyr)
library(stringr)

df_expanded <- data_proteomics_tax_prot %>%
  separate_rows(V12, sep = ",") %>%
  mutate(KO = stringr::str_remove(V12, "^ko:"))

df_joined <- df_expanded %>%
  left_join(ferroplasma_rtca_markers, by = "KO")
df_joined_filtered <- df_joined[-which(is.na(df_joined$importance)),]

library(dplyr)

# df_heat <- df_joined_filtered %>%
#   mutate(avg_LFQ = rowMeans(
#     cbind(log10_LFQ_intensity_Cave_Biofilm_1,
#           log10_LFQ_intensity_Cave_Biofilm_2),
#     na.rm = TRUE
#   ))

library(dplyr)
library(tidyr)
library(ggplot2)

# 1️⃣ Compute average LFQ (as before)
df_heat <- df_joined_filtered %>%
  rowwise() %>%
  mutate(
    avg_LFQ = mean(c(
      log10_LFQ_intensity_Cave_Biofilm_1,
      log10_LFQ_intensity_Cave_Biofilm_2
    ), na.rm = TRUE)
  ) %>%
  ungroup()
# 2️⃣ Create a combined label for plotting
df_heat <- df_heat %>%
  mutate(Genus_MAGID = paste(Genus, tax_id, sep = "-"),KO_role = paste(pathway_block," (",KO,")",sep = ""))

# 3️⃣ Create a status column for coloring
# 'measured' = LFQ exists, 'not_detected' = LFQ is NA
df_heat <- df_heat %>%
  mutate(
    status = ifelse(!is.na(avg_LFQ), "measured", "not_detected")
  )
df_heat <- df_heat %>%
  mutate(
    Genus_MAGID = case_when(
      Genus_MAGID == "Ferroplasma-KNPMNEEE" ~ "Ferroplasma c.",
      Genus_MAGID == "Ferroplasma-SCMAG00006" ~ "Ferroplasma MAG6",
      Genus_MAGID == "Mycobacterium-KDJLIKBO" ~ "M. Methanotrophicum",
      TRUE ~ Genus_MAGID
    )
  )
rTCA <- ggplot(df_heat, aes(x = KO_role, y = Genus_MAGID)) +
  # fill measured LFQ values
  geom_tile(aes(fill = avg_LFQ), color = NA) +
  
  # grey for KOs not detected
  geom_tile(
    data = subset(df_heat, status == "not_detected"),
    aes(x = KO_role, y = Genus_MAGID),
    fill = "grey80"
  ) +
  
  scale_fill_viridis_c(
    na.value = "grey80",
    name = "Avg log10 LFQ"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "KO",
    y = "Taxonomic ID",
    title = ""
  )
rTCA
# ggsave(
#   filename = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/rTCA.pdf",
#   plot = rTCA,
#   width = 12,        # in inches
#   height = 6,
#   units = "in",
#   dpi = 600,         # high resolution (even though vector)
#   device = cairo_pdf # ensures good font embedding
# )
# 
# # 2️⃣ Save as high-resolution PNG (for presentations)
# ggsave(
#   filename = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/rTCA.png",
#   plot = rTCA,
#   width = 12,
#   height = 6,
#   units = "in",
#   dpi = 600
# )
# #methane


cave_metabolism_markers_with_KO <- data.frame(
  metabolism = c(
    # methanotrophy
    "methanotrophy", "methanotrophy", "methanotrophy",
    
    # methanogenesis
    "methanogenesis", "methanogenesis", "methanogenesis",
    
    # heterotrophy carbohydrate
    "heterotrophy_carbohydrate", "heterotrophy_carbohydrate",
    
    # heterotrophy TCA / organic acids
    "heterotrophy_TCA_organic_acids", "heterotrophy_TCA_organic_acids",
    "heterotrophy_TCA_organic_acids", "heterotrophy_TCA_organic_acids",
    "heterotrophy_TCA_organic_acids", "heterotrophy_TCA_organic_acids",
    
    # heterotrophy acetate
    "heterotrophy_acetate"
  ),
  
  gene = c(
    # methanotrophy
    "pmoA", "pmoB", "pmoC",
    
    # methanogenesis
    "mcrA", "mcrB", "mcrG",
    
    # heterotrophy carbohydrate
    "glk", "pyk",
    
    # heterotrophy TCA / organic acids
    "acnA/acnB", "icd", "fum", "mdh", "sucC/sucD", "sucC/sucD",
    
    # heterotrophy acetate
    "acs"
  ),
  
  KO = c(
    # methanotrophy
    "K16157", "K16158", "K16159",
    
    # methanogenesis
    "K00399", "K00401", "K00402",
    
    # heterotrophy carbohydrate
    "K00845", "K00873",
    
    # heterotrophy TCA / organic acids
    "K01681", "K00031", "K01679", "K00024", "K01902", "K01903",
    
    # heterotrophy acetate
    "K01895"
  ),
  
  enzyme = c(
    # methanotrophy
    "Particulate methane monooxygenase subunit A",
    "Particulate methane monooxygenase subunit B",
    "Particulate methane monooxygenase subunit C",
    
    # methanogenesis
    "Methyl-coenzyme M reductase alpha subunit",
    "Methyl-coenzyme M reductase beta subunit",
    "Methyl-coenzyme M reductase gamma subunit",
    
    # heterotrophy carbohydrate
    "Glucokinase",
    "Pyruvate kinase",
    
    # heterotrophy TCA / organic acids
    "Aconitate hydratase / aconitase",
    "Isocitrate dehydrogenase",
    "Fumarate hydratase",
    "Malate dehydrogenase",
    "Succinyl-CoA synthetase beta subunit",
    "Succinyl-CoA synthetase alpha subunit",
    
    # heterotrophy acetate
    "Acetyl-CoA synthetase"
  ),
  
  pathway_role = c(
    # methanotrophy
    "Core methane oxidation marker",
    "Core methane oxidation marker",
    "Core methane oxidation marker",
    
    # methanogenesis
    "Core methane-producing marker",
    "Core methane-producing marker",
    "Core methane-producing marker",
    
    # heterotrophy carbohydrate
    "Sugar uptake/catabolism",
    "Glycolysis end point",
    
    # heterotrophy TCA / organic acids
    "TCA cycle",
    "TCA cycle",
    "TCA cycle",
    "TCA cycle",
    "TCA cycle",
    "TCA cycle",
    
    # heterotrophy acetate
    "Acetate assimilation"
  ),
  
  importance = c(
    # methanotrophy
    "diagnostic", "diagnostic", "diagnostic",
    
    # methanogenesis
    "diagnostic", "diagnostic", "diagnostic",
    
    # heterotrophy carbohydrate
    "supportive", "supportive",
    
    # heterotrophy TCA / organic acids
    "supportive", "supportive", "supportive", "supportive", "supportive", "supportive",
    
    # heterotrophy acetate
    "supportive"
  ),
  
  cave_context_interpretation = c(
    # methanotrophy
    "If present, indicates direct methane oxidation",
    "If present, indicates direct methane oxidation",
    "If present, indicates direct methane oxidation",
    
    # methanogenesis
    "If present, strongly indicates methane production",
    "If present, strongly indicates methane production",
    "If present, strongly indicates methane production",
    
    # heterotrophy carbohydrate
    "Supports sugar use if organics are available",
    "Supports glycolytic conversion to pyruvate",
    
    # heterotrophy TCA / organic acids
    "Central TCA backbone, also relevant to reductive variants",
    "Central TCA backbone, also relevant to reductive variants",
    "Central TCA backbone, also relevant to reductive variants",
    "Central TCA backbone, also relevant to reductive variants",
    "Central TCA backbone, also relevant to reductive variants",
    "Central TCA backbone, also relevant to reductive variants",
    
    # heterotrophy acetate
    "Supports acetate scavenging / assimilation"
  ),
  
  stringsAsFactors = FALSE
)
# Load packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("KEGGREST")
library(KEGGREST)
library(dplyr)
library(tidyr)

# 1️⃣ Get KOs in the methane metabolism pathway
methane_pathway <- keggGet("map00680")[[1]]  # "Methane metabolism"

# Extract the KOs listed in the PATHWAY entry
# They are usually under "ORTHOLOGY"
kegg_kos <- methane_pathway$ORTHOLOGY

# kegg_kos is a named list: names = KO IDs, values = descriptions
df_kos <- data.frame(
  KO = names(kegg_kos),
  Name = unname(kegg_kos),
  stringsAsFactors = FALSE
)

# 2️⃣ Get KEGG modules
# Find all modules that include any of these KOs
modules <- methane_pathway$MODULE
modules_ids <- names(methane_pathway$MODULE)

ko_module_list <- list()

for (mod_id in modules_ids) {
  mod <- keggGet(mod_id)[[1]]
  
  if (!is.null(mod$DEFINITION)) {
    # Extract KOs in module
    kos_in_mod <- names(mod$ORTHOLOGY)
    kos_in_mod_all <- mod$ORTHOLOGY
    
    if (length(kos_in_mod) > 0) {
      ko_module_list[[length(ko_module_list)+1]] <- data.frame(
        KO = kos_in_mod,
        KO_names = kos_in_mod_all,
        Module = mod_id,
        Module_name = mod$NAME,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Combine into a single table
df_ko_modules <- do.call(rbind, ko_module_list)
df_ko_expanded <- df_ko_modules %>%
  # unify separators
  mutate(KO = str_replace_all(KO, "\\+", ",")) %>%
  # separate multiple KOs into rows
  separate_rows(KO, sep = ",") %>%
  # trim whitespace
  mutate(KO = str_trim(KO)) %>%
  # keep only the KO ID
  mutate(KO = str_extract(KO, "^K[0-9]{5}")) %>%
  # keep only relevant columns (Module, Module_name)
  select(KO, Module, Module_name) %>%
  # group by KO and concatenate modules
  group_by(KO) %>%
  summarise(
    Modules = paste(unique(Module), collapse = ";"),
    Module_names = paste(unique(Module_name), collapse = ";"),
    .groups = "drop"
  )

# View result
df_ko_expanded
# 2️⃣ Result


#cave_metabolism_markers_with_KO
df_joined <- df_expanded %>%
  left_join(df_ko_expanded, by = "KO")
df_joined_filtered <- df_joined[-which(is.na(df_joined$Modules)),]

library(dplyr)
library(tidyr)
library(ggplot2)

# 1️⃣ Compute average LFQ (as before)
df_heat <- df_joined_filtered %>%
  rowwise() %>%
  mutate(
    avg_LFQ = mean(c(
      log10_LFQ_intensity_Cave_Biofilm_1,
      log10_LFQ_intensity_Cave_Biofilm_2
    ), na.rm = TRUE)
  ) %>%
  ungroup()
# 2️⃣ Create a combined label for plotting
df_heat <- df_heat %>%
  mutate(Genus_MAGID = paste(Genus, tax_id, sep = "-"),KO_role = paste(Modules," (",KO,")",sep = ""))

# 3️⃣ Create a status column for coloring
# 'measured' = LFQ exists, 'not_detected' = LFQ is NA
df_heat <- df_heat %>%
  mutate(
    status = ifelse(!is.na(avg_LFQ), "measured", "not_detected")
  )
df_heat <- df_heat %>%
  mutate(
    Genus_MAGID = case_when(
      Genus_MAGID == "Ferroplasma-KNPMNEEE" ~ "Ferroplasma c.",
      Genus_MAGID == "Ferroplasma-SCMAG00006" ~ "Ferroplasma MAG6",
      Genus_MAGID == "Mycobacterium-KDJLIKBO" ~ "M. Methanotrophicum",
      TRUE ~ Genus_MAGID
    )
  )
methane_general <- ggplot(df_heat, aes(x = KO_role, y = Genus_MAGID)) +
  # fill measured LFQ values
  geom_tile(aes(fill = avg_LFQ), color = NA) +
  
  # grey for KOs not detected
  geom_tile(
    data = subset(df_heat, status == "not_detected"),
    aes(x = KO_role, y = Genus_MAGID),
    fill = "grey80"
  ) +
  
  scale_fill_viridis_c(
    na.value = "grey80",
    name = "Avg log10 LFQ"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    x = "Modules (KO)",
    y = "Taxonomic ID",
    title = ""
  )
methane_general

#############################################
# Lingyi: Add categorical strip labels to methane_general
# Add a category column based on Module_names content

length(df_heat$KO_role)
length(df_heat_new$KO_role)
length(unique(df_heat_new$KO_role))

# How many unique KO_role per category?
df_heat_new %>%
  select(KO_role, category) %>%
  distinct() %>%
  count(category)

df_heat_new <- df_heat %>%
  mutate(
    category = case_when(
      grepl("methanotroph|methane oxidation", Module_names, ignore.case = TRUE) ~ "Methanotroph related",
      grepl("methanogenesis|methane production", Module_names, ignore.case = TRUE) ~ "Methanogenesis related",
      TRUE ~ "Other"
    )
  )

methane_general <- ggplot(df_heat_new, aes(x = KO_role, y = Genus_MAGID)) +
  geom_tile(aes(fill = avg_LFQ), color = NA) +
  geom_tile(
    data = subset(df_heat, status == "not_detected"),
    aes(x = KO_role, y = Genus_MAGID),
    fill = "grey80"
  ) +
  scale_fill_viridis_c(na.value = "grey80", name = "Avg log10 LFQ") +
  # facet groups KOs by category with label strip at the bottom
  facet_grid(. ~ category, scales = "free_x", space = "free_x",
             switch = "x") +          # <-- moves strip to bottom
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y      = element_text(size = 8),
    strip.text.x     = element_text(size = 10, face = "bold"),
    strip.placement  = "outside",     # <-- strip outside axis labels
    panel.spacing    = unit(0.3, "lines")
  ) +
  labs(x = NULL, y = "Taxonomic ID", title = "")

methane_general

# ggsave(
#   filename = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/methane_general.pdf",
#   plot = methane_general,
#   width = 12,        # in inches
#   height = 6,
#   units = "in",
#   dpi = 600,         # high resolution (even though vector)
#   device = cairo_pdf # ensures good font embedding
# )
# 
# # 2️⃣ Save as high-resolution PNG (for presentations)
# ggsave(
#   filename = "/home/chrats/Desktop/Projects/Mycobacterium_sulfur_cave/FIGURES/methane_general.png",
#   plot = methane_general,
#   width = 12,
#   height = 6,
#   units = "in",
#   dpi = 600
# )
