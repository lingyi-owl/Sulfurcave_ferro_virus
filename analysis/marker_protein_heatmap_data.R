# Shared data preparation for Fig. S4 (FigS04_rTCA_heatmap.R) and Fig. S5
# (FigS05_methane_heatmap.R). Builds the marker-protein table `df_cell` and related
# objects. Formerly plotting/analysis/marker_protein_heatmap_data.R; its ggplot2 draft
# plot (figures/supplementary/FigS04/drafts/) was removed.
#
# =============================================================================
# rTCA / TCA marker-protein heatmap  (supplementary version of Figure03_5.R)
#
# Reproduces the reductive-TCA proteomics heatmap, but with the rows relabelled
# to "<genus> MAG <n>" (MAG number from the study bin ID) instead of the
# original "<genus>-<tax_id>" labels.
#
# Cell convention (identical to the original Figure03_5.R):
#   colour (viridis) = marker protein DETECTED & quantified in the cave biofilm
#                      metaproteome  (average log10 LFQ over the two replicates)
#   grey             = marker KO gene PRESENT in the genome but protein NOT
#                      detected proteomically
#   white (no tile)  = genome does not carry that marker KO at all
#
# Inputs (relative to the git repo root, set `repo` below):
#   proteomics/proteinGroups.xlsx                                  (metaproteomics LFQ)
#   metagenomics/Functional_annotation/sulfur_cave_ref_proteins_eggnog.emapper.annotations
#   metagenomics/info_table_DNA.csv                                (tax_id -> Genus)
#
# Row-label provenance:
#   - MAG number  : trailing integer of the study bin ID (SC_MAG_000NN / SC_Bin_000NN)
#   - genus       : Table 2 (Taxonomy assignment of MAGs); blank -> "NA"
#   - the 5 in-house reference genomes (non-SC tax_ids) are shown as
#     "<original figure label> / <genus> MAG <n>", the bin they map to per
#     proteomics/info_table_cleanids.csv
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(ggtext)
})

repo <- REPO

# ---- 1. Metaproteomics LFQ (proteinGroups sheet 2) --------------------------
ProteinGroups <- read_xlsx(file.path(repo, "proteomics/proteinGroups.xlsx"), sheet = 2)
data_proteomics <- as.data.frame(ProteinGroups[, 1:2])
data_proteomics$LFQ_intensity_Cave_Biofilm_1 <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_1`)
data_proteomics$LFQ_intensity_Cave_Biofilm_2 <- as.numeric(ProteinGroups$`LFQ intensity cave_wcl_2`)
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_1 <- log10(data_proteomics$LFQ_intensity_Cave_Biofilm_1)
data_proteomics$log10_LFQ_intensity_Cave_Biofilm_2 <- log10(data_proteomics$LFQ_intensity_Cave_Biofilm_2)
data_proteomics[data_proteomics == -Inf] <- 0
data_proteomics$ID <- gsub("_", "", data_proteomics$`Majority protein IDs`)

# ---- 2. eggNOG annotations of the reference proteome ------------------------
egg <- read.delim(file.path(repo, "metagenomics/Functional_annotation/sulfur_cave_ref_proteins_eggnog.emapper.annotations"),
                  header = FALSE, comment.char = "#")
# tax_id: SC_MAG_00014 -> SCMAG00014 ; Crude_2 -> Crude2 ; otherwise first token
egg$tax_id <- sapply(strsplit(egg$V1, "_"), function(x) {
  if (length(x) >= 3 && x[1] == "SC" && x[2] %in% c("MAG", "Bin")) paste0(x[1:3], collapse = "")
  else if (length(x) >= 2 && x[1] == "Crude")                     paste0(x[1:2], collapse = "")
  else                                                            x[1]
})
egg$ID <- gsub("_", "", egg$V1)

data_proteomics_prot <- data_proteomics %>%
  separate_rows(ID, sep = ";") %>%
  right_join(egg, by = "ID")

# ---- 3. rTCA / TCA marker KOs (24) ------------------------------------------
ferroplasma_rtca_markers <- data.frame(
  KO = c("K00169","K00170","K00171","K00172",      # porA-D          (rTCA core)
         "K00174","K00175","K00176","K00177",      # korA-D/oorA-D   (rTCA core)
         "K15230","K15231",                        # aclA/B          (citrate cleavage)
         "K15232","K15233",                        # ccs/ccl         (alt. citrate cleavage)
         "K00244","K00245","K00246","K00247",      # frdA-D          (reductive succinate branch)
         "K00031","K01681","K01676","K00024","K01902","K01903",  # backbone
         "K01006","K01007"),                       # ppsA/ppdk       (gluconeogenesis)
  pathway_block = c(rep("rTCA core ", 8), rep("Citrate cleavage", 2),
                    rep("Alternative citrate cleavage", 2),
                    rep("Reductive succinate branch", 4),
                    rep("rTCA backbone", 6), rep("Gluconeogenesis", 2)),
  importance = c(rep("core", 8), rep("diagnostic", 4), rep("supportive", 10), rep("important", 2)),
  stringsAsFactors = FALSE
)

# ---- 4. attach genus (info_table_DNA) & keep only marker KOs ----------------
info_table_DNA <- read.csv(file.path(repo, "metagenomics/info_table_DNA.csv"), sep = ";")
info_table_DNA$Alternative_ID2 <- gsub("_", "", info_table_DNA$Alternative_ID2)
data_proteomics_tax_prot <- merge(data_proteomics_prot, info_table_DNA,
                                  by.x = "tax_id", by.y = "Alternative_ID2")

df_expanded <- data_proteomics_tax_prot %>%
  separate_rows(V12, sep = ",") %>%
  mutate(KO = str_remove(V12, "^ko:"))
df_joined <- df_expanded %>% left_join(ferroplasma_rtca_markers, by = "KO")
df_joined_filtered <- df_joined[-which(is.na(df_joined$importance)), ]   # keep marker KOs only

# average LFQ over the two biofilm replicates + KO column label
df_heat <- df_joined_filtered %>%
  rowwise() %>%
  mutate(avg_LFQ = mean(c(log10_LFQ_intensity_Cave_Biofilm_1,
                          log10_LFQ_intensity_Cave_Biofilm_2), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(KO_role = paste0(pathway_block, " (", KO, ")"))
# an all-NA average comes back as NaN -> mark as not-detected
df_heat$avg_LFQ[is.nan(df_heat$avg_LFQ)] <- NA

# one cell per (tax_id, KO_role): mean of measured values (na.rm)
df_cell <- df_heat %>%
  group_by(tax_id, KO_role) %>%
  summarise(avg_LFQ = ifelse(all(is.na(avg_LFQ)), NA_real_, mean(avg_LFQ, na.rm = TRUE)),
            .groups = "drop") %>%
  mutate(status = ifelse(is.na(avg_LFQ), "not_detected", "measured"))

# =============================================================================
# 5. ROW LABELS :  <genus> MAG <n>   (see header for provenance)
# =============================================================================
# study bin ID for each tax_id (SC ids map to themselves; references via cleanids)
tax2bin <- c(
  SCMAG00001="SC_MAG_00001", SCMAG00002="SC_MAG_00002", SCMAG00003="SC_MAG_00003",
  SCMAG00004="SC_MAG_00004", SCMAG00005="SC_MAG_00005", SCMAG00006="SC_MAG_00006",
  SCMAG00007="SC_MAG_00007", SCMAG00009="SC_MAG_00009", SCMAG00010="SC_MAG_00010",
  SCMAG00011="SC_MAG_00011", SCMAG00012="SC_MAG_00012", SCMAG00013="SC_MAG_00013",
  SCMAG00014="SC_MAG_00014", SCMAG00015="SC_MAG_00015", SCMAG00017="SC_MAG_00017",
  SCBin00023="SC_Bin_00023", SCBin00024="SC_Bin_00024", SCBin00025="SC_Bin_00025",
  SCBin00027="SC_Bin_00027", SCBin00028="SC_Bin_00028", SCBin00029="SC_Bin_00029",
  SCBin00030="SC_Bin_00030", SCBin00032="SC_Bin_00032", SCBin00033="SC_Bin_00033",
  KDJLIKBO ="SC_MAG_00016", KNPMNEEE ="SC_MAG_00008", EKPBOPOH ="SC_MAG_00018",
  Crude2   ="SC_MAG_00019", Crude3   ="SC_Bin_00030")   # Crude3 = SC_Bin_00022/SC_Bin_00030 -> MAG 30

mag_num <- function(tid) {
  b <- tax2bin[[tid]]
  as.integer(sub("^0+", "", sub(".*_", "", b)))
}

# genus per MAG number, from Table 2 (blank = unassigned -> "NA")
genus_by_mag <- c(
  `1`="Paracoccus", `2`="Lawsonella", `3`="", `4`="Acidithiobacillus", `5`="Cutibacterium",
  `6`="Ferroplasma", `7`="PALSA-457", `8`="Ferroplasma", `9`="Acidiphilium", `10`="Ferrimicrobium",
  `11`="", `12`="Palsa-601", `13`="Acidithiobacillus", `14`="Sulfobacillus", `15`="Micrococcus",
  `16`="Mycobacterium", `17`="Corynebacterium", `18`="Cuniculiplasma", `19`="Mycobacterium",
  `23`="Moraxella", `24`="B-DKE", `25`="Micrarchaeum", `27`="Janibacter",
  `28`="Staphylococcus", `29`="Streptococcus", `30`="Mycobacterium", `32`="Micrarchaeum", `33`="")
italic_genera <- setdiff(unique(unlist(genus_by_mag)), "")   # italicise every genus name

# base markdown label ( *genus* italic where Table 2 italicises it )
base_label <- function(tid) {
  n <- mag_num(tid); g <- genus_by_mag[[as.character(n)]]
  if (is.null(g) || g == "") return(paste0("NA MAG ", n))
  if (g %in% italic_genera)  return(paste0("*", g, "* MAG ", n))
  paste0(g, " MAG ", n)
}

# per-tax_id markdown labels (ggtext markdown: *...* = italic)
label_md <- setNames(vapply(df_cell$tax_id %>% unique(), base_label, character(1)),
                     df_cell$tax_id %>% unique())
# dual labels for the 5 in-house reference genomes:  <original figure label> / <new label>
label_md["KDJLIKBO"] <- "*Ca*. M. Methanotrophicum/*Mycobacterium* MAG 16"
label_md["EKPBOPOH"] <- "*Cuniculiplasma*-EKPBOPOH/*Cuniculiplasma* MAG 18"
label_md["KNPMNEEE"] <- "*Ferroplasma* MAG 8"
label_md["Crude2"]   <- "*Mycobacterium* bin 2/*Mycobacterium* MAG 19"
label_md["Crude3"]   <- "*Mycobacterium* bin 3/*Mycobacterium* MAG 30"

df_cell$row_label <- label_md[df_cell$tax_id]

# row order: ascending MAG number; the Crude3 draft sits just after the SC_Bin_00030 bin
df_cell$order_key <- vapply(df_cell$tax_id,
                            function(t) mag_num(t) + ifelse(t == "Crude3", 0.5, 0), numeric(1))
row_levels <- df_cell %>% distinct(row_label, order_key) %>% arrange(order_key) %>% pull(row_label)
df_cell$row_label <- factor(df_cell$row_label, levels = rev(row_levels))  # ggplot: first level at bottom

# column order: pathway-block order, then KO, keeping only KOs present in the data
col_levels <- ferroplasma_rtca_markers %>%
  mutate(KO_role = paste0(pathway_block, " (", KO, ")")) %>%
  filter(KO_role %in% df_cell$KO_role) %>% pull(KO_role) %>% unique()
df_cell$KO_role <- factor(df_cell$KO_role, levels = col_levels)
