#!/usr/bin/env Rscript
# =============================================================================
# Proteome allocation to translation / growth machinery — PER SAMPLE
# Two biofilm whole-cell-lysate samples + one in-vitro culture, shown separately.
# Lineage = bar colour; sample = fill pattern.
#
# Culture bars are flagged unreliable by ROBUST identification count in that
# sample (identified "By MS/MS" AND >=1 unique peptide), not by a raw protein
# count: SC_MAG_00006 culture = 16 robust IDs (weak, hatched);
#         SC_MAG_00008 culture = 50 robust IDs (solid, kept);
#         Mycobacterium culture = well supported (solid).
#
# Inputs (edit the CONFIG block):
#   proteinGroups.xlsx, sulfur_cave_ref_proteins_eggnog.emapper.annotations,
#   ko_category_map.csv
# Outputs:
#   allocation_by_sample.csv, allocation_split_by_sample.png / .pdf
# =============================================================================
devtools::install_github("IMNMV/ClaudeR")
library(ClaudeR)
install_cli(tools = "claude")
claudeAddin()   # starts the bridge; keep RStudio open

suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(ggplot2); library(ggpattern)
  library(ggtext)
})

# ---- CONFIG -----------------------------------------------------------------
SULFURCAVE <- "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git"
PROT_XLSX <- file.path(SULFURCAVE, "proteomics", "proteinGroups.xlsx")
EGGNOG    <- file.path(SULFURCAVE, "metagenomics", "Functional_annotation",
                       "sulfur_cave_ref_proteins_eggnog.emapper.annotations")
# KO -> category lookup (see proteomics/ko_map_provenance.txt for how it was built)
KO_MAP    <- file.path(SULFURCAVE, "proteomics", "ko_category_map.csv")

# three samples: two biofilm whole-cell-lysate + one culture (sup fractions dropped)
samples <- c(Biofilm_1 = "cave_wcl_1",
             Biofilm_2 = "cave_wcl_2",
             Culture   = "invitro")
lfq_col <- function(s) paste0("LFQ intensity ", s)
uniq_col<- function(s) paste0("Unique peptides ", s)
idt_col <- function(s) paste0("Identification type ", s)

cats6 <- c("Ribosomal proteins", "Translation factors",
           "Aminoacyl-tRNA synthetases", "Ribosome biogenesis",
           "DNA replication", "Cell division")
priority <- c("Ribosomal proteins", "Aminoacyl-tRNA synthetases",
              "Translation factors", "Ribosome biogenesis",
              "DNA replication", "Cell division")
ROBUST_MIN <- 30   # a lineage x sample is "reliable" if >= this many robust IDs

# ---- 1. Load + QC-filter ----------------------------------------------------
pg <- read_excel(PROT_XLSX, sheet = "ProteinGroups") %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`),
         is.na(`Only identified by site`),
         !str_starts(`Protein IDs`, "CON__"))

# ---- 2. eggNOG ORF -> KO(s) -------------------------------------------------
egg <- read_tsv(EGGNOG, comment = "#", col_names = FALSE,
                show_col_types = FALSE, progress = FALSE)
orf_ko_list <- egg %>%
  transmute(orf = X1, ko = str_extract_all(X12, "K\\d{5}")) %>%
  filter(lengths(ko) > 0) %>% { setNames(.$ko, .$orf) }
# ORF -> COG functional category (col 7); category V = defense mechanisms,
# CRISPR-Cas dominated. Used for the CRISPR-Cas/defense allocation overlay.
orf2cog <- setNames(egg$X7, egg$X1)

# ---- 3. KO -> category ------------------------------------------------------
ko_map <- read_csv(KO_MAP, show_col_types = FALSE)
ko2cat <- setNames(ko_map$category, ko_map$KO)

# ---- 4. Normalize IDs (proteinGroups compact -> eggNOG underscored) ---------
norm_id <- function(tok) {
  tok <- str_replace(tok, "^SCMAG(\\d+)_(c_.*)$", "SC_MAG_\\1_\\2")
  tok <- str_replace(tok, "^SCBin(\\d+)_(c_.*)$", "SC_Bin_\\1_\\2")
  tok <- str_replace(tok, "^Crude(\\d+)_(\\d+)$", "Crude_\\1_\\2")
  tok
}

# ---- 5. Per-group KOs, single category, lineage -----------------------------
group_kos <- function(majID) {
  toks <- norm_id(str_split(majID, ";")[[1]])
  unique(unlist(orf_ko_list[toks[toks %in% names(orf_ko_list)]]))
}
classify <- function(kos) {
  if (length(kos) == 0) return(NA_character_)
  hits <- unique(ko2cat[kos[kos %in% names(ko2cat)]])
  for (p in priority) if (p %in% hits) return(p)
  NA_character_
}
prefixes_of <- function(majID)
  str_extract(str_split(majID, ";")[[1]], "^(SCMAG\\d+|SCBin\\d+|Crude\\d+|[A-Z]+)")
FERRO6 <- "SCMAG00006"; FERRO8 <- "KNPMNEEE"
# Mycobacterium = the complete "Mycobacterium Methanotrophicum" genome only
# (KDJLIKBO = SC_MAG_00016); the other conspecific bins (Crude2/Crude3/SCBin00030)
# are excluded at the user's request.
MYCO   <- c("KDJLIKBO")
assign_lineage <- function(majID) {
  p <- prefixes_of(majID)
  f6 <- FERRO6 %in% p; f8 <- FERRO8 %in% p; my <- any(MYCO %in% p)
  if ((f6 | f8) & my) return(NA_character_)
  if (f6 & f8) return(NA_character_)
  if (f6) return("SC_MAG_00006"); if (f8) return("SC_MAG_00008")
  if (my) return("Mycobacterium"); NA_character_
}
# defense (CRISPR-Cas): group carries an ORF in COG category V
group_defense <- function(majID) {
  toks <- norm_id(str_split(majID, ";")[[1]])
  any(str_detect(orf2cog[toks], "V"), na.rm = TRUE)
}
pg <- pg %>%
  mutate(category = map_chr(map(`Majority protein IDs`, group_kos), classify),
         lineage  = map_chr(`Majority protein IDs`, assign_lineage),
         defense  = map_lgl(`Majority protein IDs`, group_defense))

# ---- 6. Per-sample allocation + robust-ID count -----------------------------
lin_levels <- c("SC_MAG_00006", "SC_MAG_00008", "Mycobacterium")
alloc_one <- function(sname) {
  col <- lfq_col(samples[[sname]]); uc <- uniq_col(samples[[sname]]); ic <- idt_col(samples[[sname]])
  # Use ONLY robustly-sequenced groups: identified By MS/MS in this sample AND
  # carrying >=1 unique peptide there (and with a positive LFQ value). Both the
  # category numerator and the lineage-total denominator are built from this set.
  all_det <- pg %>% filter(!is.na(lineage), .data[[col]] > 0)
  d <- all_det %>% filter(.data[[ic]] == "By MS/MS", .data[[uc]] > 0)
  # allocation fractions (robust set only)
  frac <- d %>% group_by(lineage) %>% mutate(tot = sum(.data[[col]])) %>%
    group_by(lineage, category) %>%
    summarise(frac = sum(.data[[col]]) / dplyr::first(tot), .groups = "drop") %>%
    filter(category %in% cats6)
  # CRISPR-Cas / defense allocation (COG V): summed LFQ of defense groups over the
  # same robust lineage-total denominator. Added as a pseudo-category row so it
  # flows through pivot_wider alongside the six KO categories.
  defense <- d %>% group_by(lineage) %>%
    summarise(frac = sum(.data[[col]][defense]) / sum(.data[[col]]), .groups = "drop") %>%
    mutate(category = "Defense (CRISPR-Cas)")
  # counts: n_detected (any LFQ>0) for reference, n_robust (the set used above)
  counts <- all_det %>%
    mutate(rob = (.data[[ic]] == "By MS/MS") & (.data[[uc]] > 0)) %>%
    group_by(lineage) %>%
    summarise(n_detected = n(), n_robust = sum(rob, na.rm = TRUE), .groups = "drop")
  bind_rows(frac, defense) %>% left_join(counts, by = "lineage") %>% mutate(sample = sname)
}
by_cat <- map_dfr(names(samples), alloc_one)

summary_tbl <- by_cat %>%
  pivot_wider(names_from = category, values_from = frac, values_fill = 0)
for (c in cats6) if (!c %in% names(summary_tbl)) summary_tbl[[c]] <- 0
summary_tbl <- summary_tbl %>%
  mutate(ribosomal = `Ribosomal proteins`,
         translation_assoc = `Ribosomal proteins` + `Translation factors` +
           `Aminoacyl-tRNA synthetases` + `Ribosome biogenesis`,
         growth_machinery  = translation_assoc + `DNA replication` + `Cell division`,
         defense_crispr    = `Defense (CRISPR-Cas)`,
         low_conf = n_robust < ROBUST_MIN)
write_csv(summary_tbl, "allocation_by_sample.csv")

# ---- 7. Figure --------------------------------------------------------------
pal   <- c(SC_MAG_00006 = "#fc8d62", SC_MAG_00008 = "#ffd92f", Mycobacterium = "#66c2a5")
# genus italicised (scientific convention); strain/MAG suffixes stay upright.
# plotmath expressions: genus italicised, strain/MAG suffix upright. This is
# device-independent (does not depend on the graphics backend resolving an
# italic font face, unlike ggtext markdown). Positional list; order matches
# lin_levels = c(SC_MAG_00006, SC_MAG_00008, Mycobacterium).
leg   <- list(
  expression(italic("Ferroplasma") * " MAG 6"),
  expression(italic("Ferroplasma") * " c. MAG 8"),
  expression(italic("Mycobacterium"))
)
# ggplot has no per-bar hatching; encode low-confidence with an outline + alpha instead.
mk_df <- function(metrics, labels) {
  # species-major bar order: MAG6{B1,B2,Cult}, MAG8{B1,B2,Cult}, Myco{B1,B2,Cult}
  bar_levels <- as.vector(t(outer(lin_levels, names(samples), paste, sep = "|")))
  summary_tbl %>%
    select(lineage, sample, low_conf, all_of(metrics)) %>%
    filter(!low_conf) %>%   # drop low-confidence bars entirely (kept in the CSV)
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "frac") %>%
    mutate(pct = 100 * frac,
           lineage = factor(lineage, levels = lin_levels),
           sample  = factor(sample, levels = names(samples)),
           bar     = factor(paste(lineage, sample, sep = "|"), levels = bar_levels),
           metric  = factor(metric, levels = metrics, labels = labels),
           # Culture = solid species-coloured fill; biofilms = white (stripes carry colour)
           fill_col = ifelse(sample == "Culture", pal[as.character(lineage)], "white")) %>%
    arrange(metric, bar)
}
cat_labs <- c("Ribosomal\nproteins","Translation\nfactors","aaRS",
              "Ribosome\nbiogenesis","DNA\nreplication","Cell\ndivision")
dfA <- mk_df(cats6, cat_labs)

# Sample -> pattern: Biofilm 1 = single diagonal stripe, Biofilm 2 = crosshatch
# (both crossing stripes), Culture = empty (hollow / no fill). Species stays the
# fill colour.
pat_map    <- c(Biofilm_1 = "stripe", Biofilm_2 = "crosshatch", Culture = "none")
pat_angle  <- c(Biofilm_1 = 45,       Biofilm_2 = 45,           Culture = 0)
pat_labs   <- c(Biofilm_1 = "Biofilm 1",
                Biofilm_2 = "Biofilm 2",
                Culture   = "Culture")

# Allocation to growth-machinery categories.
# White bar interior throughout; species identity is carried by the OUTLINE and
# the STRIPE colour so the Culture bar (no pattern) reads as an empty bar with a
# coloured border. Biofilm 1 = single stripe (45 deg), Biofilm 2 = crosshatch
# (both crossing stripes), Culture = empty.
pA <- ggplot(dfA, aes(metric, pct, group = bar)) +
  geom_col_pattern(
    aes(colour = lineage,                    # species outline colour
        fill = fill_col,                     # Culture solid; biofilms white
        pattern = sample,                    # stripe / crosshatch / none
        pattern_angle = sample,
        pattern_fill = lineage,              # stripes coloured by species
        pattern_colour = lineage,
        linetype = low_conf),                # dotted outline = unreliable culture
    position        = position_dodge2(width = 0.9, padding = 0.15,
                                       preserve = "single"),
    linewidth       = 0.4,
    pattern_density = 0.25, pattern_spacing = 0.030,
    pattern_key_scale_factor = 0.5) +
  scale_colour_manual(values = pal, labels = leg, name = "Lineage") +
  scale_fill_identity(guide = "none") +
  scale_pattern_manual(values = pat_map, labels = pat_labs, name = "Sample") +
  scale_pattern_angle_manual(values = pat_angle, guide = "none") +
  scale_pattern_fill_manual(values = pal, guide = "none") +
  scale_pattern_colour_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dotted"),
                        guide = "none") +
  labs(x = NULL, y = "% of detected proteome (LFQ)", title = NULL) +
  theme_minimal(base_size = 20, base_family = "Arial") +
  theme(text = element_text(family = "Arial", size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
        axis.text.y = element_text(angle = 0, size = 20),
        axis.title  = element_text(size = 20),
        
        legend.position = "right",
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        
        # Gap between individual legend keys
        legend.key.spacing.y = unit(0.15, "cm"),
        
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
        axis.line.x  = element_line(colour = "black", linewidth = 0.6),
        axis.line.y  = element_line(colour = "black", linewidth = 0.6)) +
  guides(
    # species legend: show solid colour swatches
    colour  = guide_legend(override.aes = list(fill = unname(pal[lin_levels]),
                                               pattern = "none")),
    # sample legend: neutral grey; Culture key shows a solid grey fill
    pattern = guide_legend(override.aes = list(colour = "grey40",
                                               fill = c("white", "white", "grey40"),
                                               pattern_fill = "grey40",
                                               pattern_colour = "grey40")))

# Render to the active graphics device (RStudio Plots panel) instead of saving.
print(pA)

# print(summary_tbl %>% select(lineage, sample, n_detected, n_robust, low_conf,
#                              ribosomal, growth_machinery))
# message("Done.")
