#!/usr/bin/env Rscript
# =============================================================================
# Fig. S6 — Proteome allocation to growth-machinery categories + CRISPR-Cas
# (per sample). The six growth categories are mutually exclusive (one KO priority
# class per protein group); CRISPR-Cas is a separate metric (DefenseFinder Cas
# loci, matching the Figure 4 definition) shown as an additional bar, and may
# overlap the six. Note: for SC_MAG_00008 a few very highly expressed Cas
# proteins give a large summed allocation, so the CRISPR-Cas bar can exceed the
# growth-machinery bars.
# Two biofilm whole-cell-lysate samples + one lab CH4 (in-vitro) culture, shown
# separately. Lineage = bar colour/outline; sample = fill pattern.
#
# Only robustly-sequenced protein groups are used (identified "By MS/MS" AND
# >=1 unique peptide in that sample), for both numerator and denominator; a
# lineage x sample with < ROBUST_MIN robust IDs is dropped from the figure
# (SC_MAG_00006 culture = 16 robust IDs -> dropped; kept in the CSV).
#
# Inputs (CONFIG block): proteinGroups.xlsx, eggNOG annotations, ko_category_map.csv
# Output: figures/supplementary/FigS06/FigS06_proteome_allocation.svg (Fig. S6)
# =============================================================================
# Optional library path for packages installed in a project-local library.
extra_lib <- Sys.getenv("SULFURCAVE_R_LIB", "")
if (nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(purrr); library(ggplot2); library(ggpattern)
})

# ---- CONFIG: Rscript FigS06_proteome_allocation.R [project_root] [output_dir]
argv <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else getwd()
SULFURCAVE <- if (length(argv) >= 1) normalizePath(argv[1]) else dirname(script_dir)
OUTDIR <- if (length(argv) >= 2) argv[2] else file.path(SULFURCAVE, "GrowthRate", "reproduced_allocation_R")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
PROT_XLSX <- file.path(SULFURCAVE, "proteomics", "proteinGroups.xlsx")
EGGNOG <- file.path(SULFURCAVE, "metagenomics", "Functional_annotation", "sulfur_cave_ref_proteins_eggnog.emapper.annotations")
KO_MAP <- file.path(SULFURCAVE, "proteomics", "ko_category_map.csv")
OUTFILE <- if (length(argv) >= 2) file.path(OUTDIR, "FigS06_proteome_allocation.svg") else
  file.path(SULFURCAVE, "figures", "supplementary", "FigS06", "FigS06_proteome_allocation.svg")  # Fig. S6

FONT <- "Arial"; FS <- 20   # all figure text: Arial, size 20, black

# three samples: two biofilm whole-cell-lysate + one lab CH4 culture
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
ROBUST_MIN <- 50   # a lineage x sample is "reliable" if >= this many robust IDs
# Set to 50 to match the reported methods. The figure is identical for any value
# in 17..50: the robust-ID counts are 16, 50, 294, 333, 465, 599, 704, 847, 951,
# so only SC_MAG_00006 in the culture (16) is ever dropped. Note SC_MAG_00008 in
# the culture sits exactly at 50 and is therefore KEPT (test is n_robust < MIN);
# use 51 if that bar should also be dropped.

# ---- 1. Load + QC-filter ----------------------------------------------------
pg <- read_excel(PROT_XLSX, sheet = "ProteinGroups", guess_max = 100000L) %>%
  mutate(excel_row = row_number() + 1L) %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`),
         is.na(`Only identified by site`),
         !str_starts(`Protein IDs`, "CON__"))

# ---- 2. eggNOG ORF -> KO(s) -------------------------------------------------
egg <- read_tsv(EGGNOG, comment = "#", col_names = FALSE,
                show_col_types = FALSE, progress = FALSE)
orf_ko_list <- egg %>%
  transmute(orf = X1, ko = str_extract_all(X12, "K\\d{5}")) %>%
  filter(lengths(ko) > 0) %>% { setNames(.$ko, .$orf) }
# CRISPR-Cas proteins from DefenseFinder (subtype starting with "CAS"), matching
# the manuscript's Figure 4 definition. eggNOG COG categories do NOT reliably
# classify Cas proteins, so a COG-based proxy misses them; the Cas loci are read
# from the per-genome DefenseFinder system tables and matched by locus tag.
DF_DIR <- file.path(SULFURCAVE, "metagenomics", "defencefinder")
read_cas <- function(d) {
  f <- file.path(DF_DIR, d, "defense_finder_systems.tsv")
  if (!file.exists(f)) stop("Missing DefenseFinder input: ", f)
  s <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
  unlist(strsplit(s$protein_in_syst[grepl("^CAS", s$subtype)], ","))
}
cas_raw <- unique(unlist(lapply(c("Ferro", "MAG6", "Myco"), read_cas)))
# proteinGroups uses bare KNPMNEEE_ tags (SC_MAG_00008) and SCMAG00006_-prefixed
# tags (SC_MAG_00006); cover both token forms.
cas_set <- unique(c(cas_raw, paste0("SCMAG00006_", cas_raw[startsWith(cas_raw, "c_")])))

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
# (KDJLIKBO = SC_MAG_00016).
MYCO   <- c("KDJLIKBO")
assign_lineage <- function(majID) {
  p <- prefixes_of(majID)
  f6 <- FERRO6 %in% p; f8 <- FERRO8 %in% p; my <- any(MYCO %in% p)
  if ((f6 | f8) & my) return(NA_character_)
  if (f6 & f8) return(NA_character_)
  if (f6) return("SC_MAG_00006"); if (f8) return("SC_MAG_00008")
  if (my) return("Mycobacterium"); NA_character_
}
# CRISPR-Cas flag: protein group carries a DefenseFinder Cas locus tag
group_cas <- function(majID) any(str_split(majID, ";")[[1]] %in% cas_set)
pg <- pg %>%
  mutate(category = map_chr(map(`Majority protein IDs`, group_kos), classify),
         lineage  = map_chr(`Majority protein IDs`, assign_lineage),
         cas      = map_lgl(`Majority protein IDs`, group_cas))

# ---- 6. Per-sample allocation + robust-ID count -----------------------------
lin_levels <- c("SC_MAG_00006", "SC_MAG_00008", "Mycobacterium")
alloc_one <- function(sname) {
  col <- lfq_col(samples[[sname]]); uc <- uniq_col(samples[[sname]]); ic <- idt_col(samples[[sname]])
  all_det <- pg %>% filter(!is.na(lineage), .data[[col]] > 0)
  d <- all_det %>% filter(.data[[ic]] == "By MS/MS", .data[[uc]] > 0)
  frac <- d %>% group_by(lineage) %>% mutate(tot = sum(.data[[col]])) %>%
    group_by(lineage, category) %>%
    summarise(frac = sum(.data[[col]]) / dplyr::first(tot), .groups = "drop") %>%
    filter(category %in% cats6)
  # CRISPR-Cas allocation (DefenseFinder Cas proteins): summed LFQ of Cas groups
  # over the same robust lineage-total denominator. Independent of the six
  # mutually-exclusive KO categories, so it is shown as a separate bar.
  cas_row <- d %>% group_by(lineage) %>%
    summarise(frac = sum(.data[[col]][cas]) / sum(.data[[col]]), .groups = "drop") %>%
    mutate(category = "CRISPR-Cas")
  counts <- all_det %>%
    mutate(rob = (.data[[ic]] == "By MS/MS") & (.data[[uc]] > 0)) %>%
    group_by(lineage) %>%
    summarise(n_detected = n(), n_robust = sum(rob, na.rm = TRUE), .groups = "drop")
  bind_rows(frac, cas_row) %>% left_join(counts, by = "lineage") %>% mutate(sample = sname)
}
by_cat <- map_dfr(names(samples), alloc_one)

summary_tbl <- by_cat %>%
  pivot_wider(names_from = category, values_from = frac, values_fill = 0)
for (c in c(cats6, "CRISPR-Cas"))
  if (!c %in% names(summary_tbl)) summary_tbl[[c]] <- 0
summary_tbl <- summary_tbl %>% mutate(low_conf = n_robust < ROBUST_MIN)

# ---- Reproducible tables and verification ------------------------------------
summary_tbl <- summary_tbl %>% mutate(
  total_growth_fraction = rowSums(pick(all_of(cats6))),
  ribosomal_pct = 100 * `Ribosomal proteins`,
  growth_machinery_pct = 100 * total_growth_fraction,
  CRISPR_Cas_pct = 100 * `CRISPR-Cas`, included = !low_conf)
write_csv(summary_tbl, file.path(OUTDIR, "allocation_summary.csv"))
write_csv(by_cat, file.path(OUTDIR, "allocation_by_category.csv"))
audit <- map_dfr(names(samples), function(sname) {
  ss <- samples[[sname]]
  pg %>% transmute(excel_row, protein_group_id = id,
    majority_protein_ids = `Majority protein IDs`, lineage, category, cas,
    sample = sname, LFQ = .data[[lfq_col(ss)]],
    unique_peptides = .data[[uniq_col(ss)]],
    identification_type = .data[[idt_col(ss)]],
    robust = !is.na(lineage) & !is.na(LFQ) & LFQ > 0 &
      !is.na(unique_peptides) & unique_peptides >= 1 &
      !is.na(identification_type) & identification_type == "By MS/MS")
})
write_csv(audit, file.path(OUTDIR, "protein_group_audit.csv"))
checks <- list()
for (rel in c("allocation_by_sample.csv", "proteomics/allocation_by_sample.csv", "plotting/allocation_by_sample.csv")) {
  f <- file.path(SULFURCAVE, rel)
  if (!file.exists(f)) next
  old <- read_csv(f, show_col_types = FALSE)
  for (metric in c(cats6, "n_detected", "n_robust")) {
    x <- summary_tbl %>% select(lineage, sample, reproduced = all_of(metric)) %>%
      left_join(old %>% select(lineage, sample, previous = all_of(metric)), by = c("lineage", "sample")) %>%
      mutate(reference = rel, metric = metric, difference = reproduced - previous,
             matches = !is.na(previous) & abs(difference) <= 1e-12 * pmax(1, abs(previous)))
    checks[[length(checks) + 1L]] <- x
  }
}
checks <- bind_rows(checks)
write_csv(checks, file.path(OUTDIR, "comparison_with_saved_results.csv"))
stopifnot(all(checks$matches), all(summary_tbl$total_growth_fraction <= 1),
          sum(summary_tbl$low_conf) == 1,
          summary_tbl$n_robust[summary_tbl$lineage == "SC_MAG_00006" & summary_tbl$sample == "Culture"] == 16,
          summary_tbl$n_robust[summary_tbl$lineage == "SC_MAG_00008" & summary_tbl$sample == "Culture"] == 50)
inputs <- c(PROT_XLSX, EGGNOG, KO_MAP, file.path(DF_DIR, c("Ferro", "MAG6", "Myco"), "defense_finder_systems.tsv"))
write_csv(tibble(path = inputs, md5 = unname(tools::md5sum(inputs))), file.path(OUTDIR, "input_checksums.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "R_sessionInfo.txt"))
message("Verified ", nrow(checks), " saved-result comparisons; all passed.")
print(summary_tbl %>% select(lineage, sample, n_robust, included, ribosomal_pct, growth_machinery_pct, CRISPR_Cas_pct))

# ---- 7. Figure --------------------------------------------------------------
pal   <- c(SC_MAG_00006 = "#fc8d62", SC_MAG_00008 = "#ffd92f", Mycobacterium = "#66c2a5")
# genus italicised (convention); strain/MAG suffix upright. Order matches lin_levels.
leg   <- list(
  expression(italic("Ferroplasma") ~ "MAG 6"),
  expression(italic("Ferroplasma") ~ "c. MAG 8"),
  expression(italic("Mycobacterium"))
)
mk_df <- function(metrics, labels) {
  bar_levels <- as.vector(t(outer(lin_levels, names(samples), paste, sep = "|")))
  summary_tbl %>%
    select(lineage, sample, low_conf, all_of(metrics)) %>%
    filter(!low_conf) %>%                     # drop low-confidence bars (kept in CSV)
    pivot_longer(all_of(metrics), names_to = "metric", values_to = "frac") %>%
    mutate(pct = 100 * frac,
           lineage = factor(lineage, levels = lin_levels),
           sample  = factor(sample, levels = names(samples)),
           bar     = factor(paste(lineage, sample, sep = "|"), levels = bar_levels),
           metric  = factor(metric, levels = metrics, labels = labels),
           fill_col = "white") %>%   # all bars white-filled; pattern carries the sample
    arrange(metric, bar)
}
# Single faceted layout: panel a = six growth-machinery categories, panel b =
# CRISPR-Cas. Faceting with ggh4x space = "free_x" makes each panel's width
# proportional to its number of x categories (a = 6, b = 1), so EVERY bar has the
# SAME width across both panels; independent = "y" gives each panel its own
# y-scale (a ~0-6 %, b ~0-25 %). One plot -> one shared legend, and identical
# pattern parameters render identically because bar width and plot height match.
# The two are still visually separate panels: the six growth categories are a
# mutually-exclusive partition, whereas CRISPR-Cas is an independent overlay
# (DefenseFinder Cas loci) on a much larger scale.
library(ggh4x)
cat_labs6 <- c("Ribosomal\nproteins","Translation\nfactors","aaRS",
               "Ribosome\nbiogenesis","DNA\nreplication","Cell\ndivision")
dfA <- mk_df(cats6, cat_labs6);              dfA$panel <- "a"   # growth machinery
dfB <- mk_df("CRISPR-Cas", "CRISPR-\nCas");  dfB$panel <- "b"   # CRISPR-Cas
dfC <- bind_rows(dfA, dfB)
dfC$panel  <- factor(dfC$panel,  levels = c("a", "b"))
dfC$metric <- factor(dfC$metric, levels = c(cat_labs6, "CRISPR-\nCas"))

# Sample -> pattern (ggpattern): Biofilm 1 = horizontal stripes, Biofilm 2 =
# 45 deg diagonal stripes, Culture (lab CH4) = 45 deg crosshatch. All bars are
# white-filled and every pattern is drawn in the lineage colour, so lineage =
# colour/outline and sample = pattern throughout.
# Spacing and density are mapped per sample (not fixed) because the three
# patterns need different scales to stay distinguishable at the narrow dodged
# bar width; the crosshatch in particular collapses to a near-solid bar at
# coarser spacing.
pat_map   <- c(Biofilm_1 = "stripe", Biofilm_2 = "stripe", Culture = "crosshatch")
pat_angle <- c(Biofilm_1 = 0,        Biofilm_2 = 45,       Culture = 45)
pat_space <- c(Biofilm_1 = 0.020,    Biofilm_2 = 0.045,    Culture = 0.015)
pat_dens  <- c(Biofilm_1 = 0.50,     Biofilm_2 = 0.50,     Culture = 0.52)

# ---- panel-b pattern-scale correction ---------------------------------------
# ggpattern's pattern_spacing is a fraction of the SHORTER side of each panel's
# viewport ("snpc") -- it is not an absolute length, and it does not track bar
# width. With space = "free_x" panel a is ~6x wider than panel b while both
# share the same height, so the shorter side is the HEIGHT in panel a but the
# WIDTH in panel b. The same nominal spacing therefore renders ~PANEL_K times
# finer in panel b (uncorrected, Biofilm 1's fine stripes collapsed to a
# near-solid bar there). Multiplying panel b's spacing by PANEL_K restores a
# matching absolute scale; pattern_density is a dimensionless duty cycle and
# needs no correction.
#   PANEL_K = min(w_a, h) / min(w_b, h), measured on the 14 x 6 in output at
#   300 dpi: panel a 2180 x 1127 px, panel b 421 x 1127 px. After correction
#   the measured stripe periods agree across panels to within ~4 %.
# NOTE: PANEL_K depends on the output size (14 x 6 in) and on the 6:1 split of
# x categories between the panels; recompute it if either changes.
PANEL_K <- 1127 / 421
dfC$pat_space_adj <- pat_space[as.character(dfC$sample)] *
                     ifelse(dfC$panel == "b", PANEL_K, 1)
pat_labs  <- c("Biofilm 1", "Biofilm 2", "Lab CH4")   # culture = "Lab CH4", no subscript

fig <- ggplot(dfC, aes(metric, pct, group = bar)) +
  geom_col_pattern(
    aes(colour = lineage, fill = fill_col,
        pattern = sample, pattern_angle = sample,
        pattern_spacing = pat_space_adj, pattern_density = sample,
        pattern_fill = lineage, pattern_colour = lineage,
        linetype = low_conf),
    position        = position_dodge2(width = 0.9, padding = 0.15, preserve = "single"),
    linewidth       = 0.4,
    pattern_key_scale_factor = 0.5) +
  scale_colour_manual(values = pal, labels = leg, name = "Lineage") +
  scale_fill_identity(guide = "none") +
  scale_pattern_manual(values = pat_map, labels = pat_labs, name = "Sample") +
  # same name + labels as scale_pattern_manual so all four pattern scales merge
  # into ONE "Sample" legend whose keys show the true angle/spacing/density
  scale_pattern_angle_manual(values = pat_angle, labels = pat_labs, name = "Sample") +
  scale_pattern_spacing_identity() +   # per-row values (panel-b corrected)
  scale_pattern_density_manual(values = pat_dens, labels = pat_labs, name = "Sample") +
  scale_pattern_fill_manual(values = pal, guide = "none") +
  scale_pattern_colour_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dotted"), guide = "none") +
  ggh4x::facet_grid2(cols = vars(panel), scales = "free",
                     independent = "y", space = "free_x") +
  labs(x = NULL, y = "% of robust protein-group LFQ", title = NULL) +
  theme_minimal(base_size = FS, base_family = FONT) +
  # all text: Arial, size 20, black; panel labels (strips) plain
  theme(text         = element_text(family = FONT, size = FS, colour = "black"),
        axis.text.x  = element_text(angle = 45, hjust = 1, size = FS, colour = "black"),
        axis.text.y  = element_text(angle = 0,  size = FS, colour = "black"),
        axis.title   = element_text(size = FS, colour = "black"),
        legend.text  = element_text(size = FS, colour = "black"),
        legend.title = element_text(size = FS, colour = "black"),
        legend.key.spacing.y = unit(0.15, "cm"),
        legend.position  = "right",
        strip.text       = element_text(family = FONT, size = FS, colour = "black", hjust = 0),
        strip.background = element_blank(),
        panel.spacing    = unit(1.4, "lines"),
        panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.6),
        axis.line.x  = element_line(colour = "black", linewidth = 0.6),
        axis.line.y  = element_line(colour = "black", linewidth = 0.6)) +
  guides(
    colour  = guide_legend(override.aes = list(fill = unname(pal[lin_levels]),
                                               pattern = "none")),
    # spacing is an identity scale now, so the keys need it stated explicitly
    pattern = guide_legend(override.aes = list(colour = "grey40",
                                               fill = "white",
                                               pattern_spacing = unname(pat_space),
                                               pattern_fill = "grey40",
                                               pattern_colour = "grey40")))

dir.create(dirname(OUTFILE), showWarnings = FALSE, recursive = TRUE)
ggsave(OUTFILE, fig, width = 14, height = 6, device = svglite::svglite)
message("wrote ", OUTFILE)
