# =============================================================================
# General methane-metabolism marker-protein heatmap  --  ComplexHeatmap version
#
# The methane counterpart of FigS04_rTCA_heatmap.R. Same design
# principle, applied to the methane figure from plotting/Figure03_5.R:
#   - proteins are kept if their KEGG ortholog (KO) belongs to a KEGG
#     methane-metabolism MODULE (map00680, fetched live via KEGGREST);
#   - columns are labelled  "<Module(s)> (KO)"  and split into three panels:
#       Methanotroph related | Methanogenesis related | Other
#     (from the module names, exactly as the original methane_general plot);
#   - rows are the genomes carrying >=1 methane-module KO (29 genomes; this set
#     includes Babela MAG 20, NA MAG 21, NA MAG 31 which the rTCA marker set
#     excluded).
#
# Styling matches the rTCA ComplexHeatmap supplementary figure:
#   * rows relabelled "<genus> MAG <n>" (genus italic where Table 2 italicises;
#     blank genus -> "NA"; the 5 in-house reference genomes get dual
#     "<original label> / <genus> MAG <n>" labels);
#   * three-state cells  colour(viridis)=detected / grey=present-not-detected /
#     white=gene absent;
#   * rows grouped by genus, pathway group titles under the x axis (0 deg),
#     all text Arial 20, no plot title;
#   * legends: "Mean protein intensity (LFQ; log 10)" and "Gene presence"
#     (Yes=grey / No=white).
#
# Depends on analysis/marker_protein_heatmap_data.R (sourced) for the shared data pipeline
# (df_expanded, tax2bin, mag_num, genus_by_mag, italic_genera).  Set `repo`.
# =============================================================================

# ---- Paths (see README, "Running the scripts") ----
REPO     <- Sys.getenv("SULFURCAVE_REPO", unset = here::here())  # repository root

suppressPackageStartupMessages({
  library(KEGGREST); library(dplyr); library(tidyr); library(stringr)
  library(ComplexHeatmap); library(circlize)
  library(grid); library(gridtext); library(viridisLite); library(svglite)
})

repo   <- REPO
outdir <- file.path(repo, "plotting")
FIGOUT_S05 <- file.path(repo, "figures/supplementary/FigS05")  # figure output folder
dir.create(FIGOUT_S05, recursive = TRUE, showWarnings = FALSE)
FONT <- "Arial"; FS <- 20

# ---- shared pipeline (df_expanded: proteins x KO with genus + LFQ) ----------
source(file.path(repo, "analysis", "marker_protein_heatmap_data.R"), echo = FALSE)

# ---- 1. methane-metabolism KO -> module map (KEGG map00680) -----------------
methane_pathway <- keggGet("map00680")[[1]]
modules_ids <- names(methane_pathway$MODULE)
ko_module_list <- list()
for (mod_id in modules_ids) {
  mod <- keggGet(mod_id)[[1]]
  if (!is.null(mod$DEFINITION)) {
    kos <- names(mod$ORTHOLOGY)
    if (length(kos) > 0)
      ko_module_list[[length(ko_module_list) + 1]] <-
        data.frame(KO = kos, Module = mod_id, Module_name = mod$NAME, stringsAsFactors = FALSE)
  }
}
df_ko_expanded <- do.call(rbind, ko_module_list) %>%
  mutate(KO = str_replace_all(KO, "\\+", ",")) %>%
  separate_rows(KO, sep = ",") %>%
  mutate(KO = str_trim(KO), KO = str_extract(KO, "^K[0-9]{5}")) %>%
  select(KO, Module, Module_name) %>%
  group_by(KO) %>%
  summarise(Modules      = paste(unique(Module), collapse = ";"),
            Module_names = paste(unique(Module_name), collapse = ";"), .groups = "drop")

# ---- 2. keep proteins whose KO is in a methane module -----------------------
df_j  <- df_expanded %>% left_join(df_ko_expanded, by = "KO")
df_jf <- df_j[-which(is.na(df_j$Modules)), ]

# ---- 3. one cell per (genome, Module(KO)) + methane category ----------------
df_heat_m <- df_jf %>%
  rowwise() %>%
  mutate(avg_LFQ = mean(c(log10_LFQ_intensity_Cave_Biofilm_1,
                          log10_LFQ_intensity_Cave_Biofilm_2), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(avg_LFQ  = ifelse(is.nan(avg_LFQ), NA_real_, avg_LFQ),
         KO_role  = paste0(Modules, " (", KO, ")"),
         category = case_when(
           grepl("methanotroph|methane oxidation", Module_names, ignore.case = TRUE) ~ "Methanotroph related",
           grepl("methanogenesis|methane production", Module_names, ignore.case = TRUE) ~ "Methanogenesis related",
           TRUE ~ "Other"))

df_cell <- df_heat_m %>%
  group_by(tax_id, KO_role, category) %>%
  summarise(avg_LFQ = ifelse(all(is.na(avg_LFQ)), NA_real_, mean(avg_LFQ, na.rm = TRUE)), .groups = "drop") %>%
  mutate(status = ifelse(is.na(avg_LFQ), "not_detected", "measured"))

# ---- 4. row labels "<genus> MAG <n>" (rTCA principle, extended) -------------
tax2bin_m <- c(tax2bin, SCBin00021 = "SC_Bin_00021", SCBin00031 = "SC_Bin_00031",
               SCMAG00020 = "SC_MAG_00020")
mag_num_m <- function(tid) as.integer(sub("^0+", "", sub(".*_", "", tax2bin_m[[tid]])))
# genera for methane-only MAGs (Table 2): MAG 20 = Babela; MAG 21, 31 blank
genus_by_mag[["20"]] <- "Babela"; genus_by_mag[["21"]] <- ""; genus_by_mag[["31"]] <- ""
italic_genera <- setdiff(unique(unlist(genus_by_mag)), "")   # include methane-only genera (Babela)

base_label_m <- function(tid) {
  n <- mag_num_m(tid); g <- genus_by_mag[[as.character(n)]]
  if (is.null(g) || g == "") return(paste0("NA MAG ", n))
  if (g %in% italic_genera)  return(paste0("*", g, "* MAG ", n))
  paste0(g, " MAG ", n)
}
taxa <- unique(df_cell$tax_id)
lab  <- setNames(vapply(taxa, base_label_m, character(1)), taxa)
lab["KDJLIKBO"] <- "*Ca*. M. Methanotrophicum/*Mycobacterium* MAG 16"
lab["EKPBOPOH"] <- "*Cuniculiplasma*-EKPBOPOH/*Cuniculiplasma* MAG 18"
lab["KNPMNEEE"] <- "*Ferroplasma* MAG 8"
lab["Crude2"]   <- "*Mycobacterium* bin 2/*Mycobacterium* MAG 19"
lab["Crude3"]   <- "*Mycobacterium* bin 3/*Mycobacterium* MAG 30"
df_cell$row_label <- lab[df_cell$tax_id]
df_cell$order_key <- vapply(df_cell$tax_id, function(t) mag_num_m(t) + ifelse(t == "Crude3", 0.5, 0), numeric(1))
rl_top <- df_cell %>% distinct(row_label, order_key) %>% arrange(order_key) %>% pull(row_label)

# ---- 5. matrices (rows top->bottom = MAG ascending) -------------------------
# column groups: three functional groups (assigned from each KO's KEGG module):
#   Methane oxidation (methanotrophy)          <- M00174
#   Methanogenesis (+ cofactors)               <- M00356/357/563/567 + M00358/378/935
#   Central carbon & C1 assimilation           <- M00020/344/345/346/608
grp_methanogen <- c("M00356","M00357","M00563","M00567","M00358","M00378","M00935")
grp3 <- function(kr){
  mods <- strsplit(sub(" \\(K[0-9]{5}\\)$", "", kr), ";")[[1]]
  if ("M00174" %in% mods) "Methane oxidation (methanotrophy)"
  else if (any(mods %in% grp_methanogen)) "Methanogenesis (+ cofactors)"
  else "Central carbon & C1 assimilation (redox / biosynthesis)"
}
block_levels <- c("Methane oxidation (methanotrophy)",
                  "Methanogenesis (+ cofactors)",
                  "Central carbon & C1 assimilation (redox / biosynthesis)")
col_tab <- df_cell %>% distinct(KO_role, category) %>%
  mutate(ko        = sub("^.*\\((K[0-9]+)\\)$", "\\1", KO_role),
         first_mod = sub("^(M[0-9]{5}).*", "\\1", sub(" \\(K[0-9]{5}\\)$", "", KO_role)),
         block     = vapply(KO_role, grp3, character(1)),
         block     = factor(block, levels = block_levels)) %>%
  arrange(block, first_mod, KO_role)
cl        <- col_tab$KO_role
col_block <- col_tab$block
col_ko    <- col_tab$ko
# KO IDs cited in the manuscript text -> drawn bold on the x-axis
mentioned_ko <- c("K01006","K01007","K00024","K00031","K01681","K01902","K01903","K00174","K00175",
                  "K03841","K01624","K00600","K00058","K00831","K16157","K16158","K16159","K16161",
                  "K08093","K08094","K21071","K01595","K01689")
ko_face <- ifelse(col_ko %in% mentioned_ko, "bold", "plain")

# panel titles = the functional-group names, wrapped for a compact bottom strip
block_titles <- vapply(levels(col_block), function(s)
  paste(strwrap(s, width = 22), collapse = "\n"), character(1))

mat  <- matrix(NA_real_, length(rl_top), length(cl), dimnames = list(rl_top, cl))
stat <- matrix("absent",  length(rl_top), length(cl), dimnames = list(rl_top, cl))
for (i in seq_len(nrow(df_cell))) {
  a <- df_cell$row_label[i]; b <- df_cell$KO_role[i]
  mat[a, b] <- df_cell$avg_LFQ[i]; stat[a, b] <- df_cell$status[i]
}

# ---- 6. row grouping by genus ----------------------------------------------
row_tax <- sapply(rl_top, function(l) df_cell$tax_id[match(l, df_cell$row_label)])
row_gen <- sapply(row_tax, function(t){ g <- genus_by_mag[[as.character(mag_num_m(t))]]
                                        if (is.null(g) || g == "") "NA" else g })
mn <- tapply(sapply(row_tax, mag_num_m), row_gen, min)
gen_levels <- names(sort(mn)); gen_levels <- c(setdiff(gen_levels, "NA"), if ("NA" %in% gen_levels) "NA")
row_split  <- factor(unname(row_gen), levels = gen_levels)
# genus group titles: italic; bold-italic if the genus is cited in the manuscript text
mentioned_genera <- c("Ferroplasma", "Acidithiobacillus", "Mycobacterium")
gen_title  <- vapply(gen_levels, function(g) {
  if (g == "NA") return("NA")
  if (g %in% mentioned_genera) paste0("***", g, "***") else paste0("*", g, "*")
}, character(1))

# ---- 7. colours, legends, cell painter -------------------------------------
lfq_max <- max(mat, na.rm = TRUE)
col_fun <- colorRamp2(seq(0, lfq_max, length.out = 5), viridis(5))
cell_fun <- function(j, i, x, y, w, h, fill) {
  s <- stat[i, j]
  fc <- if (s == "measured") col_fun(mat[i, j]) else if (s == "not_detected") "grey80" else "white"
  grid.rect(x, y, w, h, gp = gpar(fill = fc, col = "white", lwd = 0.6))
}
ttl_gp <- gpar(fontfamily = FONT, fontsize = FS, fontface = "plain", col = "black")
lab_gp <- gpar(fontfamily = FONT, fontsize = FS)
lgd_val  <- Legend(col_fun = col_fun, title = "Mean protein intensity (LFQ; log 10)",
                   at = round(seq(0, lfq_max, length.out = 5), 1),
                   title_gp = ttl_gp, labels_gp = lab_gp,
                   grid_width = unit(6, "mm"), legend_height = unit(45, "mm"))
lgd_grey <- Legend(labels = c("Yes", "No"), title = "Gene presence",
                   legend_gp = gpar(fill = c("grey80", "white")), border = "grey40",
                   title_gp = ttl_gp, labels_gp = lab_gp,
                   grid_height = unit(8, "mm"), grid_width = unit(8, "mm"),
                   row_gap = unit(6, "mm"))

row_md <- gt_render(rl_top, gp = gpar(fontfamily = FONT, fontsize = FS))
ht <- Heatmap(
  mat, name = "LFQ", col = col_fun, na_col = "white",
  rect_gp = gpar(type = "none"), cell_fun = cell_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_order = rl_top, column_order = cl,
  row_split = row_split, cluster_row_slices = FALSE,
  row_title = gt_render(gen_title, gp = gpar(fontfamily = FONT, fontsize = FS)),
  row_title_rot = 0, row_title_side = "right", row_gap = unit(1.2, "mm"),
  column_split = col_block, cluster_column_slices = FALSE,
  column_title = block_titles, column_title_side = "bottom",
  column_title_gp = gpar(fontfamily = FONT, fontsize = FS), column_title_rot = 0,
  column_labels = col_ko,
  column_names_gp = gpar(fontfamily = FONT, fontsize = FS, fontface = ko_face), column_names_rot = 45,
  row_labels = row_md,
  row_names_gp = gpar(fontfamily = FONT, fontsize = FS), row_names_side = "left",
  show_heatmap_legend = FALSE, border = TRUE, column_gap = unit(1.5, "mm"))

draw_it <- function() draw(ht, annotation_legend_list = list(lgd_val, lgd_grey),
                           annotation_legend_side = "right", legend_gap = unit(10, "mm"),
                           padding = unit(c(2, 165, 2, 2), "mm"))  # bottom, left, top, right

W <- 36; H <- 13
cairo_pdf(file.path(FIGOUT_S05, "FigS05_methane_heatmap.pdf"), width = W, height = H); draw_it(); dev.off()
png(file.path(FIGOUT_S05, "FigS05_methane_heatmap.png"), width = W, height = H, units = "in", res = 300, type = "cairo"); draw_it(); dev.off()
svglite(file.path(FIGOUT_S05, "FigS05_methane_heatmap.svg"), width = W, height = H); draw_it(); dev.off()
