library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)

# --- 1. DATA LOADING & PREPARATION ---
no_mags_file <- '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/vcat_iphop_abundance_202602/vcat_iphop_abundance_no_mags.txt'
with_mags_file <- '/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/vcat_iphop_abundance_202602/vcat_iphop_abundance_with_mags.txt'
  
df_no_mags <- read.csv2(no_mags_file, sep = '\t')
df_with_mags <- read.csv2(with_mags_file, sep = '\t')

df <- df_no_mags

df <- df %>%
  mutate(across(starts_with("Sample"), ~ as.numeric(as.character(.x))))

# Convert from wide to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with("Sample"), 
    names_to = "sample", 
    values_to = "abundance"
  ) %>%
  mutate(sample = gsub("_abundance", "", sample)) 

# --- 2. TRANSFORMATION TO LODES FORMAT ---
lodes_all <- df_long %>%
  transmute(contig_id, Virus = vcat_class, Host = iphop_genus, weight = abundance, sample) %>%
  ggalluvial::to_lodes_form(
    axes  = c("Virus", "Host"),
    key   = "axis",
    value = "stratum",
    id    = "contig_id"
  ) %>%
  mutate(
    stratum = as.character(stratum),
    # Handle missing values as "Unclassified"
    stratum = ifelse(is.na(stratum) | stratum == "", "Unclassified", stratum),
    stratum = factor(stratum)
  )

# --- 3. AESTHETICS (MAPS & COLORS) ---
name_map <- c(
  "Sample_ERR10036468" = "Biofilm 1",
  "Sample_ERR10036469" = "Biofilm 2",
  "Sample_ERR10036470" = "Lab CH4"
)

my_colors <- c(
  # --- Virus Taxa ---
  "Caudoviricetes" = "#e7298a",
  "Faserviricetes" = "#2ca02c",
  "Megaviricetes"  = "#d62728",
  "Naldaviricetes" = "#1f77b4",
  "Below 1%"       = "#8da0cb",   # virus block
  "Unclassified"   = "gray",      # virus block
  
  # --- Host Genera ---
  "Acidithiobacillus"     = "#fc8d62",   # LOCKED (group 2)
  "Below 1%"              = "#8da0cb",    # LOCKED (group 2)
  "Cuniculiplasma"        = "#e5c494",    # LOCKED (group 2)
  "Ferroplasma"           = "#ffd92f",    # LOCKED (group 2)
  "Mycobacterium"         = "#66c2a5",    # LOCKED (group 2)
  "Unclassified"          = "#b3b3b3",    # LOCKED (group 2)
  
  "Escherichia"           = "#C49C94FF",  # kept (was fine)
  "JAKAFX01"              = "#9EDAE5FF",   # kept (was fine)
  "Sulfobacillus"         = "#98DF8AFF",   # kept (was fine)
  
  "Igneacidithiobacillus" = "#4d004b",    # CHANGED (was #decbe4, clashed w/ Marinobacter)
  "Marinobacter"          = "#005a32",    # CHANGED (was #fddaec, clashed w/ Igneacidithiobacillus)
  "Lachnospira"           = "#525252",    # CHANGED (was #e5d8bd, clashed w/ Cuniculiplasma/Streptomyces/Mobilitalea)
  "Mobilitalea"           = "#02818a",    # CHANGED (was #DBDB8DFF, clashed w/ Streptomyces/Lachnospira)
  "Streptomyces"          = "#993404"     # CHANGED (was #ffffcc, clashed w/ Mobilitalea/Lachnospira)
)

# --- 4. PLOTTING ---
p <- ggplot(lodes_all,
            aes(x = axis, stratum = stratum, alluvium = contig_id, 
                y = weight * 100, fill = stratum)) +
  # Flows
  geom_alluvium(alpha = 0.6, knot.pos = 0.4) +
  # Vertical bars
  geom_stratum(width = 0.3, color = "grey30") +
  # Scales & Colors
  scale_x_discrete(limits = c("Virus", "Host"), expand = c(.12, .12)) +
  scale_fill_manual(values = my_colors, name = "Virus / Host") +
  # Multi-panel display by Sample
  facet_wrap(~ sample, nrow = 1, labeller = labeller(sample = name_map)) +
  # Labels
  labs(
    x = "Virus class --> Host genus", y = "% mapped reads",
    # title = "Virus and predicted host abundance",
    fill = "Taxa / Hosts"
  ) +
  # Theme Customization
  theme_minimal(base_family = "Arial", base_size = 20) +
  theme(
    legend.position = "none",
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    panel.grid.major.x = element_blank() # Clean up vertical lines
  )

p

# ggsave(
#   filename = "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/figures/virus_taxa_host_abundance_202602/no_mags_color_legend.svg", 
#   plot = p, 
#   device = "svg", 
#   width = 16, 
#   height = 10, 
#   units = "in",
#   limitsize = FALSE
# )
# legend
library(ggplot2)
library(cowplot)

virus_colors <- c(
  "Caudoviricetes" = "#e7298a",
  "Faserviricetes" = "#2ca02c",
  "Megaviricetes"  = "#d62728",
  "Naldaviricetes" = "#1f77b4",
  "Below 1%"       = "#8da0cb",
  "Unclassified"   = "#b3b3b3"
)
host_colors <- c(
  "Acidithiobacillus"     = "#fc8d62",
  "Below 1%"              = "#8da0cb",
  "Cuniculiplasma"        = "#e5c494",
  "Escherichia"           = "#C49C94FF",
  "Ferroplasma"           = "#ffd92f",
  "Igneacidithiobacillus" = "#4d004b",
  "JAKAFX01"              = "#9EDAE5FF",
  "Lachnospira"           = "#525252",
  "Marinobacter"          = "#005a32",
  "Mobilitalea"           = "#02818a",
  "Mycobacterium"         = "#66c2a5",
  "Streptomyces"          = "#993404",
  "Sulfobacillus"         = "#98DF8AFF",
  "Unclassified"          = "#b3b3b3"
)

# entries that are NOT taxon names and should stay upright
non_taxa <- c("Below 1%", "Unclassified")

# build italic expression labels for taxon names, plain for the rest
make_labels <- function(cols, italicize) {
  nm <- names(cols)
  if (!italicize) return(nm)
  setNames(
    lapply(nm, function(x) if (x %in% non_taxa) x else bquote(italic(.(x)))),
    nm
  )
}

legend_only <- function(cols, title, italicize = FALSE) {
  d <- data.frame(cat = factor(names(cols), levels = names(cols)), x = 1, y = 1)
  ggplot(d, aes(x, y, fill = cat)) +
    geom_col() +
    scale_fill_manual(values = cols, name = title,
                      labels = make_labels(cols, italicize),
                      guide = guide_legend(ncol = 1, byrow = TRUE)) +
    theme_void(base_family = "Arial", base_size = 20) +
    theme(
      legend.justification = c(0, 1),
      legend.title = element_text(face = "bold", hjust = 0),
      legend.text  = element_text(),
      legend.key.size = unit(0.65, "cm"),
      legend.key.spacing.y = unit(0.35, "cm")
    )
}

leg_v <- get_legend(legend_only(virus_colors, "Virus class"))                 # virus stays upright
leg_h <- get_legend(legend_only(host_colors,  "Host genus", italicize = TRUE)) # host genera italic
combined <- plot_grid(leg_v, leg_h, ncol = 2, rel_widths = c(0.9, 1))
combined
# ggsave("legends_virus_host.png", combined, width = 7.5, height = 5, dpi = 300, bg = "white")



