setwd("/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/Sulfurcave_ferro_virus_git/")
allocation_df <- read.delim("proteomics/allocation_by_sample.csv", sep = ',', check.names = F)

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpattern)

cats <- c("Aminoacyl-tRNA synthetases", "Cell division", "DNA replication",
          "Ribosomal proteins", "Ribosome biogenesis", "Translation factors")

plot_df <- allocation_df %>%
  pivot_longer(cols = all_of(cats), names_to = "category", values_to = "value")

ggplot(plot_df, aes(x = category, y = value, fill = lineage)) +
  geom_col(position = position_dodge(preserve = "single"), colour = "black", linewidth = 0.3) +
  facet_wrap(~ sample, ncol = 1) +                 # one panel per sample
  labs(x = NULL, y = "Value", fill = "Lineage") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
