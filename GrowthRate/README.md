# R reproduction of proteome allocation and Figure S6

`../plotting/FigS06_proteome_allocation.R` replaces the Python workflow and performs both
the calculation and figure generation in R. It uses the existing Figure S6 (numbered S15 in an earlier draft)
ggplot2/ggpattern/ggh4x code, retaining colours, patterns, Arial typography,
two panels and separate y scales. The y-axis wording is clarified to
`% of robust protein-group LFQ`. SVG is the publication master; PNG is a preview.

## Run

From this GrowthRate folder, with the packages available in your R library:

```sh
Rscript ../plotting/FigS06_proteome_allocation.R
```

Arguments optionally specify the project root and output directory:

```sh
Rscript ../plotting/FigS06_proteome_allocation.R /path/to/project /path/to/output
```

Packages: readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2, ggpattern,
ggh4x, svglite. Install missing packages once using `install.packages()`.
The script itself does not install packages or make network calls.

To use packages installed in a separate R library, set `SULFURCAVE_R_LIB`:

```sh
SULFURCAVE_R_LIB=/path/to/R-library Rscript ../plotting/FigS06_proteome_allocation.R
```

## Inputs and definitions

Reads `proteomics/proteinGroups.xlsx` (original `ProteinGroups` sheet), eggNOG
annotations, `proteomics/ko_category_map.csv`, and the three DefenseFinder
system tables (Ferro, MAG6, Myco). Inputs are resolved from the project root.
Missing DefenseFinder tables cause an error, rather than an apparent zero.
The workbook type-guessing window covers all rows so sparse QC flags remain
visible; this fixes import warnings in the original R script.

After removing reverse, contaminant and only-by-site groups, a group is retained
per sample if LFQ > 0, identification type is `By MS/MS`, and unique peptides
>= 1. The denominator includes all retained groups for that lineage/sample,
including groups without a growth annotation. There is no shared-detection
requirement across the two biofilms. At least 50 robust groups are required
for inclusion in the figure; MAG 6 culture (16) is excluded and MAG 8 culture
(50) is retained. Excluded values remain in the output tables.

Growth categories use the project's fixed KO map and priority:
ribosomal > aaRS > translation factors > ribosome biogenesis > DNA replication
> cell division. Each group has at most one growth category; total growth
allocation includes ribosomal proteins. The map is BRITE-derived with a curated
cell-division category, rather than six exclusively BRITE-defined categories.

Majority protein IDs determine KO annotations and lineage, following the newer
`SupplementaryPlotting/FigS15_proteome_allocation_to_growth.R`. Groups spanning
more than one focal lineage are excluded; other non-focal source memberships
are handled as in that existing script. This is a reproduction, not a revised
protein-inference analysis.

CRISPR-Cas allocation independently uses DefenseFinder subtypes prefixed CAS,
matching locus tags to Majority protein IDs. It uses the same robust denominator
and can overlap growth categories. It is not the older COG-V defence proxy.

## Outputs and checks

The default output folder for the tables is `GrowthRate/reproduced_allocation_R/`:

- The figure itself is written to `figures/supplementary/FigS06/FigS06_proteome_allocation.svg`
  (or to the output folder, if one is given as the second argument).
- `allocation_summary.csv`: category fractions (0..1), counts, inclusion flags,
  and ribosomal, summed growth and Cas percentages (0..100).
- `allocation_by_category.csv`: per-category fractions and identification counts.
- `protein_group_audit.csv`: source Excel row, group IDs, category, lineage,
  sample LFQ, peptide evidence and robust flag.
- `comparison_with_saved_results.csv`: 216 category/count comparisons, all passed.
- `input_checksums.csv`, `R_sessionInfo.txt`: input hashes and R package versions.

The R figures were visually inspected. All six growth fractions and two counts
match all three existing result exports to 1e-12 tolerance. The Cas percentages
also match the preceding independent Python reproduction. Old Python-generated
tables remain in `reproduced_allocation` as a prior audit; new R outputs are in
the separately named folder above.

Environmental Ferroplasma ranges are **2.2--5.1% ribosomal** and **6.4--16.5%
total growth-associated allocation**. MAG 6 has 15.3--16.5% growth allocation
versus MAG 8's 6.4--8.0%. Cas allocation is 1.6--2.4% in MAG 6 and 13.2--25.2%
in MAG 8. These are relative LFQ allocations, not measured growth rates.
