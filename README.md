# Sulfurcave_ferro_virus

Code and processed data for the multi-omic characterisation of the ultra-acidic
(pH ~1), gas-driven biofilm of the Sulfur Cave (Puturosu, Romania).

> Wu L.-Y., Khot V.M., Gallant J., Flot J.-F., Gherman V.D., Bodelier P.,
> Sarbu S.M., Bitter W., Dutilh B.E., Melkonian C. *<manuscript title>*.
> In submission (2026). <!-- TODO: final title, journal, DOI -->

The analyses focus on the metabolic interdependence between the
sulfur-oxidising archaeon *Ferroplasma* and the methanotroph
*Ca.* Mycobacterium methanotrophicum, and on the physicochemical adaptations
of the cave virome.

- **Genome-resolved metagenomics**: MAG recovery, taxonomy and pangenomics (anvi'o).
- **Metaproteomics**: functional profiling and proteome allocation of the dominant lineages.
- **Viromics and host linkage**: viral contigs, CRISPR spacer–protospacer matching.
- **Structural proteomics of phage proteins**: pI, GRAVY, aromaticity and amino-acid
  composition of cave versus IMG/VR phage proteins.

## Repository layout

| Folder | Contents |
|---|---|
| `metagenomics/` | MAG summaries, coverage, taxonomy (CAT, Kaiju), eggNOG annotations, DefenseFinder output |
| `pangenomics/` | anvi'o pangenome summaries, KO/OG matrices, enrichment tables, GapMind inputs |
| `proteomics/` | MaxQuant `proteinGroups.xlsx`, CD-HIT cluster files, KO category map, derived tables |
| `viromics/` | viral contig features, CRISPRTarget report, cave vs IMG/VR comparisons |
| `GrowthRate/` | proteome-allocation tables and checks written by the Fig. S6 script |
| `plotting/` | one script per manuscript figure (`FigN_*.R`, `FigSNN_*.R`) |
| `analysis/` | analysis and data-preparation scripts |
| `figures/` | rendered figures by manuscript number (`main/Fig1–6`, `supplementary/FigS01–S15`); see `figures/README.md` |

## Script → figure map

Figure numbers follow the manuscript; each figure has at most one script.

| Figure | Script | Output folder |
|---|---|---|
| Fig. 1 | `plotting/Fig1_taxonomy_MAG_coverage.R` | — |
| Fig. 2 (a–c) | `plotting/Fig2_pangenome_metabolism.R` | `figures/main/Fig2/` |
| Fig. 3 | `plotting/Fig3_metaproteomics_overview.R` | `figures/main/Fig3/` |
| Fig. 4 | `plotting/Fig4_proteome_pangenome_defense.R` | `figures/main/Fig4/` |
| Fig. 5 (a, b, d) | `plotting/Fig5_virus_host_CRISPR.R` (panel c is built in the Fig. 4 script) | — |
| Fig. 6 | `plotting/Fig6_tail_fiber_features.R` | `figures/main/Fig6/panels/` |
| Fig. S1 | `plotting/FigS01_MAG_coverage_LabCH4.R` | — |
| Fig. S3 | `plotting/FigS03_COG_upset.R` | — |
| Fig. S4 | `plotting/FigS04_rTCA_heatmap.R` | `figures/supplementary/FigS04/` |
| Fig. S5 | `plotting/FigS05_methane_heatmap.R` | `figures/supplementary/FigS05/` |
| Fig. S6 | `plotting/FigS06_proteome_allocation.R` (tables to `GrowthRate/reproduced_allocation_R/`) | `figures/supplementary/FigS06/` |
| Fig. S7 | `plotting/FigS07_MAG6_gene_map.R` | — |
| Fig. S8 | `plotting/FigS08_virus_host_abundance_noMAGs.R` | — |
| Fig. S10 | `plotting/FigS10_physicochemical_features.R` | `figures/supplementary/FigS10/` |
| Fig. S11 | `plotting/FigS11_aa_composition_radar.R` | `figures/supplementary/FigS11/` |
| Fig. S12 | `plotting/FigS12_physicochemical_effect_sizes.R` | `figures/supplementary/FigS12/` |
| Fig. S13 | `plotting/FigS13_aa_effect_sizes.R` | `figures/supplementary/FigS13/` |
| Fig. S14 | `plotting/FigS14_myco_TM_pI.R` | `figures/supplementary/FigS14/` |

No script exists in the repository for Fig. 2d (tree + ANI), Fig. S2 (GapMind),
Fig. S9 (gene-sharing network) or Fig. S15 (structure prediction). <!-- TODO -->

Analysis scripts (no manuscript figure of their own) are in `analysis/`:

| Script | Purpose |
|---|---|
| `analysis/Fig2_pangenome_statistics.R` | pangenome statistics and KEGG-module table (`pangenomics/result_SD_pw.csv`) behind Fig. 2 and the Results text; also contains the original draft Fig. 2 plots |
| `analysis/marker_protein_heatmap_data.R` | marker-protein table shared by Figs S4 and S5 (sourced by both) |
| `analysis/compute_physicochemical_features.py` | per-protein features for Figs 6 and S10–S13 (`viromics/cave_imgvr_features_raw.csv`) |
| `analysis/gRodon.R` | codon-usage growth-rate prediction |

## Running the scripts

All paths are built from two variables defined at the top of each script:

- `REPO`: the repository root. It is found automatically with the
  [`here`](https://here.r-lib.org) package (via the `.here` file), so you can
  run scripts from any folder inside the repository. Override it with
  `SULFURCAVE_REPO=/path/to/repo`.
- `EXTERNAL`: inputs that are too large or not ours to redistribute (see
  below). Defaults to `<repo>/external_data/` (git-ignored). Override it with
  `SULFURCAVE_EXTERNAL=/path/to/data`.

```sh
git clone https://github.com/lingyi-owl/Sulfurcave_ferro_virus.git
cd Sulfurcave_ferro_virus
Rscript plotting/FigS12_physicochemical_effect_sizes.R
```

### Inputs expected under `EXTERNAL`

<!-- TODO: deposit these (e.g. Zenodo) or move them into the repository, then list the DOI here -->

| Path under `EXTERNAL` | Used by |
|---|---|
| `mycobacterium_pangenome/` (`MYCO_Pan_gene_clusters_summary.txt`, `merged_TM_protein_properties.tsv`) | `FigS14_myco_TM_pI.R` |
| `2026_analysis/SulfurCave_imgvr_proteins/all/` | `analysis/compute_physicochemical_features.py` |
| `vcat_iphop_abundance_202602/` | `Fig5_virus_host_CRISPR.R`, `FigS08_virus_host_abundance_noMAGs.R` |
| `ferroplasma_CRISPR_annotations/`, `ferroplasma_CRISPR/`, `sulfur_cave_data/` | `Fig5_virus_host_CRISPR.R`, `analysis/gRodon.R` |
| `CRISPR/`, `GENOMES/Ferroplasma_complete_genome/` | `FigS07_MAG6_gene_map.R` |
| `Defencefinder/` | `Fig4_proteome_pangenome_defense.R` (DefenseFinder copies) |
| `PANGENOMES/` | second half of `analysis/Fig2_pangenome_statistics.R` |
| `figures/Fig4/` (intermediate tables) | `Fig4_proteome_pangenome_defense.R` |

### Software

R ≥ 4.3. Main CRAN packages: ggplot2, dplyr, tidyr, readr, readxl, stringr,
purrr, data.table, patchwork, cowplot, ggbeeswarm, ggrepel, ggsignif, ggh4x,
ggpattern, ggtext, ggalluvial, gggenes, geomtextpath, circlize, UpSetR,
ComplexUpset, svglite, here. Bioconductor: ComplexHeatmap, Biostrings,
KEGGREST. GitHub: ggradar (ricardo-bion/ggradar), gRodon2
(jlw-ecoevo/gRodon2). Python ≥ 3.9 with numpy, pandas and Biopython.
<!-- TODO: add sessionInfo() / renv.lock and versions of external tools
(anvi'o, eggNOG-mapper, Kaiju, CAT, DefenseFinder, CRISPRTarget, GapMind, MaxQuant) -->

## Data availability

Raw reads are not stored here. Biofilm metagenomes: ENA BioProject
PRJEB45004 (runs ERR10036468–ERR10036470). Enrichment-culture reads and
assembly: NCBI BioProject PRJNA837300. Proteomics: PRIDE <!-- TODO: PXD accession -->.

## Licence and citation

Code is released under the MIT licence (`LICENSE`). If you use this
repository, please cite the article above (see `CITATION.cff`).
