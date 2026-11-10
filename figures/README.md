# Figures

Rendered figures, organised by the numbering in the manuscript
(*Cave_Ferro_Virus*, 84-page draft of Nov 2026). Each figure has its own folder:

- the folder itself holds the master file(s) for that figure;
- `panels/` holds individual panels or sub-assemblies for figures assembled by hand;
- `drafts/` holds previews, earlier versions and variants that are not in the manuscript;
- `source_data/` holds data tables found next to the figure files.

`MOVE_LOG.csv` records where every file came from (old path → new path) and why it was
placed where it is. Several files were renamed: supplementary figure files had been
numbered by an earlier draft (old S9–S12 and S15 are manuscript S10–S13 and S6), and
the Fig. 3 panel letters followed an older layout (see below).

## Index

| Figure | Folder | Script | Status |
|---|---|---|---|
| Fig. 1 | — | `plotting/Fig1_taxonomy_MAG_coverage.R` | no rendered file in the repository |
| Fig. 2 | `main/Fig2/` | `plotting/Fig2_pangenome_metabolism.R` | only a draft (size-adjusted 2b); 2c writes to `main/Fig2/panels/` |
| Fig. 3 | `main/Fig3/panels/` | `plotting/Fig3_metaproteomics_overview.R` (does not save the panels) | panels a–f + legends |
| Fig. 4 | `main/Fig4/` | `plotting/Fig4_proteome_pangenome_defense.R` | only an early draft panel |
| Fig. 5 | — | `plotting/Fig5_virus_host_CRISPR.R` (c: Fig. 4 script) | no rendered file |
| Fig. 6 | `main/Fig6/panels/` | `plotting/Fig6_tail_fiber_features.R` | panels |
| Fig. S1 | — | `plotting/FigS01_MAG_coverage_LabCH4.R` | no rendered file |
| Fig. S2 | — | none | no rendered file |
| Fig. S3 | — | `plotting/FigS03_COG_upset.R` | no rendered file |
| Fig. S4 | `supplementary/FigS04/` | `plotting/FigS04_rTCA_heatmap.R` | master (ggplot draft in `drafts/`) |
| Fig. S5 | `supplementary/FigS05/` | `plotting/FigS05_methane_heatmap.R` | master |
| S4 + S5 page | `supplementary/FigS04_S05_combined/` | none (script removed; not a manuscript figure) | archive |
| Fig. S6 | `supplementary/FigS06/` | `plotting/FigS06_proteome_allocation.R` | master (older version in `drafts/`) |
| Fig. S7 | — | `plotting/FigS07_MAG6_gene_map.R` | no rendered file |
| Fig. S8 | — | `plotting/FigS08_virus_host_abundance_noMAGs.R` | no rendered file |
| Fig. S9 | — | none | no rendered file |
| Fig. S10 | `supplementary/FigS10/` | `plotting/FigS10_physicochemical_features.R` | master |
| Fig. S11 | `supplementary/FigS11/` | `plotting/FigS11_aa_composition_radar.R` | master |
| Fig. S12 | `supplementary/FigS12/` | `plotting/FigS12_physicochemical_effect_sizes.R` | master |
| Fig. S13 | `supplementary/FigS13/` | `plotting/FigS13_aa_effect_sizes.R` | master |
| Fig. S14 | `supplementary/FigS14/` | `plotting/FigS14_myco_TM_pI.R` | master |
| Fig. S15 | — | none | no rendered file |

## Fig. 3 panel letters

The files in the former `fig3_out/` used an older layout. They were renamed to the
manuscript letters:

| Old file | Manuscript panel |
|---|---|
| `fig3a` | 3a — 3D LFQ scatter |
| `fig3b` | 3c — protein groups per MAG |
| `fig3c` | 3e — LFQ vs DNA coverage |
| `fig3d` | 3b — pathways, *Ca.* M. methanotrophicum |
| `fig3e` | 3d — pathways, *Ferroplasma* c. MAG 8 |
| `fig3f` | 3f — pathways, *Ferroplasma* MAG 6 |

`*_matched` files are variants of panels b, d, e and f kept alongside the plain versions; which
of the two is in the manuscript has not been checked.

`unassigned/` holds plots that match no manuscript figure.
