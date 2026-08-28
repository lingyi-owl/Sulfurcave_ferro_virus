# Methods — Proteome allocation to translation and growth machinery

## Metaproteomic quantification and MAG assignment

Protein groups were taken from the MaxQuant `proteinGroups` output. Reverse-database
hits, potential contaminants, and proteins identified only by a modified site were
removed (leaving 3,276 protein groups, identical to the workbook's
`ProteinGroups filtered` sheet). Because iBAQ intensities were not available in this
search, label-free quantification (LFQ) intensities were used throughout.

Allocation was computed per sample for three MS runs, which were kept separate rather
than summed: two cave-biofilm whole-cell-lysate fractions (`cave_wcl_1`, `cave_wcl_2`;
"Biofilm 1" and "Biofilm 2") and one *in vitro* culture (`invitro`; "Culture"). The two
biofilm supernatant fractions (`cave_sup_1/2`) were not used.

Protein groups were assigned to source populations using the MaxQuant "Majority protein
IDs", whose locus tags carry a genome/MAG prefix. Three lineage-specific populations were
quantified: *Ferroplasma* SC_MAG_00006 (contig prefix SCMAG00006); *Ferroplasma*
SC_MAG_00008 (the complete circular genome; locus prefix KNPMNEEE); and *Mycobacterium*
(the complete *Mycobacterium Methanotrophicum* genome SC_MAG_00016; locus prefix
KDJLIKBO). Prefix-to-taxon assignments follow `info_table_cleanids.csv`. A protein group
was assigned to a population if any of its majority protein IDs mapped to that
population's genome; groups shared between two of the target populations, or shared with a
non-target taxon, were excluded from the per-population sums. Population separation was
verified before allocation: only 34 of 1,436 *Ferroplasma* groups (2.4%) were shared
between SC_MAG_00006 and SC_MAG_00008, and no protein group was shared between
*Ferroplasma* and *Mycobacterium*.

## Identification confidence

For each population in each sample, identifications were classified as *robust* when the
protein group was sequenced in that sample (MaxQuant identification type "By MS/MS") and
carried at least one unique peptide there, as opposed to being inferred only by
match-between-runs ("By matching"). A population×sample combination resting on fewer than
30 robust identifications was flagged as low-confidence (dotted outline in the figure).
This flag is material for the culture: *Ferroplasma* SC_MAG_00006 is only sparsely and
largely non-robustly detected in the *in vitro* culture (59 groups, 16 robust) and its
culture allocation is therefore unreliable, whereas *Ferroplasma* SC_MAG_00008 (56 groups,
50 robust) and *Mycobacterium* (468 groups, 465 robust) are well supported there.

## Functional classification

Open reading frames were annotated with eggNOG-mapper, and each protein group inherited
the KEGG Orthology (KO) identifiers of its constituent ORFs. Six growth-machinery
categories were defined by KO membership in KEGG BRITE hierarchies retrieved from the KEGG
REST API: ribosomal proteins (the "Ribosomal proteins" section of ko03011, all taxa,
including the archaeal set), translation factors (ko03012; initiation, elongation and
release factors), aminoacyl-tRNA synthetases (the aaRS section of ko01007), ribosome
biogenesis factors (ko03009), and DNA replication proteins (ko03032). Cell division was
defined from a curated set of 18 KEGG-verified division and septation orthologs (FtsZ,
FtsA, FtsQ, FtsI, FtsK, FtsN, FtsB, FtsW, MinC, MinD, ZapA, the ParAB partitioning
proteins, and the Mur/Ddl peptidoglycan ligases); *Ferroplasma* (Thermoplasmatales)
divides via an FtsZ-based mechanism, so the archaeal Cdv system and the shape determinant
MreB were deliberately excluded. Categories were made mutually exclusive by assigning each
KO to a single category following the priority order ribosomal proteins > aminoacyl-tRNA
synthetases > translation factors > ribosome biogenesis > DNA replication > cell division,
so that a protein contributing to more than one category was counted once. The resulting
KO-to-category map (1,176 KOs) and its full provenance are provided as
`ko_category_map.csv` and `ko_map_provenance.txt`.

## Allocation metrics

Allocation was computed using only robustly-sequenced protein groups — those identified
"By MS/MS" with at least one unique peptide in the given sample (see Identification
confidence) — for both the category numerator and the lineage-total denominator;
match-between-runs-only quantifications were excluded. For each population in each sample,
proteome allocation to a category was computed as the summed LFQ intensity of the robust
protein groups assigned to that category, divided by the total summed LFQ intensity of all
robust protein groups assigned to that population in that sample:

    allocation(category) = Σ intensity(category proteins in MAG)
                           ---------------------------------------
                           Σ intensity(all proteins in MAG)

Ribosomal allocation was defined as the ribosomal-protein fraction. Translation-associated
allocation was the sum of the ribosomal, translation-factor, aminoacyl-tRNA-synthetase and
ribosome-biogenesis fractions. Growth-machinery allocation additionally included the
DNA-replication and cell-division fractions. Intensities were never compared against an
"average protein"; the denominator was always the total lineage-specific proteomic
intensity. The fraction of each population's LFQ intensity that fell in KO-annotated
proteins was recorded as an annotation-coverage control (pooled cave-biofilm: SC_MAG_00006
79%, SC_MAG_00008 55%, *Mycobacterium* 89%). Per-sample allocation values and
identification counts are provided in `allocation_by_sample.csv`.

## Interpretation and caveats

Low growth-machinery allocation combined with comparatively high maintenance (defense/
stress) allocation is consistent with a slow, stress-adapted physiological state. This
proxy is not linearly proportional to growth rate: even near-dormant cells retain active
ribosomes for protein turnover, so ribosomal allocation should be read as a
physiological-state indicator rather than a growth-rate measurement. Because LFQ
intensities — unlike iBAQ — are not corrected for the number of theoretically observable
peptides, summed intensities carry a protein length/detectability bias; within-population
category rankings and the cross-population ordering are robust, but the exact percentages
would shift under an iBAQ-based quantification. Cross-population comparison of individual
categories should also account for the differing annotation coverage noted above, and the
*Ferroplasma* SC_MAG_00006 culture values should be disregarded on the identification-
confidence grounds noted above.

## Software

Analysis and figures were produced in R (readxl, readr, dplyr, tidyr, stringr, purrr for
data processing; ggplot2 with ggpattern and patchwork for the figure). KEGG BRITE
hierarchies and ortholog names were retrieved from the KEGG REST API
(https://rest.kegg.jp). The analysis script is `plotting/proteome_allocation_by_sample.R`.
