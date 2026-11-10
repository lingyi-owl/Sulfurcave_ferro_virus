# plotting

One script per manuscript figure, named by figure number (`FigN_*.R` main text,
`FigSNN_*.R` supplementary). The table in the top-level README lists each script and
its output folder under `figures/`. Shared data preparation lives in `analysis/`.

The CSV files here are inputs read by `Fig6_tail_fiber_features.R`
(`fig6_aa_profile.csv`, `fig6_paired_means_acd.csv`) and saved results compared
against by `FigS06_proteome_allocation.R` (`allocation_by_sample.csv`).
`cave_imgvr_features_slim.csv` is no longer read by any script.
