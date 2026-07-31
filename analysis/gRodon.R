#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# run_gRodon.R  —  predict maximal growth rate from Prokka annotations
#
# Input : a Prokka *.ffn file (nucleotide FASTA of CDS + non-coding RNA).
#         gRodon auto-filters the non-coding transcripts (rRNA/tRNA) by
#         length, and identifies "highly expressed" genes from the
#         ribosomal-protein product names in the FASTA headers.
#
# Output: a one-row-per-genome CSV with predicted doubling time (hours),
#         confidence interval, and the codon-usage diagnostics.
# ------------------------------------------------------------------

## ---- 0. One-time install (uncomment to run once) -----------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "coRdon"))
install.packages(c("matrixStats", "devtools"))
devtools::install_github("jlw-ecoevo/gRodon2")   # v2: adds metagenome_v2 + eukaryote support

suppressPackageStartupMessages({
  library(gRodon)
  library(Biostrings)
})

# ================================================================
# CASE A — a single genome / MAG
# ================================================================

ffn_path <- "/Users/wu000058/Library/Mobile Documents/com~apple~CloudDocs/Projects/SulfurCave/sulfur_cave_data/complete_genomes/Ferroplasma_complete_genome/Ferroplasma_Prokka/PROKKA_04262022.ffn"           # <-- your Prokka .ffn file

genes <- readDNAStringSet(ffn_path)

# Flag ribosomal proteins as the highly-expressed set, straight from
# the Prokka product annotations carried in the FASTA headers.
highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)

message("Total transcripts: ", length(genes),
        " | ribosomal-protein genes flagged: ", sum(highly_expressed))
# gRodon guidance: predictions get biased when < 10 highly-expressed
# genes are available, so treat low counts with caution.

# Pick ONE mode for the organism you have:
#   "full"          complete isolate genome (default)
#   "partial"       incomplete genome / MAG / SAG (varying completeness)
#   "metagenome_v2" whole-community assembly, GC-bias-corrected (gRodon2)
#   "metagenome"    community assembly, original v1 model
res <- predictGrowth(genes, highly_expressed, mode = "partial")

# Optional: correct for growth temperature of a mesophile (NOT validated
# on extremophiles — use with care):
# res <- predictGrowth(genes, highly_expressed, mode = "full", temperature = 37)

print(res)
# res$d       -> predicted minimum doubling time (hours)
# res$LowerCI / res$UpperCI -> 95% CI on the doubling time
# res$nHE     -> number of highly-expressed genes actually used
# res$CUBHE, res$ConsistencyHE, res$CPB, res$FilteredSequences -> diagnostics

as.data.frame(res) |>
  transform(genome = basename(ffn_path)) |>
  write.csv("gRodon_result.csv", row.names = FALSE)


# ================================================================
# CASE B — a metagenome assembly (community-average growth potential)
# ================================================================
# genes <- readDNAStringSet("metagenome_scaffolds.ffn")
# highly_expressed <- grepl("ribosomal protein", names(genes), ignore.case = TRUE)
#
# # Unweighted community average:
# predictGrowth(genes, highly_expressed, mode = "metagenome_v2")
#
# # Abundance-weighted (recommended): supply mean depth of coverage per ORF,
# # as a numeric vector aligned to `genes` (relative values are fine).
# # cov <- <named/ordered vector of per-ORF mean coverage>
# # predictGrowth(genes, highly_expressed, mode = "metagenome_v2",
# #               depth_of_coverage = cov)


# ================================================================
# CASE C — many genomes / MAGs in a folder (batch)
# ================================================================
run_one <- function(ffn, mode = "full") {
  g  <- readDNAStringSet(ffn)
  he <- grepl("ribosomal protein", names(g), ignore.case = TRUE)
  out <- tryCatch(as.data.frame(predictGrowth(g, he, mode = mode)),
                  error = function(e) data.frame(d = NA, note = conditionMessage(e)))
  out$genome  <- basename(ffn)
  out$n_ribo  <- sum(he)
  out
}

# ffn_files <- list.files("prokka_out", pattern = "\\.ffn$",
#                         full.names = TRUE, recursive = TRUE)
# batch <- do.call(rbind, lapply(ffn_files, run_one, mode = "partial"))  # MAGs -> "partial"
# write.csv(batch, "gRodon_batch_results.csv", row.names = FALSE)
