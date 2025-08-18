#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# gh_fasta_msa.R
# Perform multiple sequence alignment (MSA) and distance matrix calculation
# for all GH family FASTA files in the current directory.
#
# Workflow:
#   1) Detect all *.fasta files in current directory
#   2) Read sequences (DNA/RNA/AA) using Biostrings
#   3) Run MSA (msa package)
#   4) Convert alignment for distance calculation (seqinr)
#   5) Compute pairwise identity distance matrix
#   6) Save results as <basename>_matrix.tsv
#   7) Plot dendrogram (hierarchical clustering, complete linkage)
#   8) Optionally plot .dnd trees if they exist
#
# Output:
#   gh18_matrix.tsv, gh19_matrix.tsv, ...
#   (plots are shown in R graphics window by default)
#
# Authors: Sotiris & Dimitris
# Repository: https://github.com/dimitris-vlch/igem-thrace-2025
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Biostrings)
  library(msa)
  library(seqinr)
  library(ape)
})

# --- Choose reader function (set according to your data type) ----------
# Default here: DNA
read_fn <- function(path) Biostrings::readDNAStringSet(path)
# For RNA:  read_fn <- function(path) Biostrings::readRNAStringSet(path)
# For AA:   read_fn <- function(path) Biostrings::readAAStringSet(path)

# --- Locate all FASTA files --------------------------------------------
fasta_files <- list.files(pattern = "\\.fasta$", ignore.case = TRUE)
if (length(fasta_files) == 0) {
  stop("No FASTA files found in working directory: ", getwd())
}

# --- Helper function to process one FASTA ------------------------------
process_fasta <- function(fasta_path) {
  message(">> Processing ", fasta_path)

  if (!file.exists(fasta_path)) {
    warning("File not found: ", fasta_path)
    return(NULL)
  }

  # Load sequences
  hemoSeq <- read_fn(fasta_path)

  # Run multiple sequence alignment
  hemoAln  <- msa::msa(hemoSeq)

  # Convert to seqinr format
  hemoAln2 <- msa::msaConvert(hemoAln, type = "seqinr::alignment")

  # Distance matrix
  d <- seqinr::dist.alignment(hemoAln2, "identity")
  m <- as.matrix(d)
  m <- m[rowSums(is.na(m)) == 0, colSums(is.na(m)) == 0, drop = FALSE]
  d_clean <- as.dist(m)

  # Save distance matrix
  base <- tools::file_path_sans_ext(basename(fasta_path))
  out_file <- paste0(base, "_matrix.tsv")
  write.table(as.matrix(d_clean), file = out_file,
              sep = "\t", quote = FALSE, col.names = NA)
  message("   -> Wrote ", out_file)

  # Plot clustering dendrogram
  hc <- hclust(d_clean, method = "complete")
  plot(hc, hang = -1, cex = 0.4,
       main = paste("Hierarchical clustering:", basename(fasta_path)))
}

# --- Main loop over all FASTA files -----------------------------------
for (f in fasta_files) {
  process_fasta(f)
}

# --- Optional: plot precomputed trees ---------------------------------
if (file.exists("sequences.dnd")) {
  tr <- ape::read.tree("sequences.dnd")
  plot(tr, main = "Phylogenetic Tree (sequences.dnd)")
}
if (file.exists("sequences2.dnd")) {
  tr2 <- ape::read.tree("sequences2.dnd")
  plot(tr2, main = "Phylogenetic Tree (sequences2.dnd)")
}
