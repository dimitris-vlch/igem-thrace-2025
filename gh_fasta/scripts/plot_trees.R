#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# plot_trees.R
# Load existing distance matrices (<basename>_matrix.tsv) and plot both:
#   - Hierarchical clustering dendrogram (complete linkage)
#   - Neighbor-Joining phylogenetic tree
#
# Usage:
#   Rscript plot_trees.R
#
# Output:
#   Plots shown in R graphics device
#   (optionally save to PDF if you uncomment the pdf() block)
#
# Authors: Sotiris & Dimitris
# ------------------------------------------------------------------------------

library(ape)

tsv_files <- list.files(pattern = "_matrix.tsv$")

if (length(tsv_files) == 0) {
  stop("No _matrix.tsv files found in working directory: ", getwd())
}

# Uncomment these lines if you want to save everything in one PDF:
# pdf("all_trees.pdf")

for (f in tsv_files) {
  base <- tools::file_path_sans_ext(basename(f))
  m <- as.matrix(read.table(f, header = TRUE, row.names = 1))
  d <- as.dist(m)

  # Dendrogram
  hc <- hclust(d, method = "complete")
  plot(hc, hang = -1, cex = 0.5,
       main = paste("Hierarchical Clustering:", base))

  # Neighbor-Joining tree
  tree <- nj(d)
  plot(tree, main = paste("Neighbor-Joining Tree:", base))
}

# dev.off()
