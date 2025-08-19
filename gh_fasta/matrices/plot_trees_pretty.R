#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# plot_trees_pretty_labels.R
# Pretty dendrograms + NJ trees from existing *_matrix.tsv,
# enriching tip labels using the matching raw TSV (e.g. gh80.tsv).
#
# For each <base>_matrix.tsv:
#   - Loads distance matrix
#   - Loads <base>.tsv to map accession -> (accession | species/organism)
#   - Writes:
#       <base>_dendrogram.pdf
#       <base>_nj.pdf
#
# Usage:
#   Rscript plot_trees_pretty_labels.R
#
# Authors: Sotiris & Dimitris
# Repository: https://github.com/dimitris-vlch/igem-thrace-2025
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ape)
})

# ---------- helpers ----------------------------------------------------

# Try to read matrix robustly (handles weird headers, preserves names)
load_matrix <- function(path) {
  try1 <- try(read.table(path, header = TRUE, row.names = 1,
                         sep = "\t", check.names = FALSE,
                         quote = "", comment.char = "", fill = TRUE),
              silent = TRUE)
  if (!inherits(try1, "try-error")) {
    df <- try1
  } else {
    try2 <- try(read.table(path, header = TRUE,
                           sep = "\t", check.names = FALSE,
                           quote = "", comment.char = "", fill = TRUE),
                silent = TRUE)
    if (inherits(try2, "try-error")) {
      stop("Failed to read matrix: ", path, "\n", attr(try2, "condition")$message)
    }
    df <- try2
    rn <- df[[1]]
    df <- df[, -1, drop = FALSE]
    rownames(df) <- rn
  }

  if (nrow(df) != ncol(df)) {
    common <- intersect(rownames(df), colnames(df))
    df <- df[common, common, drop = FALSE]
  }

  for (j in seq_len(ncol(df))) df[[j]] <- suppressWarnings(as.numeric(df[[j]]))
  m <- as.matrix(df)

  m <- (m + t(m)) / 2
  diag(m) <- 0

  keep <- rowSums(is.na(m)) == 0 & colSums(is.na(m)) == 0
  m <- m[keep, keep, drop = FALSE]

  if (nrow(m) < 2) stop("Matrix too small after cleaning: ", path)
  m
}

extract_accession <- function(x) {
  parts <- unlist(strsplit(x, "[| ]+"))
  parts <- parts[nzchar(parts)]
  parts <- sub("\\.[0-9]+$", "", parts)
  hit <- parts[grepl("^[A-Za-z]{1,4}[0-9]{4,}$", parts)]
  if (length(hit)) return(hit[1])
  hit2 <- parts[grepl("(?=.*[A-Za-z])(?=.*[0-9])", parts, perl = TRUE)]
  if (length(hit2)) return(hit2[1])
  parts[1]
}

build_mapping <- function(meta_path) {
  if (!file.exists(meta_path)) return(NULL)
  meta <- try(read.table(meta_path, header = TRUE, sep = "\t",
                         quote = "", comment.char = "", check.names = FALSE, fill = TRUE),
              silent = TRUE)
  if (inherits(meta, "try-error")) return(NULL)

  cn <- tolower(colnames(meta))

  acc_idx <- which(cn %in% c("accession", "run_accession", "sample_accession", "study_accession"))
  if (!length(acc_idx)) {
    if (grepl("accession", cn[1])) acc_idx <- 1 else return(NULL)
  }

  acc <- sub("\\.[0-9]+$", "", as.character(meta[[acc_idx[1]]]))

  org_idx <- which(cn %in% c("scientific_name", "organism", "organism_name", "tax_name", "species"))
  org <- rep(NA_character_, length(acc))
  if (length(org_idx)) {
    org <- as.character(meta[[org_idx[1]]])
  } else {
    tx_idx <- which(grepl("tax", cn))
    if (length(tx_idx)) org <- as.character(meta[[tx_idx[1]]])
  }

  label <- ifelse(!is.na(org) & nzchar(org), paste0(acc, " | ", org), acc)
  mp <- tapply(label, acc, `[`, 1)
  return(mp)
}

make_labels <- function(raw_names, mapping = NULL, max_chars = 45) {
  accs <- vapply(raw_names, extract_accession, character(1))
  if (!is.null(mapping)) {
    lbl <- ifelse(accs %in% names(mapping), mapping[accs], accs)
  } else {
    lbl <- accs
  }
  too_long <- nchar(lbl) > max_chars
  lbl[too_long] <- paste0(substr(lbl[too_long], 1, max_chars - 1), "…")
  as.character(lbl)
}

# ---------- plotting ---------------------------------------------------

# ---------- plotting ---------------------------------------------------

# Dendrogram: keep as is (wider canvas, trimmed left, extra right space)
plot_dendro_pdf <- function(d, labels, title, file) {
  hc <- hclust(d, method = "complete"); hc$labels <- labels
  dend <- as.dendrogram(hc)

  n   <- length(labels)
  cex <- max(min(0.8 * (60 / max(60, n)), 0.8), 0.25)
  w   <- 18                      # wider canvas for long labels
  h   <- max(8, n * 0.18)

  pdf(file, width = w, height = h)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)
  par(mar = c(4, 6, 3, 10))      # bottom, left, top, right
  plot(dend, horiz = TRUE, leaflab = "perpendicular",
       main = paste("Hierarchical clustering:", title), cex = cex)
}

# Neighbor-Joining: robust scale bar placement (never clipped)
plot_nj_pdf <- function(d, labels, title, file) {
  tr <- nj(d); tr$tip.label <- labels; tr <- ladderize(tr)

  n   <- length(labels)
  cex <- max(min(0.8 * (60 / max(60, n)), 0.8), 0.25)
  w   <- 18
  h   <- max(8, n * 0.18)

  pdf(file, width = w, height = h)
  op <- par(no.readonly = TRUE); on.exit({par(op); dev.off()}, add = TRUE)

  # give some breathing room around the plot
  par(mar = c(8, 4, 3, 10))  # extra bottom helps with long trees

  plot(tr, type = "phylogram", direction = "rightwards",
       cex = cex, no.margin = FALSE,
       main = paste("Neighbor-Joining tree:", title))

  # Place the scale bar *inside* the plotting region so the PDF viewer can't clip it.
  usr <- par("usr")                        # c(xmin, xmax, ymin, ymax)
  dx  <- usr[2] - usr[1]
  dy  <- usr[4] - usr[3]
  x0  <- usr[1] + 0.05 * dx                # 5% from the left
  y0  <- usr[3] + 0.05 * dy                # 5% above the bottom
  add.scale.bar(x = x0, y = y0, cex = 0.9, lwd = 1.2)
}

# ---------- main -------------------------------------------------------

matrix_files <- list.files(pattern = "_matrix\\.tsv$")
if (!length(matrix_files)) stop("No _matrix.tsv files found in: ", getwd())

for (mf in matrix_files) {
  base <- sub("_matrix\\.tsv$", "", basename(mf))
  message(">> Plotting ", mf, " with metadata from ", paste0(base, ".tsv"))

  m <- load_matrix(mf)
  d <- as.dist(m)

  mapping <- build_mapping(paste0(base, ".tsv"))
  labs <- make_labels(rownames(m), mapping, max_chars = 45)

  plot_dendro_pdf(d, labs, paste0(base, "_matrix"), paste0(base, "_matrix_dendrogram.pdf"))
  plot_nj_pdf(d, labs, paste0(base, "_matrix"),    paste0(base, "_matrix_nj.pdf"))
}
message("Done.")
