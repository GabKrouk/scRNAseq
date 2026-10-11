# Step 65: raw counts of the Shahan et al. (2022) WT root atlas (GEO GSE152766,
# GSE152766_Root_Atlas_spliced_unspliced_raw_counts.rds.gz) -> total raw UMI counts (spliced + unspliced),
# written as binary CSC components for Python (genes x cells). No normalisation, no imputation:
# raw counts only (normalised/SCT matrices impute low-expressed TFs, see the scTARGET QC lesson).
# Data live outside Dropbox: ~/scdata/Shahan2022/
suppressPackageStartupMessages(library(Matrix))
A <- path.expand("~/scdata/Shahan2022")
OUT <- file.path(A, "derived"); dir.create(OUT, showWarnings = FALSE)
x <- readRDS(file.path(A, "raw_counts", "Root_Atlas_spliced_unspliced_raw_counts.rds"))   # GEO file is gzipped twice: gzip -dc once first
# The file is a Seurat (v3) object with two assays (spliced_RNA, unspliced_RNA). Seurat is not needed:
# the raw count matrices are read from the S4 slots directly.
as_ <- attr(x, "assays"); cat("assays:", names(as_), "\n")
mats <- lapply(as_, function(a) attr(a, "counts"))
for (n in names(mats)) cat(n, class(mats[[n]]), dim(mats[[n]]), "nnz", length(mats[[n]]@x), "\n")
g <- rownames(mats[[1]]); cl <- colnames(mats[[1]])
for (m in mats) stopifnot(identical(rownames(m), g), identical(colnames(m), cl))
md <- attr(x, "meta.data")
data.table::fwrite(data.frame(cell = cl, nCount_spliced = colSums(mats$spliced_RNA), nCount_unspliced = colSums(mats$unspliced_RNA)),
                   file.path(OUT, "cells_counts.tsv"), sep = "\t")
M <- Reduce(`+`, mats)                       # spliced + unspliced (+ ambiguous if present) = total raw UMIs
M <- as(M, "dgCMatrix")
cat("total matrix:", dim(M), " nnz:", length(M@x), " integer counts:", all(M@x == round(M@x)), "\n")
writeBin(as.integer(M@i), file.path(OUT, "counts_i.int32"), size = 4)
writeBin(as.integer(M@p), file.path(OUT, "counts_p.int32"), size = 4)
writeBin(as.integer(round(M@x)), file.path(OUT, "counts_x.int32"), size = 4)
writeLines(g, file.path(OUT, "genes.txt")); writeLines(cl, file.path(OUT, "cells.txt"))
writeLines(as.character(dim(M)), file.path(OUT, "dim.txt"))
cat("first genes:", head(g), "\nfirst cells:", head(cl), "\n")
