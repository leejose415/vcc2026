# Export a Jiang / Satija pathway Perturb-seq Seurat object's full raw RNA counts without the Seurat package.
# The object is read with readRDS and its slots through attr() (Seurat is not installed). The RNA counts are a dgCMatrix
# (genes x cells); stored column by column it IS the CSR of cells x genes, so x / i / p are written as they are.
#   Rscript scripts/extract_jiang_counts.R <Seurat_object_X_Perturb_seq.rds> <outdir>
# Writes <outdir>/x.bin (int32 counts), i.bin (int32 gene index), p.bin (int64 cell pointers), genes.txt, metadata.csv.gz,
# shape.txt (cells genes nnz), complete.txt. Every count must be a non-negative integer (stops otherwise).
suppressPackageStartupMessages(library(Matrix))
a <- commandArgs(trailingOnly = TRUE); src <- a[1]; out <- a[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cat('read', src, '\n'); flush.console()
x <- readRDS(src)
m <- attr(x, 'meta.data'); assay <- attr(x, 'assays')$RNA
C <- attr(assay, 'counts')
if (is.null(C) || !inherits(C, 'dgCMatrix')) {                 # Seurat v5 Assay5: counts in layers
  L <- attr(assay, 'layers'); C <- if (!is.null(L)) L$counts else NULL
}
stopifnot(inherits(C, 'dgCMatrix'), ncol(C) == nrow(m), identical(colnames(C), rownames(m)))
genes <- rownames(C)
rm(x, assay); invisible(gc())
for (s in seq.int(1, length(C@x), by = 50000000)) {
  v <- C@x[s:min(s + 49999999, length(C@x))]
  stopifnot(all(is.finite(v)), all(v >= 0), all(v == round(v)))
}
rm(v); invisible(gc())
writeBin(as.integer(C@x), file.path(out, 'x.bin'), size = 4, endian = 'little')
writeBin(as.integer(C@i), file.path(out, 'i.bin'), size = 4, endian = 'little')
writeBin(as.double(C@p), file.path(out, 'p.bin'), size = 8, endian = 'little')   # pointers as float64 (exact below 2^53)
writeLines(genes, file.path(out, 'genes.txt'))
writeLines(as.character(c(ncol(C), nrow(C), length(C@x))), file.path(out, 'shape.txt'))
write.csv(m, gzfile(file.path(out, 'metadata.csv.gz')), row.names = TRUE)
writeLines('ok', file.path(out, 'complete.txt'))
cat('exported', ncol(C), 'cells x', nrow(C), 'genes,', length(C@x), 'stored counts\n')
