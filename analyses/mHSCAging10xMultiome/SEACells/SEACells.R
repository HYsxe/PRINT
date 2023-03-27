# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(Seurat)
library(Matrix)
library(hdf5r)

# Retrieve LSI embedding of single cells
cellEmbedding <- readRDS("../../../data/mHSCAging10xMultiome/LSIEmbedding.rds")
scATACSeurat <- readRDS("../../../data/mHSCAging10xMultiome/scATACSeurat.rds")
seuratClusters <- scATACSeurat$seurat_clusters

# Only keep HSCs
HSC_filter <- (scATACSeurat$group == "HSC") & (seuratClusters %in% c(0,6,10))
scATACSeurat <- scATACSeurat[,HSC_filter]
cellEmbedding <- cellEmbedding[HSC_filter, ]

# Retrieve metadata
cellMetadata <- scATACSeurat@meta.data

# Get count matrix
countsIJV <- as.matrix(summary(scATACSeurat@assays$ATAC@counts))

# Save to a h5 file
h5_path <- "../../../data/mHSCAging10xMultiome/SEACells.h5"
h5file <- H5File$new(h5_path, mode = "w")
h5file[["cellMetadata"]] <- cellMetadata
h5file[["cellEmbedding"]] <- cellEmbedding
h5file[["counts"]] <- countsIJV
h5file[["UMAP"]] <- scATACSeurat@reductions$umap@cell.embeddings
h5file[["barcodes"]] <- rownames(cellEmbedding)
h5file[["peaks"]] <- rownames(scATACSeurat)
h5file$close_all()