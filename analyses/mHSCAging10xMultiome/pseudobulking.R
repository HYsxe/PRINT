# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(ggplot2)
source("../../code/getPseudoBulks.R")

###################
# Load input data #
###################

# Load single cell ATAC
scATACSeurat <- readRDS("../../data/mHSCAging10xMultiome/scATACSeurat.rds")
scATACSE <- readRDS("../../data/mHSCAging10xMultiome/scATACSE.rds")
cellEmbedding <- scATACSeurat@reductions$lsi@cell.embeddings[,2:20]
seuratClusters <- scATACSeurat$seurat_clusters
cellType <- scATACSeurat$group
age <- scATACSeurat$age

# Get read count per cell
scATACCounts <- scATACSeurat@assays$ATAC@counts
cellCounts <- colSums(scATACCounts)

# Load SEACells results
SEACells <- read.table("../../data/mHSCAging10xMultiome/SEACells.tsv",
                       sep = "\t")
colnames(SEACells) <- c("barcode", "center")
centerCells <- unique(SEACells$center)
centerCells <- data.frame(
  barcode = centerCells,
  age = sapply(centerCells, function(x){strsplit(x, "-")[[1]][2]})
)

########################
# Pseudobulk the cells #
########################

barcodeGroups <- NULL
pseudobulkCenters <- NULL
countThreshold <- 5e6
for(ageGroup in c("Young", "Old")){
  
  # Select HSC cells in one age group
  mask <- (age == ageGroup) & (cellType == "HSC") & (seuratClusters %in% c(0,6,10))
  centerBarcodes <- centerCells[centerCells$age == ageGroup,]$barcode
  
  # For each center cell, rank all cells in the age group by distance in LSI space
  groupMembers <- FNN::get.knnx(cellEmbedding[mask, ], 
                                cellEmbedding[mask, ][centerBarcodes, ], 
                                k = sum(mask))$nn.index
  
  # Generate pseudobulks by keep adding the next nearest cell from the center cell
  # Stop until we reach a certain total number of reads
  scBarcodes <- rownames(cellEmbedding)[mask]
  barcodeGp <- pbmcapply::pbmclapply(
    1:dim(groupMembers)[1],
    function(groupInd){
      scCounts <- cellCounts[mask][groupMembers[groupInd,]]
      nKeptCell <- min(which(cumsum(scCounts) > countThreshold))
      data.frame(barcode = scBarcodes[groupMembers[groupInd, 1:nKeptCell]],
                 group = paste0(ageGroup, "_", groupInd))
    }
  )
  barcodeGp <- data.table::rbindlist(barcodeGp)
  barcodeGroups <- rbind(barcodeGroups, barcodeGp)
  
  centerBarcodes <- data.frame(
    barcode = centerBarcodes,
    group = paste0(ageGroup, "_", 1:dim(groupMembers)[1])
  )
  pseudobulkCenters <- rbind(pseudobulkCenters, centerBarcodes)
}
groupIDs <- unique(barcodeGroups$group)

# Examine the depth of each pseudobulk
# Check depth of pseudobulks
pseudobulkDepths <- pbmcapply::pbmcmapply(
  function(ID){
    sum(scATACSE[, barcodeGroups$barcode[barcodeGroups$group == ID]]$depth) / 1e6
  },
  groupIDs,
  mc.cores = 8
)

system("mkdir ../../data/mHSCAging10xMultiome/heterogeneity")
write.table(barcodeGroups, "../../data/mHSCAging10xMultiome/barcodeGrouping.txt",
            row.names = F, sep = "\t", quote = F)
write.table(pseudobulkCenters, "../../data/mHSCAging10xMultiome/pseudobulkCenters.txt",
            row.names = F, sep = "\t", quote = F)

