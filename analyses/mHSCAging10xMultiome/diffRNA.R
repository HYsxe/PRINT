# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getGroupData.R")

###################
# Load input data #
###################

# Load footprinting project
project <- readRDS("../../data/mHSCAging10xMultiome/project.rds")

# Load LSI embedding of single cells
cellEmbedding <- readRDS("../../data/mHSCAging10xMultiome/LSIEmbedding.rds")

# Load barcodes for each pseudobulk
barcodeGrouping <- read.table("../../data/mHSCAging10xMultiome/barcodeGrouping.txt",
                              header = T)

# Load pseudobulk centers
pseudobulkCenters <- read.table("../../data/mHSCAging10xMultiome/pseudobulkCenters.txt",
                                header = T)
rownames(pseudobulkCenters) <- pseudobulkCenters$group

# Load single cell RNA data
scRNA <- readRDS("../../data/mHSCAging10xMultiome/scRNA.rds")

# Load barcodes for each pseudo-bulk
barcodeGroups <- read.table("../../data/mHSCAging10xMultiome/barcodeGrouping.txt", header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- gtools::mixedsort(unique(barcodeGroups$group))

#########################
# Re-number pseudobulks #
#########################

pseudobulkNames <- groups(project)
pseudobulkCenters <- pseudobulkCenters[pseudobulkNames,]
pseudobulkAges <- sapply(pseudobulkCenters$barcode, function(x){strsplit(x, "-")[[1]][2]})
pseudobulkLSI <- cellEmbedding[pseudobulkCenters$barcode, 2]

# Within each age-group, re-order pseudobulks by their LSI-2 coordinate
# We don't use the first LSI because it highly correlates with depth
reOrder <- c(order(pseudobulkLSI[pseudobulkAges == "Old"]), 
             order(pseudobulkLSI[pseudobulkAges == "Young"]) + sum(pseudobulkAges == "Old"))

#############################
# Differential RNA analysis #
#############################

RNAMatrix <- scRNA@assays$RNA@data
pseudobulkRNA <- getGroupRNA(RNAMatrix, barcodeGroups)
pseudobulkRNA <- pseudobulkRNA[, reOrder]
colnames(pseudobulkRNA) <- pseudobulkNames

# Get pseudobulk metadata
metadata <- t(sapply(pseudobulkNames, 
                     function(x){
                       strsplit(x, "_")[[1]][c(1,2)]
                     }))
metadata <- as.data.frame(metadata)
colnames(metadata) <- c("Age", "pbulk_ind")

# Generate DDS object
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = pseudobulkRNA,
  colData = metadata,
  design = ~ Age)

# Model fitting
dds <- DESeq2::DESeq(dds)

# Retrieve differential analysis results
diffRNA <- DESeq2::results(dds, contrast = c("Age","Old","Young"))
diffRNA <- diffRNA[!is.na(diffRNA$padj),]
diffRNA <- diffRNA[order(diffRNA$pvalue),]
write.table(diffRNA, "../../data/mHSCAging10xMultiome/diffRNA.tsv",
            sep = "\t", quote = F)
