# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(SummarizedExperiment)

#################################
# 1. Get genomic ranges of CREs #
#################################

# Load single cell ATAC
scATAC <- readRDS("../../data/BMMCCellType/atac.se.rds")

# Get region ranges
regions <- rowRanges(scATAC)
regions <- IRanges::resize(regions, width = 1000, fix = "center")
saveRDS(regions, "../../data/BMMCCellType/regionRanges.rds")

###############################################
# 2. Only keep motifs enriched in CRE regions #
###############################################

# Load PWM data
motifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))
filtMotifs <- filterMotifs(regions = regions, 
                           motifs = motifs,
                           genome = "hg38")
saveRDS(motifs, "../../data/BMMC/filteredMotifs.rds")

############################
# 3. Generate pseudo-bulks #
############################

# Reformat cell barcodes
scBarcodes <- colnames(scATAC)

# Pseudobulk each cell type
cellTypeLabels <- scATAC$cistopic.assign.l2.rank
cellTypes <- sort(unique(cellTypeLabels))
groups <- 1:length(cellTypes)
names(groups) <- cellTypes
barcodeGroups <- data.frame(
  barcode = scBarcodes,
  group = unname(groups[cellTypeLabels])
)
write.table(barcodeGroups, "../../data/BMMCCellType/barcodeGrouping.txt",
            row.names = F, sep = "\t", quote = F)

# Also save metadata for each pseudobulk
groupInfo <- data.frame(
  groupID = unname(groups),
  cellType = names(groups)
)
write.table(groupInfo, "../../data/BMMCCellType/groupInfo.txt",
            row.names = F, sep = "\t", quote = F)
