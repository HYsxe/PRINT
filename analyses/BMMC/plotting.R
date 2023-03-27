# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getAggregateFootprint.R")
source("../../code/getSubstructures.R")
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "BMMC"
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")

# Load the footprintingProject object
project <- readRDS(paste0(projectDataDir, "footprintingProject.rds"))

# Load pseudo-bulk metadata
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
groupUMAP(project) <- as.matrix(cbind(groupInfo$UMAP1, groupInfo$UMAP2))
groupCellType(project) <- groupInfo$cellType
groupPseudoTime(project) <- groupInfo$Pseudotime

# Get color mapping for cell types
cellTypes <- unique(groupInfo$cellType)
cellTypeColors <- groupInfo$color[match(cellTypes, groupInfo$cellType)]
names(cellTypeColors) <- cellTypes

# Load PWM data
motifs <- readRDS(paste0(projectDataDir, "filteredMotifs.rds"))

# Load TFBS prediction model
h5Path <- "../../data/TFBSPrediction/TFBS_model.h5"
TFBSModel <- keras::load_model_hdf5(h5Path)
TFBSModel <- keras::get_weights(TFBSModel)
h5file <- H5File$new(h5Path, mode="r")
footprintMean <- h5file[["footprint_mean"]]
footprintMean <- footprintMean[1:footprintMean$dims]
footprintSd <- h5file[["footprint_sd"]]
footprintSd <- footprintSd[1:footprintSd$dims]
h5file$close_all()
TFBSModel$footprintMean <- footprintMean
TFBSModel$footprintSd <- footprintSd
TFBSModel$scales <- c(10,20,30,50,80,100)
TFBindingModel(project) <- TFBSModel

################################
# Plot multi-kernel footprints #
################################

regionInd <- 100727
system("mkdir ../../data/BMMC/plots/CREMultiScaleFootprints")
pdf(paste0("../../data/BMMC/plots/CREMultiScaleFootprints/Peak_", regionInd, "_multiscale.pdf"),
    width = 8, height = 5)
plotKernelScanning(project,
                   regionInd = 100727,
                   lineageGroups = which(groupCellType(project) %in% "HSC/MPP"))
dev.off()

#################################################################
# Plot footprint/TF binding/Tn5 insertion across all cell types #
#################################################################

regionInd <- 178638
system("mkdir ../../data/BMMC/plots/TFBSHeatmap/")
cellTypeOrder <- c("NK", "CD4", "CD8",
                   "CLP", "LMPP", "HSC/MPP", "CMP", "MEP", "early-Ery", "late-Ery",
                   "GMP", "CD16mono", "CD14mono", "pDC", "DC", "Baso",
                   "NaiveB", "plasmaB", "MemoryB", "pro/pre-B")
rowOrder <- order(match(groupCellType(project), cellTypeOrder))

pdf(paste0("../../data/BMMC/plots/TFBSHeatmap/region_", regionInd, "_TFBSHeatmap.pdf"),
    width = 6, height = 5)
plotFeatureHeatmap(project,
                   regionInd = regionInd,
                   feature = "TFBS",
                   cellTypeColors = cellTypeColors,
                   cellTypeOrder = cellTypeOrder,
                   rowOrder = rowOrder)
dev.off()

# Visualize Tn5 insertion
plotFeatureHeatmap(project, 
                   regionInd = regionInd, 
                   feature = "insertion",
                   cellTypeColors = cellTypeColors,
                   cellTypeOrder = cellTypeOrder,
                   rowOrder = rowOrder)

# Visualize nucleosome positions across differentiation
plotFeatureHeatmap(project, 
                   regionInd = regionInd, 
                   feature = "footprint",
                   footprintRadius = 50,
                   cellTypeColors = cellTypeColors,
                   cellTypeOrder = cellTypeOrder,
                   rowOrder = rowOrder)

###############################################################
# Plot footprint dynamics along a lineage for a specific region #
###############################################################

LineageInd <- 2
lineageProbs <- groupInfo[, c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")]
lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                            function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})
lineageGroups <- which(lineageIdentities %in% LineageInd)

regionInd <- 9049
system("mkdir ../../data/BMMC/plots/nucleosomeTracks/")
pdf(paste0("../../data/BMMC/plots/nucleosomeTracks/region_", regionInd, "_TFBSHeatmap.pdf"),
    width = 6, height = 5)
plotPseudotimeHeatmap(
  project,
  regionInd = regionInd,
  feature = "TFBS",
  lineageGroups = lineageGroups,
  cellTypeColors = cellTypeColors
)
dev.off()

pdf(paste0("../../data/BMMC/plots/nucleosomeTracks/region_", regionInd, "_nucleosomeTracks.pdf"),
    width = 6, height = 5)
plotPseudotimeHeatmap(
  project,
  regionInd = regionInd,
  feature = "footprint",
  lineageGroups = lineageGroups,
  cellTypeColors = cellTypeColors,
  heatmapPalette = colorRampPalette(c("white",  "#9ECAE1", "#08519C", "#08306B"))(9),
  footprintRadius = 50
)
dev.off()

#########################
# Plot CRE segmentation #
#########################

regionInd <-  55465
system("mkdir ../../data/BMMC/plots/CRESegmentation/")
pdf(paste0("../../data/BMMC/plots/CRESegmentation/region_", regionInd, "_TFBSHeatmap.pdf"),
    width = 6, height = 5)
plotMultiScale(project, 
               regionInd = regionInd)
dev.off()
