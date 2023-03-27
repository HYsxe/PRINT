# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
library(Gviz)

#################################################
# Load bone marrow data needed for footprinting #
#################################################

# Initialize a footprintingProject object
projectName <- "BMMCCellType"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regions <- readRDS(paste0(projectDataDir, "regionRanges.rds"))

# Set the regionRanges slot
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getRegionBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load background dispersion model
dispModel(project, "20bp") <- readRDS("../../data/shared/dispModel/dispersionModel20bp.rds")

# Load cell type labels for pseudobulks
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
groups(project) <- groupInfo$groupID
cellTypes <- unique(groupInfo$cellType)

############################################################################
# Find CREs overlapping with BAC regions and extract BAC tagmentation data #
############################################################################

# Load BAC ranges
selectedBACs <- read.table("../../data/BAC/rawData/selectedBACs.txt")
colnames(selectedBACs) <- c("ID", "chr", "start", "end")
BACRanges <- GRanges(
  seqnames = selectedBACs$chr,
  ranges = IRanges(start = selectedBACs$start,
                   end = selectedBACs$end)
)
names(BACRanges) <- selectedBACs$ID

# Load countTensor for BAC tagmentation data
backgroundCounts <- readRDS("../../data/BAC/countTensors.rds")

# Get ATAC insertion tracks for each BAC
ATACTrackList <- list()
nReplicates <- 5
for(BACInd in 1:length(BACRanges)){
  
  # Get the countTensor for the current BAC
  BACCounts <- backgroundCounts[["all"]][[BACInd]]
  
  # Go through each replicate and extract the ATAC insertion track for the current BAC
  ATACTrackList[[BACInd]] <- lapply(
    1:nReplicates,
    function(repInd){
      BACRepCounts <- BACCounts %>% filter(group %in% repInd) 
      BACRepATAC <- rep(0, width(BACRanges)[BACInd])
      BACRepATAC[BACRepCounts$position] <- BACRepCounts$count
      BACRepATAC
    }
  )
  names(ATACTrackList[[BACInd]]) <- paste0("Observed, replicate ", 1:nReplicates)
  
  # Combine the data for all replicates
  ATACTrackList[[BACInd]][["Observed, all Replicates"]] <- Reduce("+", ATACTrackList[[BACInd]])
}

# Remove low coverage BACs
keptBACs <- sapply(ATACTrackList, function(x){mean(x[["Observed, all Replicates"]]) > 20})
BACRanges <- BACRanges[keptBACs]
ATACTrackList <- ATACTrackList[keptBACs]

# Find region-BAC overlaps
ovPairs <- findOverlaps(regions, BACRanges, type = "within")

# Get color mapping for cell types
cellTypeColorFile <- read.table("../../data/shared/cellTypeColors.txt", 
                                sep = ",", comment.char = "")
cellTypeColors <- cellTypeColorFile$V2
names(cellTypeColors) <- cellTypeColorFile$V1

##############################################################################
# Select a specific region overlapping a BAC. Retrieve ATAC and footprint data #
##############################################################################

# Calculate coverage for each overlapping region
coverage <- pbmcapply::pbmcmapply(
  function(ovInd){
    ovRegionInd <- ovPairs@from[ovInd]
    ret <- getCountData(project, ovRegionInd)
    countData <- ret[["countData"]]
    chunkRegionInd <- ret[["regionInd"]]
    sum(countData[[chunkRegionInd]]$count)
  },
  1:length(ovPairs),
  mc.cores = 16
)

# Select a region-BAC pair
ovInd <- 322#321#153
genome <- "hg38"
context_radius <- 50
context_len <- 2 * context_radius + 1

# Get the indices and ranges of the region and the overlapping BAC
ovRegionInd <- ovPairs@from[ovInd]
ovRegion <- regions[ovRegionInd] 
ovBACInd <- ovPairs@to[ovInd]
ovBAC <- BACRanges[ovBACInd]
relativePos <- start(ovRegion) - start(ovBAC) + 1
width <- width(ovRegion)
chr <- as.character(seqnames(ovRegion))

# For the selected region region get observed Tn5 insertion counts at each position in our BAC data
relativeRange <- relativePos:(relativePos + width(ovRegion) - 1)
BACInsertionTrack <- ATACTrackList[[ovBACInd]][["Observed, all Replicates"]][relativeRange]
BACInsertionMatrix <- sapply(
  1:nReplicates,
  function(i){
    ATACTrackList[[ovBACInd]][[paste0("Observed, replicate ", i)]][relativeRange]
  }
)

# Calculate ground truth Tn5 insertion bias using BAC data
obsBias <- BACInsertionTrack / conv(BACInsertionTrack, context_radius) * (context_radius * 2)

# Get model-predicted Tn5 bias
predBias <- regionBias(project)[ovRegionInd, ]

# For the selected region region, get predicted bias and observed cutting 
biasTracks <- list("predicted" = predBias,
                   "uncorrected" = rep(1, width),
                   "observed" = obsBias)

# Get position-by-pseudobulk matrix of ATAC signal for the current region
ret <- getCountData(project, ovRegionInd)
countData <- ret[["countData"]]
chunkRegionInd <- ret[["regionInd"]]
groupIDs <- groups(project)
BMMCInsertionMatrix <- getRegionATAC(countData, chunkRegionInd, groupIDs, width)

# Call footprints on bone marrow data with Tn5 bias generated using different methods
BMMCFootprints <- pbmcapply::pbmclapply(
  biasTracks,
  function(biasTrack){
    regionFootprintScoring(
      regionATAC = BMMCInsertionMatrix,
      Tn5Bias = biasTrack,
      dispersionModel = dispModel(project, "20bp"),
      cellTypeLabels = rep("BMMC", dim(BMMCInsertionMatrix)[2]),
      footprintRadius = 20, 
      flankRadius = 20,
      returnCellTypeScores = T,
      scoreType = "pval"
    )$cellTypeScores[, 1]
  },
  mc.cores = 8
)

# Call footprints on BAC data with Tn5 bias generated using different methods
BACFootprints <- pbmcapply::pbmclapply(
  biasTracks,
  function(biasTrack){
    regionFootprintScoring(
      regionATAC = BACInsertionMatrix,
      Tn5Bias = biasTrack,
      dispersionModel = dispModel(project, "20bp"),
      cellTypeLabels = rep("replicate", 5),
      footprintRadius = 20, 
      flankRadius = 20,
      returnCellTypeScores = T,
      scoreType = "pval"
    )$cellTypeScores[, 1]
  },
  mc.cores = 8
)

#####################
# Visualize results #
#####################

system("mkdir -p ../../data/BMMCCellType/plots/footprintTracks/")

maxY <- max(c(BACFootprints$predicted, BACFootprints$uncorrected))
plotRadius <- as.integer(width / 2)
positions <- (-plotRadius + 1):plotRadius
plotData <- data.frame(
  "Position" = (-plotRadius + 1):plotRadius,
  "BACCorrected" = pmin(BACFootprints$predicted, maxY),
  "BACUncorrected" = pmin(BACFootprints$uncorrected, maxY)
)

barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = BACInsertionTrack)
p1 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("BAC Tn5\ninsertion") +
  theme_classic() 

barData <- data.frame(x1 = positions - 2, x2 = positions + 2,
                      y1 = 0, y2 = obsBias)
p2 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Observed\nTn5 bias") +
  theme_classic()

barData <- data.frame(x1 = positions - 2, x2 = positions + 2,
                      y1 = 0, y2 = predBias)
p3 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Predicted\nTn5 bias") +
  theme_classic()

p4 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "BACCorrected"), 
              fill = "#69ACD5") + xlab("") + ylab("BAC corrected\nfootprints") +
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) +
  theme_classic()

p5 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "BACUncorrected"), 
              fill = "#69ACD5") + xlab("") + ylab("BAC uncorrected\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

library(patchwork)
pdf(paste0("../../data/BMMCCellType//plots/footprintTracks/BAC_footprints_", ovRegionInd, ".pdf"),
    width=15, height=8)
p1 / p2 / p3 / p4 / p5 
dev.off()

maxY <- max(BMMCFootprints$predicted)
barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = rowSums(BMMCInsertionMatrix))
plotBMTn5 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#006838") +
  scale_y_continuous(expand = c(0, 0)) + xlab("") + ylab("Bone marrow\nTn5 insertion") +
  theme_classic() 

plotPred <- ggplot() + 
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, BMMCFootprints$predicted)), 
              alpha = 0.5)
plotPred <- plotPred + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("Bone marrow\nmodel\ncorrected footprints") +
  xlab("") +
  theme_classic()

plotObs <- ggplot() + 
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, BMMCFootprints$observed)), 
              alpha = 0.5)
plotObs <- plotObs + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("Bone marrow\nground truth bias\ncorrected footprints") +
  xlab("") +
  theme_classic()

plotUncorrected <- ggplot() +
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, BMMCFootprints$uncorrected)), 
              alpha = 0.5)
plotUncorrected <- plotUncorrected + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("Bone marrow\nuncorrected footprints") + 
  xlab("") +
  theme_classic()

pdf(paste0("../../data/BMMCCellType//plots/footprintTracks/BMMC_footprints_", ovRegionInd, ".pdf"),
    width=15, height=8)
plotBMTn5 /  plotPred / plotObs / plotUncorrected 
dev.off()