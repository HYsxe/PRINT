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
projectName <- "HepG2"
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

for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Load barcodes for each pseudo-bulk
scATAC <- readRDS(paste0(projectDataDir, "scATAC.rds"))
barcodeGroups <- data.frame(
  barcode = colnames(scATAC),
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

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

##############################################################################
# Select a specific region overlapping a BAC. Retrieve ATAC and footprint data #
##############################################################################

# Calculate coverage of each overlapping region
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
ovInd <- 253#213#321#153
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
regionStart <- start(ovRegion)
regionEnd <- end(ovRegion)
regionRange <- regionStart:regionEnd
relativeRange <- relativePos:(relativePos + width(ovRegion) - 1)

# For the selected region region get observed Tn5 insertion counts at each position in our BAC data
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

# Get the list of biases generated using each method
biasTracks <- list("predicted" = predBias,
                   "uncorrected" = rep(1, width),
                   "observed" = obsBias)

# Get position-by-pseudobulk matrix of ATAC signal for the current region
ret <- getCountData(project, ovRegionInd)
countData <- ret[["countData"]]
chunkRegionInd <- ret[["regionInd"]]
groupIDs <- 1
chromatinInsertionMatrix <- getRegionATAC(countData, chunkRegionInd, groupIDs, width)

# Call footprints on bone marrow data with Tn5 bias generated using different methods
chromatinFootprints <- pbmcapply::pbmclapply(
  biasTracks,
  function(biasTrack){
    regionFootprintScoring(
      regionATAC = chromatinInsertionMatrix,
      Tn5Bias = biasTrack,
      dispersionModel = dispModel(project, "20"),
      footprintRadius = 20, 
      flankRadius = 20,
      returnCellTypeScores = T,
      scoreType = "pval"
    )$cellTypeScores
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
      dispersionModel = dispModel(project, "20"),
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

maxY <- max(chromatinFootprints$predicted)
plotRadius <- as.integer(width / 2)
positions <- (-plotRadius + 1):plotRadius
plotData <- data.frame(
  "Position" = positions,
  "BACCorrected" = pmin(BACFootprints$predicted, maxY),
  "BACUncorrected" = pmin(BACFootprints$uncorrected, maxY)
)

# Plot BAC Tn5 insertion counts
barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = BACInsertionTrack)
p1 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("BAC Tn5\ninsertion") +
  theme_classic() 

# Plot observed Tn5 bias
barData <- data.frame(x1 = positions - 2, x2 = positions + 2,
                      y1 = 0, y2 = obsBias)
p2 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Observed\nTn5 bias") +
  theme_classic()

# Plot predicted Tn5 bias
barData <- data.frame(x1 = positions - 2, x2 = positions + 2,
                      y1 = 0, y2 = predBias)
p3 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Predicted\nTn5 bias") +
  theme_classic()

# Plot BAC footprints computed with bias correction
p4 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "BACCorrected"), 
              fill = "#69ACD5") + xlab("") + ylab("BAC corrected\nfootprints") +
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) +
  theme_classic()

# Plot BAC footprints computed without bias correction
p5 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "BACUncorrected"), 
              fill = "#69ACD5") + xlab("") + ylab("BAC uncorrected\nfootprints") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

library(patchwork)
system("mkdir ../../data/HepG2/plots/footprintTracks/")
pdf(paste0("../../data/HepG2/plots/footprintTracks/BAC_footprints_", ovRegionInd, ".pdf"),
    width=15, height=8)
p1 / p2 / p3 / p4 / p5 
dev.off()

# Plot chromatin Tn5 insertion counts
chromatinInsertion <- rowSums(chromatinInsertionMatrix)
barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = chromatinInsertion)
plotBMTn5 <- ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#006838") +
  scale_y_continuous(expand = c(0, 0)) + xlab("") + ylab("HepG2\nTn5 insertion") +
  theme_classic() 

# Plot observed - expected deviation of Tn5 insertion
expectedRatio <- conv(predBias, 20) / conv(predBias, 40)
observedRatio <- conv(chromatinInsertion, 20) / conv(chromatinInsertion, 40)
deviation <- observedRatio - expectedRatio
plotDeviation <- ggplot() +
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = deviation), 
              alpha = 0.5) + 
  ylab("HepG2\nObs-Exp\ndeviation") +
  xlab("") +
  theme_classic()

# Plot chromatin footprints computed with model-based bias correction
plotPred <- ggplot() +
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, chromatinFootprints$predicted)), 
              alpha = 0.5) + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("HepG2 model\ncorrected footprints") +
  xlab("") +
  theme_classic()

# Plot chromatin footprints computed with observed bias correction
plotObs <- ggplot() +
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, chromatinFootprints$observed)), 
              alpha = 0.5) + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("HepG2\nground truth\ncorrected footprints") +
  xlab("") +
  theme_classic()

# # Plot chromatin footprints computed without bias correction
plotUncorrected <- ggplot() +
  geom_ribbon(aes_string(x = positions, ymin = 0, 
                         ymax = pmin(maxY, chromatinFootprints$uncorrected)), 
              alpha = 0.5) + 
  scale_y_continuous(limits = c(0,maxY), expand = c(0, 0)) + 
  ylab("HepG2\nuncorrected\nfootprints") +
  xlab("") +
  theme_classic()

pdf(paste0("../../data/HepG2/plots/footprintTracks/chromatin_footprints_", ovRegionInd, ".pdf"),
    width=15, height=8)
plotBMTn5 /  plotDeviation / plotPred / plotObs / plotUncorrected 
dev.off()