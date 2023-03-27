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

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

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
regionRanges(project) <- regions

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getRegionBias(project, nCores = 24)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each pseudo-bulk
scATAC <- readRDS(paste0(projectDataDir, "scATAC.rds"))
barcodeGroups <- data.frame(
  barcode = colnames(scATAC),
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

for(scaleSize in 2:100){
  dispModel(project, as.character(scaleSize)) <- readRDS(paste0("../../data/shared/dispModel/dispersionModel", scaleSize, "bp.rds"))
}

# Load PWM data
motifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))

################################
# Plot multi-scale footprints #
################################

regionInd <- 69221#55812
Tn5Bias <- regionBias(project)[regionInd, ]
width <- regionWidth(project)

# Get the Tn5 insertion count tensor for the current region
ret <- getCountData(project, regionInd)
countData <- ret[["countData"]]
adjustedRegionInd <- ret[["regionInd"]]

# Get position-by-pseudobulk ATAC insertion matrix
regionATAC <- getRegionATAC(countData, adjustedRegionInd, 1, width[regionInd])

# Calculate footprint scores for each position in the region using binomial testing 
footprintRadii <- 2:100
scaleScanning <- pbmcapply::pbmcmapply(
  function(footprintRadius){
    footprintPvals <- footprintScoring(
      Tn5Insertion = rowSums(regionATAC),
      Tn5Bias = Tn5Bias,
      dispersionModel = dispModel(project, as.character(footprintRadius)),
      footprintRadius = footprintRadius,
      flankRadius = footprintRadius
    )
    footprintScores <- -log10(footprintPvals)
    footprintScores <- caTools::runmax(footprintScores, 5)
    footprintScores <- conv(footprintScores, 5) / 10
    footprintScores[1:footprintRadius] <- 0
    footprintScores[(width[regionInd] - footprintRadius):width[regionInd]] <- 0
    footprintScores
  },
  footprintRadii,
  mc.cores = 16
)

footprintColors <- colorRamp2(seq(0.5, min(quantile(scaleScanning, 0.99), 2),length.out=9),
                              colors = jdb_palette("brewer_blue"))

colnames(scaleScanning) <- rep("", length(footprintRadii))
colnames(scaleScanning)[footprintRadii %% 10 == 0] <- footprintRadii[footprintRadii %% 10 == 0] * 2

positions <- 1:width[regionInd]
rownames(scaleScanning) <- rep("", width[regionInd])
rownames(scaleScanning)[positions %% 100 == 0] <- positions[positions %% 100 == 0] - width[regionInd] / 2

plotPositions <- max(footprintRadii):(width[regionInd] - max(footprintRadii))
Heatmap(t(scaleScanning[plotPositions, rev(1:length(footprintRadii))]),
        use_raster = TRUE,
        col = footprintColors,
        cluster_rows = F,
        show_row_dend = F,
        cluster_columns = F,
        name = "Footprint\nScore",
        border = TRUE,
        column_title = "Position (bp)",
        column_title_side = "bottom",
        column_names_rot = 0)

#################
# Get ChIP data #
#################

# Get the list of TFs with ChIP data
ChIPFolder <- paste0("../../data/", projectName, "/ENCODEChIP/")
ChIPMetadata <- read.csv(paste0("../../data/",  projectName, "/metadata.tsv"), sep = "\t")
ChIPMetadata <- ChIPMetadata[ChIPMetadata$File.format == "bed narrowPeak",]
ChIPMetadata <- ChIPMetadata[ChIPMetadata$File.assembly == "GRCh38",]
ChIPTFs <- unname(sapply(
  ChIPMetadata$Experiment.target,
  function(target){
    strsplit(target, "-")[[1]][1]
  }
))

# Retrieve the genomic ranges of the TFs
TFChIPRanges <- pbmcapply::pbmclapply(
  sort(unique(ChIPTFs)),
  function(TF){
    ChIPRangeList <- lapply(
      which(ChIPTFs %in% TF),
      function(entry){
        ChIPBed <- read.csv(paste0(ChIPFolder, ChIPMetadata$File.accession[entry], ".bed.gz"), 
                            sep = "\t", header = F)
        ChIPRanges <- GRanges(
          seqnames = ChIPBed[,1],
          ranges = IRanges(start = ChIPBed[, 2], end = ChIPBed[,3])
        )
        ChIPRanges$signal <- ChIPBed[, 9]
        ChIPRanges
      }
    )
    # For TFs with more than one ChIP files, take the intersection of regions in all files
    ChIPRanges <- Reduce(subsetByOverlaps, ChIPRangeList)
    ChIPRanges
  },
  mc.cores = 16
)
names(TFChIPRanges) <- sort(unique(ChIPTFs))
saveRDS(TFChIPRanges, paste0("../../data/", projectName, "/ENCODEChIPRanges.rds"))

# Also convert the integrated list to bed format
TFChIPBed <- lapply(
  sort(unique(ChIPTFs)),
  function(TF){
    data.frame(
      chr = as.character(seqnames(TFChIPRanges[[TF]])),
      start = start(TFChIPRanges[[TF]]),
      end = end(TFChIPRanges[[TF]]),
      TF = TF
    )
  }
)
TFChIPBed <- as.data.frame(data.table::rbindlist(TFChIPBed))
write.table(TFChIPBed, paste0("../../data/", projectName, "/ENCODEChIPBed.tsv"), 
            sep = "\t", quote = F, row.names = F, col.names = F)

#####################################################
# Muti-scale footprinting of aggregate motif sites #
#####################################################

# Load TF ChIP ranges
TFChIPRanges <- readRDS(paste0("../../data/", projectName, "/ENCODEChIPRanges.rds"))

# Get region-by-position Tn5 insertion count matrix
project <- getATACTracks(project)

selectedTFs <- sort(intersect(names(TFChIPRanges), names(motifs)))
footprintRadii <- 2:100
width <- unique(regionWidth(project))
plotRadius <- 500
plotRange <- (width - plotRadius + 1):(width + plotRadius)
aggregateMultiScale <- list()
for(TF in selectedTFs){
  
  print(paste("Calculating aggregate multi-scale footprint for", TF))
  
  # Find matches of the PWM within all region regions
  motifPositions <- motifmatchr::matchMotifs(pwms = motifs[[TF]], 
                                             subject = regions, 
                                             genome = "hg38",
                                             p.cutoff = 1e-5,
                                             out = "positions")[[1]]
  moftiChIPOv <- findOverlaps(motifPositions, TFChIPRanges[[TF]])
  
  # Calculate percentage of bound motifs
  boundMotifs <- motifPositions[moftiChIPOv@from]
  percentBound <- length(boundMotifs) / length(motifPositions)
  
  # Select motifs with high ChIP signal as bound
  boundMotifs$signal <- TFChIPRanges[[TF]]$signal[moftiChIPOv@to]
  boundMotifs <- boundMotifs[boundMotifs$signal > 2]
  
  # Filter out TFs where only a small portion of motifs are bound
  percentBound <- length(boundMotifs) / length(motifPositions)
  if(percentBound < 0.2){ next }
  if(length(boundMotifs) == 0){ next }
  
  # Get aggregate Tn5 insertion profile centered around motif sites for a specific TF
  aggregateFP <- getAggregateFootprint(
    project, 
    boundMotifs
  )
  
  # For each scale size, calculate footprint scores for the aggregate profile
  aggregateMultiScale[[TF]] <- pbmcapply::pbmcmapply(
    function(footprintRadius){
      footprintPvals <- footprintScoring(
        Tn5Insertion = aggregateFP$uncorrectedATAC[plotRange],
        Tn5Bias = aggregateFP$Tn5Bias[plotRange],
        dispersionModel = dispModel(project, as.character(footprintRadius)),
        footprintRadius = footprintRadius,
        flankRadius = footprintRadius
      )
      # Remove regions affected by edge effect
      # For edges, the sum of bias and counts on the left / right flank might
      # be much lower than what the model has seen, since we're adding many zeros.
      # As a result, the model prediction of ratioSD will not be accurate
      footprintPvals[1 : (2 * footprintRadius)] <- 1
      footprintPvals[(2 * plotRadius - 2 * footprintRadius + 1):(2 * plotRadius)] <- 1
      footprintScores <- -log10(footprintPvals)
      footprintScores <- caTools::runmax(footprintScores, 5)
      footprintScores <- conv(footprintScores, 5) / 10
    },
    footprintRadii,
    mc.cores = 16
  )
  
}

saveRDS(aggregateMultiScale, paste0("../../data/", projectName, "/TFMultiScale.rds"))

TF <- "NFIA"
TFMultiScale <- aggregateMultiScale[[TF]]
footprintColors <- colorRamp2(seq(quantile(TFMultiScale, 0.01),1,length.out=9),
                              colors = jdb_palette("brewer_blue"))

colnames(TFMultiScale) <- rep("", length(footprintRadii))
colnames(TFMultiScale)[footprintRadii %% 10 == 0] <- footprintRadii[footprintRadii %% 10 == 0]

positions <- 1:dim(TFMultiScale)[1]
rownames(TFMultiScale) <- rep("", dim(TFMultiScale)[1])
rownames(TFMultiScale)[positions %% 100 == 0] <- positions[positions %% 100 == 0] - plotRadius

Heatmap(t(TFMultiScale[, rev(1:length(footprintRadii))]),
        use_raster = TRUE,
        col = footprintColors,
        cluster_rows = F,
        show_row_dend = F,
        cluster_columns = F,
        name = "Footprint\nScore",
        border = TRUE,
        column_title = "Position (bp)",
        column_title_side = "bottom",
        column_names_rot = 0)