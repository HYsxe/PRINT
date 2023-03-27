# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getTFBS.R")
source("../../code/getAggregateFootprint.R")
library(hdf5r)

###################
# Load input data #
###################

projectName <- "HepG2" # One of "K562", "GM12878", "HepG2"
ChIPDataset <- "Unibind"# One of "ENCODE", "Unibind"

# Load footprint SummarizedExperiment object
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

# Load footprints
footprintRadii <- c(10, 20, 30, 50, 80, 100)
multiScaleFootprints <- list()
for(footprintRadius in footprintRadii){
  
  print(paste0("Loading data for footprint radius = ", footprintRadius))  
  chunkDir <- paste0("../../data/", projectName, "/chunkedFootprintResults/", footprintRadius, "/")
  chunkFiles <- gtools::mixedsort(list.files(chunkDir))
  scaleFootprints <- pbmcapply::pbmclapply(
    chunkFiles,
    function(chunkFile){
      chunkData <- readRDS(paste0(chunkDir, chunkFile))
      as.data.frame(t(sapply(
        chunkData,
        function(regionData){
          regionData$aggregateScores
        }
      )))
    },
    mc.cores = 16
  )
  scaleFootprints <- data.table::rbindlist(scaleFootprints)
  multiScaleFootprints[[as.character(footprintRadius)]] <- as.matrix(scaleFootprints)
}

# Load PWM data
cisBPMotifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))

# Load TF ChIP ranges
if(ChIPDataset == "ENCODE"){
  TFChIPRanges <- readRDS(paste0("../../data/", projectName, "/ENCODEChIPRanges.rds"))
}else if(ChIPDataset == "Unibind"){
  TFChIPRanges <- readRDS(paste0("../../data/shared/unibind/", projectName, "ChIPRanges.rds"))
}

# Only keep TFs with both motif and ChIP data
keptTFs <- intersect(names(TFChIPRanges), names(cisBPMotifs))
cisBPMotifs <- cisBPMotifs[keptTFs]
TFChIPRanges <- TFChIPRanges[keptTFs]

# Find motif matches across all regions
path <- paste0(dataDir(project), "motifPositionsList.rds")
nCores <- 16
combineTFs <- F
if(file.exists(path)){
  motifMatches <- readRDS(path)
}else{
  regions <- regionRanges(project)
  # Find motif matches across all regions
  print("Getting TF motif matches within CREs")
  motifMatches <- pbmcapply::pbmclapply(
    names(cisBPMotifs),
    function(TF){
      TFMotifPositions <- motifmatchr::matchMotifs(pwms = cisBPMotifs[[TF]], 
                                                   subject = regions, 
                                                   genome = refGenome(project),
                                                   out = "positions")[[1]]
      if(length(TFMotifPositions) > 0){
        TFMotifPositions$TF <- TF
        TFMotifPositions$score <- rank(TFMotifPositions$score) / length(TFMotifPositions$score)
        TFMotifPositions <- mergeRegions(TFMotifPositions)
        TFMotifPositions
      }
    },
    mc.cores = nCores
  )
  names(motifMatches) <- names(cisBPMotifs)
  if(combineTFs){motifMatches <- Reduce(c, motifMatches)}
  saveRDS(motifMatches, path)
}

# Get ATAC tracks for each region
project <- getATACTracks(project)
regionATAC <- ATACTracks(project)

##########################################
# Get multi-scale footprints around TFBS #
##########################################

TFBSData <- getTFBSTrainingData(multiScaleFootprints,
                                motifMatches,
                                TFChIPRanges,
                                regions, 
                                percentBoundThreshold = 0.1)

#########################################
# Get Tn5 insertion profile around TFBS #
#########################################

# Write TFBS training data to a file
h5_path <- paste0("../../data/", projectName, "/TFBSData", ChIPDataset, ".h5")
if(file.exists(h5_path)){
  system(paste0("rm ../../data/", projectName, "/TFBSData", ChIPDataset, ".h5"))
}
h5file <- H5File$new(h5_path, mode="w")
h5file[["motifFootprints"]] <- TFBSData[["motifFootprints"]]
h5file[["metadata"]] <- TFBSData[["metadata"]]
h5file$close_all()

####################################
# Visualize footprints around TFBS #
####################################

h5_path <- paste0("../../data/", projectName, "/TFBSData", ChIPDataset, ".h5")
h5file <- H5File$new(h5_path, mode="r")
motifFootprints <- h5file[["motifFootprints"]]
metadata <- h5file[["metadata"]]
motifFootprints <- motifFootprints[1:motifFootprints$dims[1],]
metadata <- metadata[1:metadata$dims]
h5file$close_all()
footprintRadii <- c(10, 20, 30, 50, 80, 100)
contextRadius <- 100

library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

# Split heatmap columns by kernel size
colGroups <- Reduce(c, lapply(
  footprintRadii,
  function(r){
    paste(rep(r, contextRadius * 2 + 1), "bp")
  }
))

# Label each kernel size
colNames <- Reduce(c, lapply(
  footprintRadii,
  function(r){
    c(rep("", contextRadius),
      paste(" ", r, "bp"),
      rep("", contextRadius))
  }
))

TF <- "ARNTL"
TFFilter <- metadata[, 3] %in% TF
TFFingerprint <- motifFootprints[TFFilter, ]
sampleInd <- c(sample(which(metadata[TFFilter,1] == 1), 500, replace = T),
               sample(which(metadata[TFFilter,1] == 0), 500, replace = T))
TFBinding <- c("Unbound", "Bound")[(metadata[TFFilter,1] + 1)[sampleInd]]
rowOrder <- order(TFFingerprint[sampleInd, 502], decreasing = T)
colors <- colorRamp2(seq(0, quantile(TFFingerprint, 0.95), length.out=9),
                     colors = colorRampPalette(c(rep("white", 2),  
                                                 "#9ECAE1", "#08519C", "#08306B"))(9))

Heatmap(TFFingerprint[sampleInd[rowOrder],],
        col = colors,
        cluster_rows = F,
        show_row_dend = F,
        show_column_dend = F,
        column_split = factor(colGroups, levels = paste(footprintRadii, "bp")),
        row_split = factor(TFBinding[rowOrder], levels = c("Unbound", "Bound")),
        cluster_columns = F,
        cluster_column_slices = F,
        name = paste(TF, "\nfootprint\nscore"),
        border = TRUE,
        column_title = "Kernel sizes",
        column_title_side = "bottom",
        bottom_annotation = HeatmapAnnotation(
          text = anno_text(colNames, rot = 0, location = 0.5, just = "center")
        )
)
