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
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")
source("../../code/footprintTracking.R")

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
pathToFragGrouping <- paste0(projectDataDir, "barcodeGrouping.txt")
barcodeGroups <- read.table(pathToFragGrouping, header = T)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- 
    readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

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

# Get region-by-pseudobulk ATAC matrix 
pseudobulkATACPath <- paste0(projectDataDir, "pseudobulkATAC.rds")
if(!file.exists(pseudobulkATACPath)){
  groupATAC(project) <- getGroupATAC(project)
  groupATAC(project) <- centerCounts(groupATAC(project))
  saveRDS(groupATAC(project), pseudobulkATACPath)
}else{
  groupATAC(project) <- readRDS(pseudobulkATACPath)
}

# Get gene-by-pseudobulk RNA matrix, followed by LogNormalization
pseudobulkRNAPath <- paste0(projectDataDir, "pseudobulkRNA.rds")
if(!file.exists(pseudobulkRNAPath)){
  scRNA <- readRDS(paste0(projectDataDir, "scRNA.rds"))
  RNAMatrix <- scRNA@assays$RNA@data
  pseudobulkRNA <- getGroupRNA(RNAMatrix, barcodeGroups)
  pseudobulkRNA <- t(t(pseudobulkRNA) * 1e6 / colSums(pseudobulkRNA))
  pseudobulkRNA <- log2(pseudobulkRNA + 1)
  groupRNA(project) <- pseudobulkRNA
  saveRDS(pseudobulkRNA, pseudobulkRNAPath)
}else{
  groupRNA(project) <- readRDS(pseudobulkRNAPath)
}

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

# Track nucleosomes across pseudotime
LineageInd <- 2
lineageProbs <- groupInfo[, c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")]
lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                            function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})
lineageGroups <- which(lineageIdentities %in% LineageInd)

nucleosomeTrackingResults <- readRDS("../../data/BMMC/nucleosomeTracks.rds")

###############################################################
# Plot footprint dynamics along a lineage for a specific region #
###############################################################

regionInd <- 158631

plotPseudotimeHeatmap(
  project,
  regionInd = regionInd,
  feature = "TFBS",
  lineageGroups = lineageGroups,
  cellTypeColors = cellTypeColors
)

plotPseudotimeHeatmap(
  project,
  regionInd = regionInd,
  feature = "footprint",
  lineageGroups = lineageGroups,
  cellTypeColors = cellTypeColors,
  footprintRadius = 50
)

######################################################
# Characterize the dynamics of each nucleosome track #
######################################################

# Get pseudotime windows
pseudoTimeWindowSize <- 10
pseudoTimeWindows <- lapply(1:(length(lineageGroups) - pseudoTimeWindowSize + 1),
                            function(i){i:(i + pseudoTimeWindowSize - 1)})

# Get CRE accessibility data for each pseudotime window
lineageGroupATAC <- groupATAC(project)[,lineageGroups]
pseudoTimeWindowATAC <- sapply(
  pseudoTimeWindows,
  function(pseudoTimeWindow){
    rowSums(lineageGroupATAC[,pseudoTimeWindow])
  }
)

# Characterize the dynamics of each track
tracksInfo <- pbmcapply::pbmclapply(
  1:length(nucleosomeTrackingResults),
  function(regionInd){
    regionTracks <- nucleosomeTrackingResults[[regionInd]]
    
    if(length(regionTracks) == 0){
      return(NULL)
    }
    
    regionTrackInfo <-
      lapply(
      1:length(regionTracks),
      function(trackInd){
        track <- as.data.frame(regionTracks[[trackInd]])
        colnames(track) <- c("Position", "Intensity")
        track$Intensity <- pmin(track$Intensity, 300) # Prevent Inf values
        track$pseudoTime <- 1:dim(track)[1]
        filter <- !is.na(track$Position)
        occupancy <- sum(filter) / dim(track)[1]
        track <- track[filter,]
        ptimeIntensitySlope <- unname(coef(lm(track$Intensity ~ track$pseudoTime))[2])
        ptimeATACCor <- cor(track$pseudoTime, pseudoTimeWindowATAC[regionInd, filter])
        ptimePositionCor <- cor(track$pseudoTime, track$Position)
        ptimeIntensityCor <- cor(track$pseudoTime, track$Intensity)
        footprintATACCor <- cor(track$Intensity, pseudoTimeWindowATAC[regionInd, filter])
        maxShift <- max(track$Position) - min(track$Position)
        minIntensity <- min(track$Intensity)
        varIntensity <- var(track$Intensity)
        positionSlope <- unname(coef(lm(track$Position ~ track$pseudoTime))[2])
        if(dim(track)[1] > 0){
          trackStart <- track$Position[1]
          trackEnd <- track$Position[dim(track)[1]]
        }else{
          trackStart <- NA
          trackEnd <- NA
        }
        data.frame("occupancy" = occupancy, 
                   "ptimeIntensitySlope" = ptimeIntensitySlope,
                   "ptimeATACCor" = ptimeATACCor,
                   "footprintATACCor" = footprintATACCor,
                   "ptimePositionCor" = ptimePositionCor,
                   "ptimeIntensityCor" = ptimeIntensityCor,
                   "positionSlope" = positionSlope,
                   "maxShift" = maxShift,
                   "minIntensity" = minIntensity,
                   "varIntensity" = varIntensity,
                   "trackStart" = trackStart,
                   "trackEnd" = trackEnd,
                   "trackInd" = trackInd,
                   "regionInd" = regionInd)
      }
    )
    
    regionTrackInfo <- data.table::rbindlist(regionTrackInfo)
    
    regionTrackInfo
  },
  mc.cores = 16
)
tracksInfo <- data.table::rbindlist(tracksInfo)

# If one pseudotime window contains pseudobulks of multiple cell types, choose the cell type
# with the highest frequency as the cell type label for this window
cellTypes <- sapply(
  pseudoTimeWindows,
  function(pseudoTimeWindow){
    names(sort(table(groupCellType(project)[lineageGroups][pseudoTimeWindow]), decreasing = T)[1])
  }
)

##################################################
# Loop through individual examples and plot them #
##################################################

mode <- "binding_eviction"
if(mode == "sliding"){
  tracksInfoFilt <- tracksInfo[which((abs(tracksInfo$positionSlope) > 0.5) & 
                                       (tracksInfo$minIntensity > 1)), ]
}else if(mode == "binding_eviction"){
  tracksInfoFilt <- tracksInfo[which((abs(tracksInfo$ptimeIntensitySlope) > 0.1) & 
                                       (abs(tracksInfo$ptimeATACCor) < 0.2) & # This is to prevent confounding effect of CRE opening/closing
                                       (abs(tracksInfo$ptimeIntensityCor) > 0.5)), ]
}

# Visualize nucleosome dynamics
ptime <- 1:length(pseudoTimeWindows)
smoothRadius <- 10
plotDf <- sapply(
  1:dim(tracksInfoFilt)[1],
  function(i){
    regionInd <- tracksInfoFilt$regionInd[i]
    trackInd <- tracksInfoFilt$trackInd[i]
    
    if(mode == "sliding"){
      position <- nucleosomeTrackingResults[[regionInd]][[trackInd]][,1]
      displacement <- position - position[1]
      values <- displacement
    }
    
    if(mode == "binding_eviction"){
      intensity <- nucleosomeTrackingResults[[regionInd]][[trackInd]][,2]
      intensityChange <- intensity - intensity[1]
      values <- intensityChange
    }
    
    values <- sapply(
      1:length(pseudoTimeWindows),
      function(j){
        mean(values[max(1, j - smoothRadius) : min(j + smoothRadius, length(values))])
      }
    )
    values
  }
)

# Generate color map for heatmap colors
dataColors <- colorRamp2(seq(quantile(plotDf, 0.01), quantile(plotDf, 0.99),length.out=9),
                         colors = jdb_palette("solar_extra"))

# Generate side color bar indicating pseudotime
ptimeColors <- colorRamp2(seq(0, dim(plotDf)[1],length.out=9),
                          colors = jdb_palette("solar_rojos"))
topAnno <- HeatmapAnnotation("Pseudotime" = 1:dim(plotDf)[1], 
                          border=TRUE,
                          gap = unit(2, "points"),
                          col=list(Pseudotime = ptimeColors),
                          show_annotation_name=FALSE)

# Visualize curves of nucleosome dynamics
system(paste0("mkdir -p ../../data/BMMC/plots/nucleosomeTracks/", mode))
if(mode == "sliding"){
  pdf("../../data/BMMC/plots/nucleosomeTracks/slidingHeatmap.pdf", width = 10, height = 8)
  print(Heatmap(t(plotDf), 
                row_title = "Nucleosome tracks",
                row_title_side = "left",
                col = dataColors,
                cluster_rows = T, 
                cluster_columns = F,
                name = "Displacement",
                top_annotation = topAnno))
  dev.off()
}else if(mode == "binding_eviction"){
  pdf("../../data/BMMC/plots/nucleosomeTracks/bindingEvictionHeatmap.pdf", width = 10, height = 8)
  print(Heatmap(t(plotDf), 
                row_title = "Nucleosome tracks",
                row_title_side = "left",
                col = dataColors,
                cluster_rows = T, 
                cluster_columns = F,
                name = "Change in\nfootprint score",
                top_annotation = topAnno))
  dev.off()
}

# Loop through individual examples and plot them
for(regionInd in tracksInfoFilt$regionInd){
  
  print(paste("Plotting region", regionInd))
  
  system(paste0("mkdir -p ../../data/BMMC/plots/nucleosomeTracks/", mode))
  pdf(paste0("../../data/BMMC/plots/nucleosomeTracks/", mode, "/Region_", regionInd, "_TFBS.pdf"))
  print(plotPseudotimeHeatmap(
    project,
    regionInd = regionInd,
    feature = "TFBS",
    lineageGroups = lineageGroups,
    cellTypeColors = cellTypeColors
  ))
  dev.off()
  
  pdf(paste0("../../data/BMMC/plots/nucleosomeTracks/", mode, "/Region_", regionInd, "_50bp_footprints.pdf"))
  print(plotPseudotimeHeatmap(
    project,
    regionInd = regionInd,
    feature = "footprint",
    lineageGroups = lineageGroups,
    cellTypeColors = cellTypeColors,
    heatmapPalette = colorRampPalette(c("white",  "#9ECAE1", "#08519C", "#08306B"))(9),
    footprintRadius = 50
  ))
  dev.off()
  
}
