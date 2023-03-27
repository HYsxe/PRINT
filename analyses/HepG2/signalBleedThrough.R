# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getTFBS.R")
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

###################
# Load input data #
###################

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

# Load TFBS prediction model
h5Path <- "../../data/TFBSPrediction/TFBS_model_cluster_I.h5"
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

fpScoring <- function(Tn5Insertion,
                      Tn5Bias,
                      dispersionModel,
                      footprintRadius){
  footprintPvals <- footprintScoring(
    Tn5Insertion = Tn5Insertion,
    Tn5Bias = Tn5Bias,
    dispersionModel = dispersionModel,
    footprintRadius = footprintRadius,
    flankRadius = footprintRadius
  )
  footprintScores <- -log10(footprintPvals)
  footprintScores <- caTools::runmax(footprintScores, 5)
  footprintScores <- conv(footprintScores, 5) / 10
  footprintScores[1:footprintRadius] <- 0
  footprintScores[(length(Tn5Bias) - footprintRadius):length(Tn5Bias)] <- 0
  footprintScores
}

###############################
# Plot multi-scale footprints #
###############################

regionInd <- 103294#77781#29253#77781#132408#13880#11863#103294#119538#69221 #55812
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
    footprintScores <- fpScoring(Tn5Insertion = rowSums(regionATAC),
                                 Tn5Bias = Tn5Bias,
                                 dispersionModel = dispModel(project, as.character(footprintRadius)),
                                 footprintRadius = footprintRadius)
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
pdf(paste0("../../data/HepG2/plots/signalBleedThrough/region_", regionInd, "_multiScale.pdf"),
    height = 6, width = 7)
p <- Heatmap(t(scaleScanning[plotPositions, rev(1:length(footprintRadii))]),
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
p
dev.off()
p

###############################
# Plot multi-scale footprints #
###############################

TFBSModel <- TFBindingModel(project)
scales <- TFBSModel$scales # Scales are footprint radii we use for multi-scale footprinting
dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
names(dispModels) <- scales
region <- regionRanges(project)[regionInd]

# Calculate footprint scores for each region of the current batch
regionTFBS <- getRegionTFBS(regionATAC = regionATAC,
                            Tn5Bias = regionBias(project)[regionInd,],
                            region = region,
                            dispModels = dispModels,
                            TFBSModel = TFBindingModel(project),
                            contextRadius = 100,
                            tileSize = 1)

# Compute 20 bp radius footprint scores
footprintRadius <- 20
footprintScores <- fpScoring(Tn5Insertion = rowSums(regionATAC),
                             Tn5Bias = Tn5Bias,
                             dispersionModel = dispModel(project, as.character(footprintRadius)),
                             footprintRadius = footprintRadius)

# Visualize results
contextRadius <- 100
width <- length(Tn5Bias)
pltRange <- (contextRadius + 1):(width - contextRadius)
positions <- - (width / 2 - contextRadius - 1) : (width / 2 - contextRadius)
plotData <- data.frame(
  Position = positions,
  baseline = rep(0, length(pltRange)),
  barLeft = positions - 2,
  barRight = positions + 2,
  Bias = Tn5Bias[pltRange],
  Insertion = rowSums(regionATAC)[pltRange],
  footprintScores = footprintScores[pltRange],
  TFBSScores = pmax(regionTFBS$TFBSScores, 0.05)
)

p1 <- ggplot(plotData) + 
  geom_rect(aes(xmin = barLeft, xmax = barRight, ymin = baseline, ymax = Bias), 
            fill = "#E21F26") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Predicted\nTn5 bias") +
  theme_classic()

p2 <- ggplot(plotData) + 
  geom_rect(aes(xmin = barLeft, xmax = barRight, ymin = baseline, ymax = Insertion), 
            fill = "#7F3F98") +
  scale_y_continuous(expand = c(0, 0)) + ylab("Observed\nTn5 insertion") +
  theme_classic()

p3 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "positions", ymin = 0, ymax = "footprintScores"), 
              fill = "#69ACD5") + xlab("") + ylab("Footprint\n scores") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

p4 <- ggplot(plotData) + 
  geom_ribbon(aes_string(x = "Position", 
                         ymin = 0.05,
                         ymax = "TFBSScores"), alpha = 0.5) + 
  xlab("") + 
  ylab("TF binding\n score") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()

library(patchwork)
pdf(paste0("../../data/HepG2/plots/signalBleedThrough/region_", regionInd, "_tracks.pdf"),
    height = 6, width = 7)
p <- p1 / p2 / p3 / p4
p
dev.off()
p
