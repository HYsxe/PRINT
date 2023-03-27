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
projectName <- "mESC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "mm10")
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

for(footprintRadius in 2:100){
  dispModel(project, as.character(footprintRadius)) <- readRDS(paste0("../../data/shared/dispModel/dispersionModel", footprintRadius, "bp.rds"))
}

###############################
# Get nucleosome mapping data #
###############################

# Get chemically mapped nucleosome positions (mm9)
nucData <- data.table::fread("../../data/mESC/nucleosomeChemicalMapping/GSM2183909_unique.map_95pc.txt.gz",
                             showProgress = T)
nucData$V1 <- paste0("chr", nucData$V1)

# According to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2183909
# chr20-22 are X, Y and M, respectively
nucData$V1[nucData$V1 == "chr20"] <- "chrX"
nucData$V1[nucData$V1 == "chr21"] <- "chrY"
nucData$V1[nucData$V1 == "chr22"] <- "chrM"
nucRanges <- GRanges(seqnames = nucData$V1,
                     ranges = IRanges(start = nucData$V2,
                                      end = nucData$V2))

# Get ChIP region genomic ranges for each TF
pathToChain <- "../../data/shared/mm9ToMm10.over.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Lift over from hg19 to hg38
seqlevelsStyle(nucRanges) <- "UCSC"  # necessary
nucRanges <- unlist(rtracklayer::liftOver(nucRanges, ch))

# Only keep those in CRE regions
nucRanges <- subsetByOverlaps(nucRanges, regions)

###############################
# Get nucleosome mapping data #
###############################

# Get region-by-position Tn5 insertion count matrix
if(file.exists("../../data/mESC/ATACTracks.rds")){
  ATACTracks(project) <- readRDS("../../data/mESC/ATACTracks.rds")
}else{
  project <- getATACTracks(project)
  saveRDS(ATACTracks(project), "../../data/mESC/ATACTracks.rds")
}

# Get aggregate ATAC profile around selected sites
aggregateFP <- getAggregateFootprint(
  project = project, 
  sites = nucRanges[sample(1:length(nucRanges), 1e4)]
)

# Multi-scale footprinting
footprintRadii <- 2:100
width <- unique(regionWidth(project))
plotRadius <- 500
plotRange <- (width - plotRadius + 1):(width + plotRadius)
nucMultiScale <- pbmcapply::pbmcmapply(
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

footprintColors <- colorRamp2(seq(quantile(nucMultiScale, 0.01),
                                  1,length.out=9),
                              colors = jdb_palette("brewer_blue"))

colnames(nucMultiScale) <- rep("", length(footprintRadii))
colnames(nucMultiScale)[footprintRadii %% 10 == 0] <- footprintRadii[footprintRadii %% 10 == 0] * 2

positions <- 1:dim(nucMultiScale)[1]
rownames(nucMultiScale) <- rep("", dim(nucMultiScale)[1])
rownames(nucMultiScale)[positions %% 100 == 0] <- positions[positions %% 100 == 0] - plotRadius

system(paste0("mkdir -p ../../data/", projectName, "/plots/aggregateMultiScaleFootprint/"))
pdf(paste0("../../data/", projectName, "/plots/aggregateMultiScaleFootprint/nucleosome.pdf"),
    width = 6, height = 5)
Heatmap(t(nucMultiScale[, rev(1:length(footprintRadii))]),
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
dev.off()
