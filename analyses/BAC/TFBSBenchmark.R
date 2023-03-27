# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getGroupData.R")
source("../../code/getTFBS.R")

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "BAC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regions <- readRDS(paste0(projectDataDir, "tileRanges.rds"))

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

# Load barcodes for each pseudo-bulk
barcodeGroups <- data.frame(barcode = paste("rep", 1:5, sep = ""),
                            group = "BAC")
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
tileCounts <- readRDS("../../data/BAC/tileCounts.rds")
tileCounts <- tileCounts$all
tileCounts <- lapply(
  tileCounts,
  function(x){
    x <- as.data.frame(x)
    if(dim(x)[1] > 0){
      # Now we merge the 5 replicates into 1 group
      x$group <- "BAC"
    }
    x  
  }
)
countTensor(project) <- tileCounts
groupCellType(project) <- "BAC"

for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <- readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Load PWM data
motifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))

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

# Load chain file required for lift-over
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Find motif matches across all peaks
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")
motifMatches <- getMotifPositions(project, cisBPMotifs, combineTFs = F)

#######################
# Run TFBS prediction #
#######################

getTFBS(project,
        nCores = 16)
TFBindingSE <- getTFBindingSE(project)
TFBSScores <- assay(TFBindingSE)
nFootprints <- sum(TFBSScores > 0.6)
print(nFootprints)
print(paste("False positive rate", nFootprints / length(TFBSScores)))

##################################################################
# Load results from other methods and quantify called footprints #
##################################################################

# Retrieve results from HINT-ATAC
HINTFootprints <- read.table("../../data/BAC/HINT/BAC.bed")

# Retrieve results from TOBIAS
TobiasDir <- "../../data/BAC/Tobias/prediction"
TFs <- unname(sapply(list.files(TobiasDir), function(s){stringr::str_split(s, "_")[[1]][[1]]}))
TFs <- intersect(TFs, names(cisBPMotifs))
TOBIASFootprints <- data.table::rbindlist(pbmcapply::pbmclapply(
  TFs,
  function(TF){
    # Retrieve Tobias TFBS prediction
    TobiasBound <- read.table(paste0("../../data/BAC/Tobias/prediction/", TF, "_motif/beds/", TF, "_motif_BAC_bound.bed"))
    TobiasBound
  },
  mc.cores = 16
))
TOBIASFootprintRanges <- GRanges(seqnames = TOBIASFootprints$V1,
                                 ranges = IRanges(start = TOBIASFootprints$V2, end = TOBIASFootprints$V3))
TOBIASFootprintRanges$score <- TOBIASFootprints$V11
TOBIASFootprints <- mergeRegions(TOBIASFootprintRanges)

# Calculate total length of BAC regions
BACLength <- sum(width(regions)) / 1000

methods <- c("TOBIAS", "HINT-ATAC", "CNN-based")
plotData <- data.frame(
  FPR = c(length(TOBIASFootprints) / BACLength, 
          dim(HINTFootprints)[1] / BACLength, 
          nFootprints / BACLength),
  Method = factor(methods, levels = methods)
)
pdf("../../data/BAC/plots/TFBSBenchmark.pdf", width = 4, height = 4)
ggplot(plotData) +
  geom_bar(aes(x = Method, y = FPR), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("") + ylab("Number of false positive TFBS per 1kb\ndetected on BAC naked DNA") + 
  theme_classic()
dev.off()
