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

# Initialize a footprintingProject object
projectName <- "HepG2"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Load region ranges
regionBed <- read.table(paste0(projectDataDir, "peaks.bed"))
regions <- GRanges(paste0(regionBed$V1, ":", regionBed$V2, "-", regionBed$V3))
regions <- resize(regions, 1000, fix = "center")
saveRDS(regions, paste0(projectDataDir, "regionRanges.rds"))

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
scATAC <- readRDS(paste0(projectDataDir, "scATAC.rds"))
barcodeGroups <- data.frame(
  barcode = colnames(scATAC),
  group = 1L
)
barcodeGrouping(project) <- barcodeGroups
groups(project) <- mixedsort(unique(barcodeGroups$group))
groupCellType(project) <- c(projectName)

# Getting the pseudobulk-by-region-by-position counts tensor from a fragment file
pathToFrags <- paste0(projectDataDir, "all.frags.filt.tsv.gz")
project <- getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T)

# Load background dispersion model
for(kernelSize in 2:100){
  dispModel(project, as.character(kernelSize)) <-
    readRDS(paste0(projectMainDir,"data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Specify which programs to run
args <- commandArgs(trailingOnly=TRUE)
options <- rep(F, 2)
options[as.integer(args[1])] <- T
runTFBS <- options[1]
runFootprinting <- options[2]

# Run footprinting
if(runFootprinting){
  
  # Set footprint radius
  if(is.na(args[2])){
    footprintRadius <- 20
  }else{
    footprintRadius <- as.integer(args[2])
  }
  
  project <- getFootprints(
    project, 
    mode = as.character(footprintRadius),
    nCores = 16, 
    footprintRadius = footprintRadius,
    flankRadius = footprintRadius)
  
  system(paste0("mkdir ../../data/", projectName, "/footprints/"))
  saveRDS(footprints(project, as.character(footprintRadius)), 
          paste0(projectDataDir, "footprints/", footprintRadius, "bp_footprints.rds"))
}

# Run TFBS detection
if(runTFBS){
  
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
  
  getTFBS(project,
          nCores = 12)
  
}
