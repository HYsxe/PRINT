# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
library(GenomicRanges)

#######################
# Retrieve BAC ranges #
#######################

# Load genomic ranges of BACs
BACRanges <- readRDS("../../data/BAC/BACRanges.rds")

# Convert BAC genomic ranges into 1kb tiles
tileRanges <- Reduce("c", GenomicRanges::slidingWindows(BACRanges, width = 1000, step = 1000))
tileRanges <- tileRanges[width(tileRanges) == 1000]
tileBACs <- names(tileRanges)
saveRDS(tileRanges, "../../data/BAC/tileRanges.rds")

###################################################
# Get BAC predicted bias and Tn5 insertion track #
##################################################

# Initialize a footprintingProject object
projectName <- "BAC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Set the regionRanges slot
regionRanges(project) <- tileRanges

# Use our neural network to predict Tn5 bias for all positions in all regions
# Remember to make a copy of the Tn5_NN_model.h5 file in projectDataDir!!
if(file.exists(paste0(projectDataDir, "predBias.rds"))){
  regionBias(project) <- readRDS(paste0(projectDataDir, "predBias.rds"))
}else{
  project <- getRegionBias(project, nCores = 16)
  saveRDS(regionBias(project), paste0(projectDataDir, "predBias.rds"))
}

# Load barcodes for each replicate
barcodeGroups <- data.frame(barcode = paste("rep", 1:5, sep = ""),
                            group = 1:5)
groups(project) <- mixedsort(unique(barcodeGroups$group))

# Get position-by-tile-by-replicate ATAC insertion count tensor
# We go through all down-sampling rates and get a count tensor for each of them
if(!file.exists("../../data/BAC/tileCounts.rds")){
  counts <- list()
  pathToFrags <- paste0(projectDataDir, "rawData/all.fragments.tsv.gz")
  system("rm -r ../../data/BAC/chunkedCountTensor")
  counts[["all"]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  for(downSampleRate in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
    print(paste0("Getting count tensor for down-sampling rate = ", downSampleRate))
    system("rm -r ../../data/BAC/chunkedCountTensor")
    pathToFrags <- paste0("../../data/BAC/downSampledFragments/fragmentsDownsample", downSampleRate, ".tsv.gz")
    counts[[as.character(downSampleRate)]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  }
  system("rm -r ../../data/BAC/chunkedCountTensor")
  saveRDS(counts, "../../data/BAC/tileCounts.rds")
}else{
  counts <- readRDS("../../data/BAC/tileCounts.rds")
}

##################################
# Get background dispersion data #
##################################

if(!dir.exists("../../data/BAC/dispModelData")){
  system("mkdir ../../data/BAC/dispModelData")
}

for(footprintRadius in seq(2,100,1)){
  
  ############################### Get naked DNA Tn5 observations ###############################
  
  print(paste0("Generating background dispersion data for footprintRadius = ", footprintRadius))
  flankRadius <- footprintRadius
  
  # Record background observations of Tn5 bias and insertion density in the BAC data
  bgObsData <- NULL
  for(downSampleRate in names(counts)){
    bgObsData <- rbind(bgObsData, data.table::rbindlist(
      pbmcapply::pbmclapply(
        1:length(counts[[downSampleRate]]),
        function(tileInd){
          
          # Get predicted bias
          predBias <- regionBias(project)[tileInd, ]
          
          # Get Tn5 insertion
          tileCountTensor <- counts[[downSampleRate]][[tileInd]]
          tileCountTensor <- tileCountTensor %>% group_by(position) %>% summarize(insertion  = sum(count))
          Tn5Insertion <- rep(0, length(predBias))
          Tn5Insertion[tileCountTensor$position] <- tileCountTensor$insertion
          
          # Get sum of predicted bias in left flanking, center, and right flanking windows
          biasWindowSums <- footprintWindowSum(predBias, 
                                               footprintRadius, 
                                               flankRadius)
          
          # Get sum of insertion counts in left flanking, center, and right flanking windows
          insertionWindowSums <- footprintWindowSum(Tn5Insertion, 
                                                    footprintRadius, 
                                                    flankRadius)
          
          # Combine results into a data.frame
          tileObsData <- data.frame(leftFlankBias = biasWindowSums$leftFlank,
                                    leftFlankInsertion = round(insertionWindowSums$leftFlank),
                                    centerBias = biasWindowSums$center,
                                    centerInsertion = round(insertionWindowSums$center),
                                    rightFlankBias = biasWindowSums$rightFlank,
                                    rightFlankInsertion = round(insertionWindowSums$rightFlank))
          tileObsData$leftTotalInsertion <- tileObsData$leftFlankInsertion + tileObsData$centerInsertion
          tileObsData$rightTotalInsertion <- tileObsData$rightFlankInsertion + tileObsData$centerInsertion
          tileObsData$BACInd <- tileBACs[tileInd]
          
          tileObsData
        },
        mc.cores = 8
      )
    ))
  }
  
  # Filter out regions with zero reads
  bgObsData <- bgObsData[(bgObsData$leftTotalInsertion >= 1) &
                           (bgObsData$rightTotalInsertion >= 1),]
  
  # Get features of background observations for left-side testing
  bgLeftFeatures <- as.data.frame(bgObsData[,c("leftFlankBias", "centerBias", "leftTotalInsertion")])
  bgLeftFeatures <- data.frame(bgLeftFeatures)
  bgLeftFeatures$leftTotalInsertion <- log10(bgLeftFeatures$leftTotalInsertion)
  leftFeatureMean <- colMeans(bgLeftFeatures)
  leftFeatureSD <- apply(bgLeftFeatures, 2, sd)
  bgLeftFeatures <- sweep(bgLeftFeatures, 2, leftFeatureMean)
  bgLeftFeatures <- sweep(bgLeftFeatures, 2, leftFeatureSD, "/")
  
  # Get features of background observations for right-side testing
  bgRightFeatures <- as.data.frame(bgObsData[,c("rightFlankBias", "centerBias", "rightTotalInsertion")])
  bgRightFeatures <- data.frame(bgRightFeatures)
  bgRightFeatures$rightTotalInsertion <- log10(bgRightFeatures$rightTotalInsertion)
  rightFeatureMean <- colMeans(bgRightFeatures)
  rightFeatureSD <- apply(bgRightFeatures, 2, sd)
  bgRightFeatures <- sweep(bgRightFeatures, 2, rightFeatureMean)
  bgRightFeatures <- sweep(bgRightFeatures, 2, rightFeatureSD, "/")
  
  ############################### Match background KNN observations and calculate distribution of center / (center + flank) ratio ###############################
  
  # We sample 1e5 data points
  set.seed(123)
  sampleInds <- sample(1:dim(bgObsData)[1], 1e5)
  sampleObsData <- bgObsData[sampleInds, ]
  sampleLeftFeatures <- bgLeftFeatures[sampleInds, ]
  sampleRightFeatures <- bgRightFeatures[sampleInds, ]
  
  # Match KNN observations in (flank bias, center bias, count sum) 3-dimensional space 
  leftKNN <- FNN::get.knnx(bgLeftFeatures, sampleLeftFeatures, k = 500)$nn.index
  leftRatioParams <- t(sapply(
    1:dim(leftKNN)[1],
    function(i){
      KNNObsData <- bgObsData[leftKNN[i,],]
      KNNRatio <- KNNObsData$centerInsertion / KNNObsData$leftTotalInsertion
      ratioMean <- mean(KNNRatio)
      ratioSD <- sd(KNNRatio)
      c(ratioMean, ratioSD)
    }
  ))
  leftRatioParams <- as.data.frame(leftRatioParams)
  colnames(leftRatioParams) <- c("leftRatioMean", "leftRatioSD")
  
  rightKNN <- FNN::get.knnx(bgRightFeatures, sampleRightFeatures, k = 500)$nn.index
  rightRatioParams <- t(sapply(
    1:dim(rightKNN)[1],
    function(i){
      KNNObsData <- bgObsData[rightKNN[i,],]
      KNNRatio <- KNNObsData$centerInsertion / KNNObsData$rightTotalInsertion
      ratioMean <- mean(KNNRatio)
      ratioSD <- sd(KNNRatio)
      c(ratioMean, ratioSD)
    }
  ))
  rightRatioParams <- as.data.frame(rightRatioParams)
  colnames(rightRatioParams) <- c("rightRatioMean", "rightRatioSD")
  
  # Combine background observations on both sides
  dispModelData <- cbind(sampleObsData,
                         leftRatioParams, 
                         rightRatioParams)
  dispModelData$leftTotalInsertion <- log10(dispModelData$leftTotalInsertion)
  dispModelData$rightTotalInsertion <- log10(dispModelData$rightTotalInsertion)
  
  # Save background observations to file
  write.table(dispModelData, 
              paste0("../../data/BAC/dispModelData/dispModelData", footprintRadius, "bp.txt"),
              quote = F, sep = "\t")
}

#########################
# Load dispersion model #
#########################

library(keras)
for(footprintRadius in 2:100){
  
  # For each specific footprint radius, load dispersion model and training data
  dispModelData <- read.table(paste0("../../data/BAC/dispModelData/dispModelData", footprintRadius, "bp.txt"),sep = "\t")
  dispersionModel <- keras::load_model_hdf5(paste0("../../data/shared/dispModel/dispersionModel", footprintRadius, "bp.h5"))
  modelWeights <- keras::get_weights(dispersionModel)
  
  modelFeatures <- dispModelData[, c("leftFlankBias", "rightFlankBias", "centerBias", "leftTotalInsertion", "rightTotalInsertion")]
  modelTargets <- dispModelData[, c("leftRatioMean", "leftRatioSD", "rightRatioMean", "rightRatioSD")]
  
  # Organize everything we need for background dispersion modeling into a list
  dispersionModel <- list(
    
    modelWeights = modelWeights,
    featureMean = colMeans(modelFeatures),
    featureSD = apply(modelFeatures, 2, sd),
    targetMean = colMeans(modelTargets),
    targetSD = apply(modelTargets, 2, sd)
    
  )
  
  saveRDS(dispersionModel, paste0("../../data/shared/dispModel/dispersionModel", footprintRadius ,"bp.rds"))
  
}
