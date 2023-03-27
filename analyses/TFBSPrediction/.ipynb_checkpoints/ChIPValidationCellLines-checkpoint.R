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

######################################
# Load TFBS prediction and ChIP data #
######################################

# Initialize a footprintingProject object
projectName <- "K562"
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

# Get the TF binding score SummarizedExperiment object
TFBindingSE <- getTFBindingSE(project)
TFBSRangesAll <- rowRanges(TFBindingSE)
TFBSScoresAll <- assay(TFBindingSE)

######################
# Retrieve ChIP data #
######################

dataset <- "K562" # Can be one of "K562", "GM12878", "HepG2", "Partridge2020"

if(dataset == "Partridge2020"){

  # Load ChIP data from this paper: https://www.nature.com/articles/s41586-020-2023-4
  ChIPFolder <- "../../data/HepG2/Partridge2020ChIP/"
  ChIPBeds <- list.files(ChIPFolder)
  TFNames <- sapply(ChIPBeds, function(x){strsplit(x, "_")[[1]][2]})
  names(ChIPBeds) <- TFNames
  
  # Load chain file needed for liftOver
  pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
  ch <- rtracklayer::import.chain(pathToChain)
  
  TFChIPRanges <- pbmcapply::pbmclapply(
    TFNames,
    function(TF){
      # Get genomic ranges of ChIP regions
      ChIPBed <- read.table(paste0(ChIPFolder, ChIPBeds[TF]))
      hg19ChIPRanges <- GRanges(paste0(ChIPBed$V1, ":", ChIPBed$V2, "-", ChIPBed$V3))
      
      # Lift over from hg19 to hg38
      seqlevelsStyle(hg19ChIPRanges) <- "UCSC"  # necessary
      hg38ChIPRanges<- unlist(rtracklayer::liftOver(hg19ChIPRanges, ch))
      
      hg38ChIPRanges
    },
    mc.cores = 16
  )
  names(TFChIPRanges) <- TFNames
  
}else{
  
  TFChIPRanges <- readRDS(paste0("../../data/", dataset, "/ENCODEChIPRanges.rds"))
}

#########################################################
# Evaluate model performance using ChIP as ground truth #
#########################################################

threshold <- 0.75
performance <- pbmcapply::pbmclapply(
  names(TFChIPRanges),
  function(TF){
    TFBSRanges <- TFBSRangesAll[TFBSRangesAll$TF == TF]
    TFBSScores <- TFBSScoresAll[TFBSRangesAll$TF == TF, 1]
    predBoundRanges <- TFBSRanges[TFBSScores > threshold]
    nMotifs <- length(TFBSRanges)
    nTruePositive <- length(subsetByOverlaps(predBoundRanges, TFChIPRanges[[TF]]))
    nPredPositive <- length(predBoundRanges)
    nBound <- length(subsetByOverlaps(TFBSRanges, TFChIPRanges[[TF]]))
    precision <- nTruePositive / nPredPositive
    recovery <- nTruePositive
    motifBoundR <- length(subsetByOverlaps(TFBSRanges, TFChIPRanges[[TF]])) / nMotifs
    predBoundR <- nPredPositive / length(TFBSRanges)
    data.frame("nPred" = nPredPositive,
               "precision" = precision,
               "recovery" = recovery,
               "predBoundR" = predBoundR,
               "motifBoundR" = motifBoundR,
               "TF" = TF)
  },
  mc.cores = 16
)
performance <- data.table::rbindlist(performance)
performance <- as.data.frame(performance)
performance[is.na(performance)] <- 0
performance <- performance[performance$motifBoundR > 0.2, ]
median(performance$precision)
mean(performance$recovery)

# Save results to file
write.table(performance, 
            paste0("../../data/TFBSPrediction/", projectName, "Performance.txt"),
            quote = F, sep = "\t")
