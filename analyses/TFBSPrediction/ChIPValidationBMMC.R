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

#################
# Get ChIP data #
#################

# Load PWM data
cisBPMotifs <- readRDS("../../data/shared/cisBP_human_pwms_2021.rds")

# Load metadata of ChIP datasets
ChIPMetadata <- read.csv("../../data/shared/cistromeDB/human_factor_full_QC.txt", sep = "\t")

# Filter by QC metrics
# Thresholds are specified here http://cistrome.org/db/#/about
QCFilter <- (ChIPMetadata$FRiP >= 0.01) &
  (ChIPMetadata$FastQC >= 25) &
  (ChIPMetadata$UniquelyMappedRatio >= 0.6) &
  (ChIPMetadata$PeaksFoldChangeAbove10 >= 500) &
  (ChIPMetadata$PeaksUnionDHSRatio >= 0.7) &
  (ChIPMetadata$PBC >= 0.8)
QCFilter[is.na(QCFilter)] <- F

# Filter by cell type
cellTypeFilters <- list(
  "Monocyte" = ChIPMetadata$Cell_type %in% c("Monocyte"),
  "Bcell" = ChIPMetadata$Cell_type %in% c("B Lymphocyte"),
  "Erythroid" = ChIPMetadata$Cell_type %in% c("Erythroid Cell", "Erythroid Progenitor Cell",
                                              "Erythroid progenitor")
)

for(cellType in names(cellTypeFilters)){
  
  # Combine QC and cell type filters
  ChIPMetadataFilt <- ChIPMetadata[QCFilter & cellTypeFilters[[cellType]],]
  
  # If there are multiple datasets for the same factor, pick the one with highest FRIP
  ChIPMetadataFilt <- ChIPMetadataFilt %>% group_by(Factor) %>% filter(FRiP == max(FRiP)) 
  
  # Only keep TFs with motif data
  ChIPMetadataFilt <- ChIPMetadataFilt[ChIPMetadataFilt$Factor %in% names(cisBPMotifs),]
  
  ChIPRangesAll <- lapply(
    1:dim(ChIPMetadataFilt)[1],
    function(i){
      
      # Load ChIP data for the TF
      ChIPInfo <- ChIPMetadataFilt[i,]
      ChIPBed <- read.table(paste("../../data/shared/cistromeDB/human_factor/", ChIPInfo$DCid,
                                  "_sort_peaks.narrowPeak.bed", sep = ""), 
                            sep = "\t")
      
      # Conevert bed format to GRanges
      ChIPRanges <- GenomicRanges::GRanges(seqnames = ChIPBed[,1], 
                                           ranges = IRanges::IRanges(start = as.integer(ChIPBed[,2]), 
                                                                     end = as.integer(ChIPBed[,3])))
      ChIPRanges$score <- ChIPBed[, 5]
      
      ChIPRanges
    }
  )
  names(ChIPRangesAll) <- ChIPMetadataFilt$Factor
  
  # Save results
  saveRDS(ChIPRangesAll, paste0("../../data/shared/cistromeDB/", cellType, "_TF_ChIP.rds"))
}

#############################
# Load TFBS prediction data #
#############################

# Initialize a footprintingProject object
projectName <- "BMMCCellType"
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

# Specify cell type and load TF ChIP data
cellType <- "NaiveB"
ChIPFiles <- list(
  "late-Ery" = "../../data/shared/cistromeDB/Erythroid_TF_ChIP.rds",
  "NaiveB" = "../../data/shared/cistromeDB/Bcell_TF_ChIP.rds",
  "CD14mono" = "../../data/shared/cistromeDB/Monocyte_TF_ChIP.rds"
)
TFChIPRanges <- readRDS(ChIPFiles[[cellType]])

# Load cell type annotations for all pseudo-bulks
# In this case each pseudobulk is a different cell type
groupInfo <- read.table("../../data/BMMCCellType/groupInfo.txt", header = T)
cellTypeInd <- which(groupInfo$cellType == cellType)

# Find motif matches
motifPositions <- getMotifPositions(project, cisBPMotifs, combineTFs = F)

#########################################################
# Evaluate model performance using ChIP as ground truth #
#########################################################

threshold <- 0.75
performance <- pbmcapply::pbmclapply(
  names(TFChIPRanges),
  function(TF){
    tileMotifOv <- findOverlaps(TFBSRangesAll, resize(motifPositions[[TF]], 1))
    TFBSRanges <- TFBSRangesAll[tileMotifOv@from]
    TFBSScores <- TFBSScoresAll[tileMotifOv@from, cellTypeInd]
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
performance <- performance[performance$motifBoundR > 0.1, ]
median(performance$precision)
mean(performance$recovery)

# Save results to file
write.table(performance, 
            paste0("../../data/TFBSPrediction/BMMC", cellType, "Performance.txt"),
            quote = F, sep = "\t")

#########################################
# Visualize results with genomic tracks #
#########################################

pltRange <- GRanges("chr1:25000000-29000000")
pltTFs <- c("ATF2", "BHLHE40", "YY1", "ETS1")
colors <- as.character(MetBrewer::met.brewer("Monet", length(pltTFs), type = "discrete"))
tracks <- list()
for(TFInd in 1:length(pltTFs)){
  
  TF <- pltTFs[TFInd]
  
  # Find predicted binding sites of the TF
  tileMotifOv <- findOverlaps(TFBSRangesAll, resize(motifPositions[[TF]], 1))
  TFBSRanges <- TFBSRangesAll[tileMotifOv@from]
  TFBSScores <- TFBSScoresAll[tileMotifOv@from, cellTypeInd]
  predBoundRanges <- TFBSRanges[TFBSScores > 0.75]
  predBoundRanges <- subsetByOverlaps(predBoundRanges, pltRange)
  
  # Visualize predicted TF binding sites
  tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(predBoundRanges, pltRange, 
                                                        name = paste(TF, "\nPredicted"), stacking = "dense",
                                                        fill = colors[TFInd])
  
  # Only keep ChIP peaks with a matched motif
  ChIPRanges <- subsetByOverlaps(TFChIPRanges[[TF]], TFBSRanges)
  
  # Visualize ChIP peaks
  tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(subsetByOverlaps(ChIPRanges, pltRange), 
                                                        name = paste(TF, "\nChIP"), stacking = "dense",
                                                        fill = colors[TFInd])
  
}

gtrack <- Gviz::GenomeAxisTrack()
tracks <- c(gtrack, tracks)
sizes <- rep(0.5, length(tracks))
system(paste0("mkdir -p ../../data/", projectName, "/plots/ChIPValidation/"))
pdf(paste0("../../data/", projectName, "/plots/ChIPValidation/Tracks.pdf"),
    width = 6, height = 5)
Gviz::plotTracks(tracks, 
                 sizes = sizes,
                 from = start(pltRange), 
                 to = end(pltRange), 
                 chromosome = seqnames(pltRange),
                 col.axis = "black", 
                 fontcolor.legend = "black",
                 lwd=.3,
                 title.width = 1,
                 fontcolor="black",
                 background.title="transparent",
                 col.title="black",
                 col.sampleNames = "black",
                 fontsize=10,
                 cex = 0.8)
dev.off()
