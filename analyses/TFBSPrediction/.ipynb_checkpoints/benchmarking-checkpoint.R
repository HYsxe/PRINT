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
library(hdf5r)
require(PRROC)
library(GenomicRanges)

###################
# Load input data #
###################

projectName <- "K562"

# Load TF ChIP ranges
TFChIPRanges <- readRDS(paste0("../../data/shared/unibind/", projectName, "ChIPRanges.rds"))

# Load chain file required for lift-over
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Load TF cluster membership
TFClustering <- read.csv("../../data/TFBSPrediction/clusterLabelsAllTFs.txt", sep = "\t")

# Load multi-scale TFBS predictions
multiScalePred <- read.table(paste0("../../data/TFBSPrediction/", projectName, "_pred_data.tsv"))
multiScaleRanges <- GRanges(multiScalePred$range)
strand(multiScaleRanges) <- "*"
multiScaleRanges$score <- multiScalePred$predScore
multiScaleRanges$TF <- multiScalePred$TF

# Get the list of TFs
TobiasFolder <- paste0(paste0("../../../SHARE_footprinting_v4/data/", projectName, "/Tobias/prediction"))
TobiasTFs <- unname(sapply(list.files(TobiasFolder), function(x){strsplit(x, "_")[[1]][1]}))
TFs <- Reduce(intersect, list(TobiasTFs, unique(multiScaleRanges$TF), names(TFChIPRanges)))
cluster1TFs <- intersect(TFClustering$TF[TFClustering$cluster == 1], TFs)

###############################
# Evalutate model performance #
###############################

HINTBed <- read.table(paste0("../../../SHARE_footprinting_v4/data/", projectName, "/HINT/footprints.bed"))
HINTRanges <- GRanges(seqnames = HINTBed$V1,
                      ranges = IRanges(start = HINTBed$V2, end = HINTBed$V3))
HINTRanges$score <- HINTBed$V5
seqlevelsStyle(HINTRanges) <- "UCSC"  # necessary
HINTRanges <- unlist(rtracklayer::liftOver(HINTRanges, ch))

benchmarkResults <- t(pbmcapply::pbmcmapply(
  function(TF){
    
    # Load Tobias prediction scores
    TobiasBed <- read.table(paste0("../../../SHARE_footprinting_v4/data/", projectName, "/Tobias/prediction/", TF, 
                                   "_motif/beds/", TF, "_motif_all.bed"))
    TFTobiasRanges <- GRanges(seqnames = TobiasBed$V1,
                              ranges = IRanges(start = TobiasBed$V2, end = TobiasBed$V3))
    TFTobiasRanges$score <- TobiasBed$V10
    
    # Lift over from hg19 to hg38
    seqlevelsStyle(TFTobiasRanges) <- "UCSC"  # necessary
    TFTobiasRanges <- unlist(rtracklayer::liftOver(TFTobiasRanges, ch))
    
    # Retrieve multi-scale prediction for the current TF
    TFMultiScaleRanges <- multiScaleRanges[multiScaleRanges$TF == TF]
    
    # Only keep the intersection of multi-scale motif sites and Tobais motif sites
    motifRanges <- subsetByOverlaps(subsetByOverlaps(TFMultiScaleRanges, TFTobiasRanges), HINTRanges)
    TFMultiScaleRanges <- subsetByOverlaps(TFMultiScaleRanges, motifRanges)
    TFTobiasRanges <- subsetByOverlaps(TFTobiasRanges, motifRanges)
    TFHINTRanges <- subsetByOverlaps(HINTRanges, motifRanges)
    
    # Sort motif sites by TF binding score
    TFMultiScaleRanges <- TFMultiScaleRanges[order(TFMultiScaleRanges$score, decreasing = T)]
    TFTobiasRanges <- TFTobiasRanges[order(TFTobiasRanges$score, decreasing = T)]
    TFHINTRanges <- TFHINTRanges[order(TFHINTRanges$score, decreasing = T)]
    
    # Keep the same number of top sites for both methods for fair comparison
    nPred <- as.integer(length(TFMultiScaleRanges) * 0.1)
    TFMultiScaleRanges <- TFMultiScaleRanges[1:nPred]
    TFTobiasRanges <- TFTobiasRanges[1:nPred]
    TFHINTRanges <- TFHINTRanges[1:nPred]
    
    # Calculate precision of the two methods
    TobiasPrecision <- length(subsetByOverlaps(TFTobiasRanges, TFChIPRanges[[TF]]))/length(TFTobiasRanges)
    multiScalePrecision <- length(subsetByOverlaps(TFMultiScaleRanges, TFChIPRanges[[TF]]))/length(TFMultiScaleRanges)
    HINTPrecision <- length(subsetByOverlaps(TFHINTRanges, TFChIPRanges[[TF]]))/length(TFHINTRanges)
    motifPrecision <-  length(subsetByOverlaps(motifRanges, TFChIPRanges[[TF]]))/length(motifRanges)
    c(multiScalePrecision, HINTPrecision, TobiasPrecision, motifPrecision, nPred) 
    
  },
  TFs,
  mc.cores = 16
))
colnames(benchmarkResults) <- c("MultiScale", "HINT", "Tobias", "Motif", "nPred")
colMeans(benchmarkResults)
colMeans(benchmarkResults[cluster1TFs,])
colMedians(benchmarkResults)
colMedians(benchmarkResults[cluster1TFs,])

#####################
# Visualize results #
#####################

methods = c("MultiScale", "HINT", "Tobias", "Motif")
plotData <- data.frame(
  precision = colMedians(as.matrix(benchmarkResults[, methods])),
  method = factor(methods, levels = methods))
pdf("../../data/TFBSPrediction/plots/benchmarkHabitationModelAllTFs.pdf", width = 3, height = 4)
ggplot(plotData) +
  geom_col(aes(x = method, y = precision), width = 0.75,
           position = position_stack(reverse = TRUE)) +
  xlab("") + ylab("Precision") + ggtitle("All TFs") +
  theme_classic()
dev.off()

benchmarkResults <- as.data.frame(benchmarkResults)
benchmarkResults$diff <- benchmarkResults$MultiScale - benchmarkResults$Tobias
benchmarkResults$TF <- rownames(benchmarkResults)
rownames(TFClustering) <- TFClustering$TF
benchmarkResults$cluster <- as.character(TFClustering[benchmarkResults$TF,]$cluster)
benchmarkResults <- benchmarkResults[order(-benchmarkResults$diff),]
benchmarkResults$TF <- factor(benchmarkResults$TF, levels = benchmarkResults$TF)
ggplot(benchmarkResults, aes(x = TF, y = diff, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  ylab("Improvement") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

####################################
# Combine results and save to file #
####################################

# Combine results
performanceCompare <- data.frame(
  TF = rownames(benchmarkResults),
  TobiasPrecision = benchmarkResults$Tobias,
  MotifPrecision = benchmarkResults$Motif,
  HINTPrecision = benchmarkResults$HINT,
  nPred = benchmarkResults$nPred,
  multiScalePrecision = benchmarkResults$MultiScale,
  TFCluster = TFClustering[rownames(benchmarkResults),]$cluster
)
rownames(performanceCompare) <- performanceCompare$TF
classIPerformance <- performanceCompare[cluster1TFs,]
improvement <- classIPerformance$multiScalePrecision - classIPerformance$TobiasPrecision
classIPerformance <- classIPerformance[order(-improvement),]
plotData <- data.frame(
  TF = factor(rep(classIPerformance$TF, 4),
              levels = classIPerformance$TF),
  Precision = c(classIPerformance$multiScalePrecision,
                classIPerformance$HINTPrecision,
                classIPerformance$TobiasPrecision,
                classIPerformance$MotifPrecision),
  Method = factor(rep(c("Multi-scale", "HINT", "Tobias", "Motif"), 
                      each = dim(classIPerformance)[1]),
                  levels = c("Multi-scale", "HINT", "Tobias", "Motif"))
)
ggplot(data = plotData, aes(x = TF, y = Precision, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.75) +
  scale_fill_manual(values=c("#716148B2", "#8491B4B2", "darkorange", "#00A087FF")) +
  theme_classic()

write.table(performanceCompare, paste0("../../data/TFBSPrediction/benchmark",  projectName, "UnibindAllTFModel.txt"),
            quote = F, sep = "\t")
