# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(GenomicRanges)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

source("../../code/utils.R")
source("../../code/getCounts.R")

####################################
# Get genomic ranges of BAC clones #
####################################

# Load list of BACs
selectedBACs <- read.table("../../data/BAC/rawData/selectedBACs.txt")
colnames(selectedBACs) <- c("ID", "chr", "start", "end")

# Get genomic ranges of BACs
BACRanges <- GRanges(
  seqnames = selectedBACs$chr,
  ranges = IRanges(start = selectedBACs$start,
                   end = selectedBACs$end)
)
names(BACRanges) <- selectedBACs$ID
saveRDS(BACRanges, "../../data/BAC/BACRanges.rds")

# Write BAC ranges to bed file
write.table(selectedBACs[,c("chr", "start", "end", "ID"),], "../../data/BAC/BACs.bed",
            quote = F, row.names = F, col.names = F, sep = "\t")

################################################################
# Get position-by-BAC-by-replicate ATAC insertion count tensor #
################################################################

# Initialize a footprintingProject object
projectName <- "BAC"
project <- footprintingProject(projectName = projectName, 
                               refGenome = "hg38")
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")
dataDir(project) <- projectDataDir
mainDir(project) <- projectMainDir

# Set the regionRanges slot
regionRanges(project) <- BACRanges

# Load barcodes for each replicate
barcodeGroups <- data.frame(barcode = paste("rep", 1:5, sep = ""),
                            group = 1:5)

# Down-sample fragments data
frags <- data.table::fread(paste0(projectDataDir, "rawData/all.fragments.tsv.gz"))
nFrags <- dim(frags)[1]
if(!dir.exists("../../data/BAC/downSampledFragments/")){
  system("mkdir ../../data/BAC/downSampledFragments")
}
for(downSampleRate in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
  print(paste0("Downsampling rat: ", downSampleRate))
  downSampleInd <- sample(1:nFrags, as.integer(nFrags * downSampleRate))
  downSampledFrags <- frags[downSampleInd, ]
  gz <- gzfile(paste0("../../data/BAC/downSampledFragments/fragmentsDownsample", downSampleRate, ".tsv.gz"), "w")
  write.table(downSampledFrags, gz, quote = F, row.names = F, col.names = F, sep = "\t")
  close(gz)
}

# Get position-by-BAC-by-replicate ATAC insertion count tensor for undownsampled and downsampled fragments
if(!file.exists("../../data/BAC/countTensors.rds")){
  counts <- list()
  pathToFrags <- paste0(projectDataDir, "rawData/all.fragments.tsv.gz")
  if(dir.exists("../../data/BAC/chunkedCountTensor")){
    system("rm -r ../../data/BAC/chunkedCountTensor")
  }
  counts[["all"]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  for(downSampleRate in c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){
    system("rm -r ../../data/BAC/chunkedCountTensor")
    pathToFrags <- paste0("../../data/BAC/downSampledFragments/fragmentsDownsample", downSampleRate, ".tsv.gz")
    counts[[as.character(downSampleRate)]] <- countTensor(getCountTensor(project, pathToFrags, barcodeGroups, returnCombined = T))
  }
  system("rm -r ../../data/BAC/chunkedCountTensor")
  saveRDS(counts, "../../data/BAC/countTensors.rds")
}else{
  counts <- readRDS("../../data/BAC/countTensors.rds")
}

######################################
# Get ground truth observed Tn5 bias #
######################################

# Load reference genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Parameters for bias estimation
localRadius <- 50
contextRadius <- 50
contextLen <- 2 * contextRadius + 1
localAverageThreshold <- 20

biasDf <- NULL
for(BACInd in 1:length(BACRanges)){
  
  print(paste0("Processing BAC No.", BACInd))
  
  observedBias <- list()
  localAverage <- list()
  ATACTracks  <- list()
  for(downSampleRate in names(counts)){
    
    BACCounts <- counts[[downSampleRate]][[BACInd]]
    BACCounts <- BACCounts %>% group_by(position) %>% summarize(insertion = sum(count))
    ATACTrack <- rep(0, width(BACRanges)[BACInd])
    ATACTrack[BACCounts$position] <- BACCounts$insertion
    
    # Calculate local average Tn5 insertion
    localAverage[[downSampleRate]] <- pbmcapply::pbmcmapply(
      function(i){
        mean(ATACTrack[max(1, i - localRadius):min(length(ATACTrack), i + localRadius)])
      },
      1:length(ATACTrack),
      mc.cores = 12
    )
    
    # Calculate observed Tn5 bias = observed insertion / local coverage
    observedBias[[downSampleRate]] <- ATACTrack / localAverage[[downSampleRate]]
    
    ATACTracks[[downSampleRate]] <- ATACTrack
  }
  
  localAverage <- Reduce(cbind, localAverage)
  observedBias <- Reduce(cbind, observedBias)
  ATACTracks <- Reduce(cbind, ATACTracks)
  colnames(localAverage) <- paste0("localAverage_", names(counts))
  colnames(observedBias) <- paste0("obsBias_", names(counts)) 
  colnames(ATACTracks) <- paste0("ATAC_", names(counts)) 
  
  # Only select positions with enough numbers of observations
  selectedPositions <- localAverage[, "localAverage_all"] > localAverageThreshold
  
  # Get sequence of the peak (inclueding left and right padding rgions)
  seqContext <- Biostrings::getSeq(hg38, as.character(seqnames(BACRanges[BACInd])), 
                                   start = start(BACRanges[BACInd]) - contextRadius, 
                                   width = width(BACRanges[BACInd]) + contextLen - 1,
                                   as.character = T)
  
  # Get sequence context for every position
  positions <- 1:width(BACRanges[BACInd])
  context <- substring(seqContext, positions, positions + contextLen - 1)
  
  # Return results
  if(sum(selectedPositions) > 0){
    BACResults <- Reduce(cbind, list(data.frame(context = context[selectedPositions]), 
                                     as.data.frame(observedBias[selectedPositions,]),
                                     as.data.frame(localAverage[selectedPositions,]),
                                     as.data.frame(ATACTracks[selectedPositions,]),
                                     data.frame(BACInd = BACInd)))
    biasDf <- rbind(biasDf, BACResults) 
  }
  
}

write.table(biasDf, 
            paste0("../../data/BAC/obsBias.tsv"), 
            quote = F, row.names = F, sep = "\t")

#################################
# Evaluate bias reproducibility #
#################################

nReplicates <- 5
coverageThreshold <- 20

replicateBiasMtx <- NULL
for(BACInd in 1:length(BACRanges)){
  
  print(paste0("Processing BAC No.", BACInd))
  
  BACCounts <- counts[["all"]][[BACInd]]
  
  # Get ATAC track for the current BAC and each replicate
  replicateTracks <- sapply(
    1:nReplicates,
    function(replicateInd){
      replicateBACCounts <- BACCounts %>% filter(group == replicateInd)
      ATACTrack <- rep(0, width(BACRanges)[BACInd])
      ATACTrack[replicateBACCounts$position] <- replicateBACCounts$count
      ATACTrack
    }
  )
  
  # Calculate local average ATAC
  replicateCoverage <- sapply(
    1:nReplicates,
    function(replicateInd){
      conv(replicateTracks[, replicateInd], localRadius) / (2 * localRadius)
    }
  )
  
  # Calculate observed Tn5 bias for each replicate
  replicateBias <- replicateTracks / replicateCoverage
  
  # Filter out low-coverage regions
  replicateBiasFilt <- replicateBias[rowSums(replicateCoverage) > coverageThreshold,]
  
  replicateBiasMtx <- rbind(replicateBiasMtx, replicateBiasFilt)
}
saveRDS(replicateBiasMtx, "../../data/BAC/replicateBiasMtx.rds")
replicateBiasMtx <- readRDS("../../data/BAC/replicateBiasMtx.rds")

# Visualize reproducibility between two replicates
plotData <- data.frame("Replicate1" = replicateBiasFilt[, 1],
                       "Replicate2" = replicateBiasFilt[, 2])
ggpubr::ggscatter(data = plotData,
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  x = "Replicate1", xlab = "ATAC insertion count (replicate 1)",
                  y = "Replicate2", ylab = "ATAC insertion count (replicate 2)",
                  title = paste0("BAC ", selectedBACs$ID[BACInd]),
                  size = 0.1)

# Visualize reproducibility among all 5 replicates
corMat <- cor(replicateBiasMtx)
col_fun = colorRamp2(c(0.97, 0.99), c("white", "red"))
rownames(corMat) <- paste("Rep", 1:nReplicates)
colnames(corMat) <- paste("Rep", 1:nReplicates)
Heatmap(corMat,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", corMat[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        use_raster = TRUE,
        cluster_rows = F,
        cluster_columns = F,
        name = "Correlation",
        border = TRUE,
        column_names_rot = 0,
        column_title_side = "bottom")