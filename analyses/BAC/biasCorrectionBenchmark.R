# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getGroupData.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
library(GenomicRanges)
library(rtracklayer)
library(ggpubr)

#################
# Load BAC data #
#################

# Load genomic ranges of BACs
BACRanges <- readRDS("../../data/BAC/BACRanges.rds")

# Convert BAC genomic ranges into 1kb tiles
tileRanges <- Reduce("c", GenomicRanges::slidingWindows(BACRanges, width = 1000, step = 1000))
tileRanges <- tileRanges[width(tileRanges) == 1000]
tileBACs <- names(tileRanges)

# Load reference genome
hg38 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

# Load BAC tile count tensor
countData <- readRDS("../../data/BAC/tileCounts.rds")$all
names(countData) <- 1:length(countData)

###############################################################
# Get CNN-predicted Tn5 bias and observed Tn5 insertion track #
###############################################################

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

###############
# Get Tn5 PWM #
###############

# One-hot encoding of a single DNA sequence
onehotEncode <- function(seq){
  onehotMapping <- list("A" = 1, "C" = 2, "G" = 3, "T" = 4, "N" = 0)
  len <- nchar(seq)
  onehotInd <- sapply(
    substring(seq, 1:len , 1:len),
    function(x){onehotMapping[[x]]})
  onehot <- matrix(0, nrow = 4, ncol = len)
  undeterminedBases <- onehotInd == 0
  onehot[cbind(unname(onehotInd), 1:len)[!undeterminedBases,]] <- 1
  onehot
}

# Use PWM to score a sequence
PWMScoring <- function(seq, PWM){
  onehotSeq <- onehotEncode(seq)
  sum(PWM * onehotSeq)
}

# Get position frequency matrix of Tn5 insertion
PWMRadius <- 10
tilePFM <- pbmcapply::pbmclapply(
  1:length(tileRanges),
  function(tileInd){
    tileRange <- tileRanges[tileInd]
    tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
    contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
    contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
    ATACTrack <- rowSums(getRegionATAC(countData, 
                                       regionInd = as.character(tileInd), 
                                       groupIDs = 1:5, 
                                       width = 1000))
    
    # One-hot encoding of each sequence context that appeared in the BAC
    # This is used to quantify the background sequence frequencies
    backgroundOnehot <- lapply(
      1:length(contextSeq), 
      function(i){onehotEncode(contextSeq[i])})
    backgroundPFM <- Reduce("+", backgroundOnehot)
    
    # Since each context is cut with different frequencies, we calculate
    # the foreground frequencies by weighing the background with cutting density 
    foregroundOnehot <- lapply(
      1:length(contextSeq), 
      function(i){backgroundOnehot[[i]] * ATACTrack[i]})
    foregroundPFM <- Reduce("+", foregroundOnehot)
    
    rownames(foregroundPFM) <- c("A", "C", "G", "T")
    rownames(backgroundPFM) <- c("A", "C", "G", "T")
    
    list(foregroundPFM, backgroundPFM)
  }
)
fgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[1]]}))
bgPFM <- Reduce("+", lapply(tilePFM, function(x){x[[2]]}))

# Get background base frequencies
background <- rowSums(bgPFM)
background <- background / sum(background)

# Get PPM 
PPM <- t(t(fgPFM) / colSums(fgPFM))
adjustedPPM <- PPM / background
PWM <- log2(adjustedPPM)

library(Logolas)
system("mkdir ../../data/BAC/plots")
pdf("../../data/BAC/plots/Tn5PWM.pdf", width = 8, height = 6)
logomaker(adjustedPPM, type = "Logo", return_heights = TRUE)
dev.off()

# Save the PWM to file
saveRDS(PWM, "../../data/BAC/Tn5PWM.rds")

############################
# Get k-mer Tn5 bias model #
############################

kmerBiasList <- list()
for(k in c(3,5,7)){
  
  # Generate all possible kmer sequences
  bases <- c("A", "C", "G", "T")
  kmers <- as.matrix(eval(parse(text = paste("expand.grid(", paste(rep("bases", k), collapse = ", "), ")"))))
  kmers <- sapply(1:dim(kmers)[1], function(x){paste(kmers[x,], collapse = "")})
  
  # Go through each genomic tile and record insertion numbers for kmers
  tileKmer <- pbmcapply::pbmclapply(
    1:length(tileRanges),
    function(tileInd){
      tileRange <- tileRanges[tileInd]
      tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
      contextRanges <- IRanges::resize(tilePositions, fix = "center", width = k)
      contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
      ATACTrack <- rowSums(getRegionATAC(countData, 
                                         regionInd = as.character(tileInd), 
                                         groupIDs = 1:5, 
                                         width = 1000))
      kmerFG <- rep(0, length(kmers))
      names(kmerFG) <- kmers
      kmerBG <- rep(0, length(kmers))
      names(kmerBG) <- kmers
      kmerFG[contextSeq] <- ATACTrack
      kmerBG[contextSeq] <- 1
      list(kmerFG, kmerBG)
    }
  )
  
  # Summarize results for all tiles and calculate k-mer bias
  kmerFG <- rowSums(sapply(tileKmer,function(x){x[[1]]}))
  kmerBG <- rowSums(sapply(tileKmer,function(x){x[[2]]}))
  kmerBias <- kmerFG / kmerBG
  kmerBiasList[[as.character(k)]] <- kmerBias / mean(kmerBias)
}

saveRDS(kmerBiasList, "../../data/BAC/kmerBiasList.rds")

#############################
# Benchmark on all BAC data #
#############################

kmerBiasList <- readRDS("../../data/BAC/kmerBiasList.rds")
PWM <- readRDS("../../data/BAC/Tn5PWM.rds")
TobiasBigWig <- import.bw('../../data/BAC/Tobias/BAC_bias.bw')
PWMRadius <- 10

benchmark <- t(pbmcapply::pbmcmapply(
  function(tileInd){
    
    tileRange <- tileRanges[tileInd]
    tileInsertion <- rowSums(getRegionATAC(countData, regionInd = as.character(tileInd), 
                                           groupIDs = 1:5, width = 1000))
    tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
    coverage <- sum(tileInsertion)
    
    biasBenchmark <- list()
    
    # Retrieve Tn5 bias calculated by kmer
    for(k in names(kmerBiasList)){
      contextRanges <- IRanges::resize(tilePositions, fix = "center", width = as.integer(k))
      contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
      kmerBias <-  unname(kmerBiasList[[k]][contextSeq])
      biasBenchmark[[paste0(k, "-mer")]] <- cor(kmerBias, tileInsertion)
    }
    
    # Retrieve Tn5 bias calculated by PWM
    contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
    contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
    PWMBias  <-  2 ^ sapply(contextSeq, function(seq){PWMScoring(seq, PWM)})
    biasBenchmark[["PWM"]] <- cor(PWMBias, tileInsertion)
    
    # Retrieve Tn5 bias estimated by Tobias
    TobiasTn5Bias <- subsetByOverlaps(TobiasBigWig, tileRanges[tileInd])$score
    biasBenchmark[["dinucleotidePWM"]] <- cor(2 ^ TobiasTn5Bias, tileInsertion)
    
    # Calculate accuracy by convolution neural net
    biasBenchmark[["ConvNet"]] <- cor(regionBias(project)[tileInd,], tileInsertion / conv(tileInsertion, 50))
    
    biasBenchmark <- Reduce(c, biasBenchmark)
    
    c(coverage, biasBenchmark)
    
  },
  1:length(tileRanges),
  mc.cores = 16
))

# Filter out entries with NA values
benchmark <- benchmark[rowSums(is.na(benchmark)) == 0, ]

# Filter out low-coverage tiles
benchmark <- benchmark[benchmark[, 1] > 20000, ]

methods <- c("3-mer", "5-mer", "7-mer", "PWM", "dinucleotide PWM", "CNN")
plotData <- data.frame(
  Correlation = colMeans(benchmark)[2:7],
  Method = methods
)
pdf("../../data/BAC/plots/benchmark.pdf", width = 8, height = 8)
ggplot(plotData) +
  geom_bar(aes(x = Method, y = Correlation), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Bias correction method") + ylim(0, 1) + 
  ylab("Correlation between predicted \nand observed Tn5 bias") + 
  scale_x_discrete(limits = methods) +
  theme_classic()
dev.off()

########################
# Model interpretation #
########################

# Get Tn5 bias predicted by PWM and CNN
predBias <- t(pbmcapply::pbmclapply(
  1:length(tileRanges),
  function(tileInd){
    
    tileRange <- tileRanges[tileInd]
    tileInsertion <- rowSums(getRegionATAC(countData, regionInd = as.character(tileInd), 
                                           groupIDs = 1:5, width = 1000))
    
    if(mean(tileInsertion) > 20){
      tilePositions <- IRanges::tile(tileRange, width = 1)[[1]]
      localCoverage <- conv(tileInsertion, 50) / 100
      obsBias <- tileInsertion / localCoverage
      
      # Retrieve Tn5 bias calculated by PWM
      contextRanges <- IRanges::resize(tilePositions, fix = "center", width = PWMRadius * 2 + 1)
      contextSeq <- Biostrings::getSeq(hg38, contextRanges, as.character = T)
      PWMBias  <-  2 ^ sapply(contextSeq, function(seq){PWMScoring(seq, PWM)})
      
      # Retrieve Tn5 bias estimated by CNN
      CNNBias <- regionBias(project)[tileInd,]
      
      data.frame(PWMBias = PWMBias, CNNBias = CNNBias, obsBias = obsBias, context = contextSeq)
    }else{
      NULL
    }
    
  },
  mc.cores = 16
))
predBias <- data.table::rbindlist(predBias)
predBias <- predBias[!is.na(predBias$obsBias), ]

# Visual inspection
sampleInd <- sample(1:dim(predBias)[1], 10000)
qplot(predBias$PWMBias[sampleInd], predBias$obsBias[sampleInd]) + xlim(0,50) + ylim(0,50)

# Calculate errors by the two methods on individual positions
PWMError <- log2(pmax(predBias$PWMBias, 1e-3)) - log2(pmax(predBias$obsBias, 1e-3))
CNNError <- log2(pmax(predBias$CNNBias, 1e-3)) - log2(pmax(predBias$obsBias, 1e-3))
PWMAbsError <- abs(PWMError)
CNNAbsError <- abs(CNNError)
improvement <- PWMAbsError - CNNAbsError

qplot(PWMError[sampleInd], improvement[sampleInd], )
hist(PWMError[order(improvement, decreasing = T)][1:2000])

diffContext <- predBias$context[order(improvement, decreasing = T)][1:2000]
simContext <- predBias$context[order(abs(improvement))][1:2000]
diffCGContent <- sapply(diffContext, function(x){stringr::str_count(x, "[CG]")}) / (2 * PWMRadius + 1)
simCGContent <- sapply(simContext, function(x){stringr::str_count(x, "[CG]")}) / (2 * PWMRadius + 1)

plotData <- diffCGContent
data.frame(plotData) %>%
  ggplot(aes(x = plotData)) +
  geom_histogram(bins=20, fill = '#9C5D41') +
  xlab("GC content") + ylab("Number of examples") +
  theme_classic()
