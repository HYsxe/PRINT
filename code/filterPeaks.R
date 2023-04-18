source("../../code/utils.R")
library("data.table")
library(BuenRTools)
library(GenomicRanges)

# Get reference genome from command line arguments
args <-  commandArgs(trailingOnly=TRUE)
genome <- args[1]
print(genome)

# Path to summit file(s)
peakCallingFiles <- list.files("peakCalling")
summitFile <- peakCallingFiles[stringr::str_detect(peakCallingFiles, "summits")]

# Path to blacklist bed file that specifies regions to be filtered out.(for particular reference genome) 
# Also get chromosome size file
if(genome == "mm10"){
  blacklist <- "../shared/BlacklistFiles/mm10.full.blacklist.bed"
  chrSizes <- "../shared/mm10.chrom.sizes"
}else if(genome == "hg19"){
  blacklist <- "../shared/BlacklistFiles/hg19.full.blacklist.bed"
  chrSizes <- "../shared/hg19.chrom.sizes"
}else if(genome == "hg38"){
  blacklist <- "../shared/BlacklistFiles/hg38-blacklist.v2.bed"
  chrSizes <- "../shared/hg38.chrom.sizes"
}

# Call summitsToCleanPeaks function for peak clean-up
summitsToCleanPeaks(summitFiles = paste0("peakCalling/", summitFile),
                    peakWidth = 800, # Window to use for peak summit padding
                    blackList = blacklist,
                    chromSizes = chrSizes,
                    topNPeaks = NULL, # Filter top N peaks after?
                    useQuantileRanks = FALSE, # Use normalized ranks instead of raw MACS score?
                    FDRcutoff = 0.01, # MACS peak FDR cut-off used for all peaks
                    outDir = "peakCalling/", # Make sure this directory exists first
                    expName= "filtered", # Name to give as prefix to resulting filtered peak bed file
                    resizeTo = 1000 # Re-size peaks after filtering?
)

system("mv peakCalling/filtered.fixedwidthpeaks_800bp_1000reSized.bed peaks.bed")