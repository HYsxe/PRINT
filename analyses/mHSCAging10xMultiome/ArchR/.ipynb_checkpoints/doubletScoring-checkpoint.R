# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(Seurat)
library(Signac)
library(SummarizedExperiment)
library(ggplot2)
library(ArchR)

###########################################
# Compute doublet scores for single cells #
###########################################

addArchRThreads(threads = 16) 
addArchRGenome("mm10")

# Reformat 
reformatFragmentFiles(fragmentFiles = "../../../data/mHSCAging10xMultiome/all.frags.filt.tsv.gz")
inputFiles <- "../../../data/mHSCAging10xMultiome/all.frags.filt-Reformat.tsv.gz"
names(inputFiles) <- "aging"

# Create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 300, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  threads = 1
)

# Calculate doublet scores for each cell
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1,
  threads = 16
)

# Save results
doubletEnrichment <- doubScores[[1]]$doubletEnrich
names(doubletEnrichment) <- sapply(names(doubletEnrichment), function(x){strsplit(x, "#")[[1]][2]})
saveRDS(doubletEnrichment, "../../../data/mHSCAging10xMultiome/doubletEnrichment.rds")
