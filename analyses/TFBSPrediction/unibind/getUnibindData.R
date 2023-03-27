# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(GenomicRanges)

dataset <- "GM12878"

# Get dataset metadata
ChIPDir <- "../../../data/shared/unibind/damo_hg38_all_TFBS/"
ChIPSubDirs <- list.files(ChIPDir)
ChIPMeta <- as.data.frame(t(sapply(
  ChIPSubDirs, 
  function(ChIPSubDir){
    strsplit(ChIPSubDir, "\\.")[[1]]
  }
)))
rownames(ChIPMeta) <- NULL
colnames(ChIPMeta) <- c("ID", "dataset", "TF")
ChIPMeta$Dir <- ChIPSubDirs

# Select cell types
datasetName <- list("K562" = "K562_myelogenous_leukemia",
                    "HepG2" = "HepG2_hepatoblastoma",
                    "GM12878" = "GM12878_female_B-cells_lymphoblastoid_cell_line")
ChIPMeta <- ChIPMeta[ChIPMeta$dataset == datasetName[[dataset]], ]
ChIPTFs <- ChIPMeta$TF

# Extract TFBS ChIP ranges
unibindTFBS <- pbmcapply::pbmclapply(
  sort(unique(ChIPTFs)),
  function(TF){
    
    # Retrieve the list of ChIP Ranges
    ChIPRangeList <- lapply(
      which(ChIPTFs %in% TF),
      function(entry){
        ChIPSubDir <- ChIPMeta$Dir[entry]
        ChIPFiles <- list.files(paste0(ChIPDir, ChIPSubDir))
        ChIPRanges <- lapply(
          ChIPFiles,
          function(ChIPFile){
            ChIPBed <- read.table(paste0(ChIPDir, ChIPSubDir, "/", ChIPFile))
            ChIPRange <- GRanges(seqnames = ChIPBed$V1, ranges = IRanges(start = ChIPBed$V2, end = ChIPBed$V3))
          }
        )
        ChIPRanges <- Reduce(subsetByOverlaps, ChIPRanges)
        ChIPRanges
      }
    )
    
    # For TFs with more than one ChIP files, take the intersection of regions in all files
    ChIPRanges <- Reduce(subsetByOverlaps, ChIPRangeList)
    ChIPRanges
  },
  mc.cores = 16
)
names(unibindTFBS) <- sort(unique(ChIPTFs))

# Save results to file
saveRDS(unibindTFBS, paste0("../../../data/shared/unibind/", dataset, "ChIPRanges.rds"))
