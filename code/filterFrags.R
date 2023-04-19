source("../../code/utils.R")
library("data.table")
library(GenomicRanges)
library(dplyr)
library(SummarizedExperiment)

# Get ATAC peak ranges
peakBed <- read.table("peaks.bed", sep = "\t", header = F)
peaks <- GRanges(seqnames = peakBed$V1,
                 ranges = IRanges(start = peakBed$V2, end = peakBed$V3))

# Get peak-by-cell count matrix
scATAC <- getCountsFromFrags("all.frags.tsv.gz", peaks)

# Filter by FRIP and depth
FRIPThreshold <- 0.3
depthThreshold <- 100
scATACFilt <- scATAC[, (scATAC$FRIP > FRIPThreshold) & (scATAC$FRIP < 1) & (scATAC$depth > depthThreshold)]

# Save scATAC data after filtering
saveRDS(scATACFilt, "scATAC.rds")
write.table(colnames(scATACFilt), "scBarcodes.txt", sep = "\t", quote = F,
            row.names = F, col.names = F)

# Filter the fragments file
frags <- data.table::fread("all.frags.tsv.gz", sep = "\t", showProgress = TRUE) %>% data.frame() 
fragsFilt <- frags %>% filter(V4 %in% colnames(scATACFilt))

# Save results to file
fragsGz <- gzfile("all.frags.filt.tsv.gz", "w")
write.table(fragsFilt, fragsGz, quote = F, sep = "\t", col.names = F, row.names = F)
close(fragsGz)