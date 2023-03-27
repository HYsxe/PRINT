library(BuenRTools)
library(GenomicRanges)
library(dplyr)

# Label the original hg19 files with "hg19"
system("cp peaks.bed hg19Peaks.bed")
system("mv all.frags.filt.tsv.gz all.frags.filt.hg19.tsv.gz")

# Load hg19 peaks
hg19PeakBed <- read.csv("hg19Peaks.bed", sep = "\t", header = F)
hg19Peaks <- GRanges(seqnames = hg19PeakBed$V1,
                     ranges = IRanges(start =  hg19PeakBed$V2,
                                      end =  hg19PeakBed$V3))

# Load chain file for liftOver
pathToChain <- "../../data/shared/hg19ToHg38.over.tab.chain"
ch <- rtracklayer::import.chain(pathToChain)

# Lift over from hg19 to hg38
seqlevelsStyle(hg19Peaks) <- "UCSC"  # necessary
hg38Peaks <- unlist(rtracklayer::liftOver(hg19Peaks, ch))

# Remove fragmented peaks and resize all peaks to 1kb
hg38Peaks <- hg38Peaks[width(hg38Peaks) >= 1000]
hg38Peaks <- resize(hg38Peaks, width = 1000, fix = "center")

# Get fragments aligned to hg19
hg19Frags <- fragsToRanges("all.frags.filt.hg19.tsv.gz")
seqlevelsStyle(hg19Frags) <- "UCSC"  # necessary
hg38Frags <- unlist(rtracklayer::liftOver(hg19Frags, ch))

# fragsToRanges converts fragments file (0-based indexing) to GRanges (1-based indexing)
# Now we convert it back
start(hg38Frags) <- start(hg38Frags) - 1 
hg38Frags <- data.frame(chr = as.character(seqnames(hg38Frags)),
                        start = start(hg38Frags),
                        end = end(hg38Frags),
                        barcode = hg38Frags$V4)

# Save results to files
hg38PeaksBed <- data.frame(chr = as.character(seqnames(hg38Peaks)),
                           start = start(hg38Peaks),
                           end = end(hg38Peaks))
write.table(hg38PeaksBed, "peaks.bed", quote = F, sep = "\t", col.names = F, row.names = F)
saveRDS(hg38Peaks, "peakRanges.rds")
fragsGz <- gzfile("all.frags.filt.tsv.gz", "w")
write.table(hg38Frags, fragsGz, quote = F, sep = "\t", col.names = F, row.names = F)
close(fragsGz)