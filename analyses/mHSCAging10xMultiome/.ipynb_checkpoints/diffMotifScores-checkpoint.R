# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getTFBS.R")
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(ggrepel)

###################
# Load input data #
###################

# Load footprinting project
project <- readRDS("../../data/mHSCAging10xMultiome/project.rds")

regions <- regionRanges(project)

# Load LSI embedding of single cells
cellEmbedding <- readRDS("../../data/mHSCAging10xMultiome/LSIEmbedding.rds")

# Load pseudobulk centers
pseudobulkCenters <- read.table("../../data/mHSCAging10xMultiome/pseudobulkCenters.txt",
                                header = T)
rownames(pseudobulkCenters) <- pseudobulkCenters$group

# Get substructure-by-pseudobulk SummarizedExperiment object
substructurePath <- "../../data/mHSCAging10xMultiome/substructureSE.rds"
if(file.exists(substructurePath)){
  substructureSE <- readRDS("../../data/mHSCAging10xMultiome/substructureSE.rds")
}else{
  substructureSE <- getSubstructureSE(project)
  saveRDS(substructureSE, substructurePath)
}
substructureRanges <- rowRanges(substructureSE)
substructureMat <- assay(substructureSE)

# Load results from differential RNA testing
diffRNA <- read.table("../../data/mHSCAging10xMultiome/diffRNA.tsv")

# Load differential TFBS and CRE results
diffCompare <- read.table("../../data/mHSCAging10xMultiome/diff_CRE_substructure_compare.tsv",
                          sep = "\t", header = T)

# Load single cell ATAC data
scATACSeurat <- readRDS("../../data/mHSCAging10xMultiome/scATACSeurat.rds")

# Load TF motifs
cisBPMotifs <- readRDS("../..//data/shared/cisBP_mouse_pwms_2021.rds")

#########################
# Re-number pseudobulks #
#########################

pseudobulkNames <- groups(project)
pseudobulkCenters <- pseudobulkCenters[pseudobulkNames,]
pseudobulkAges <- sapply(pseudobulkCenters$barcode, function(x){strsplit(x, "-")[[1]][2]})
pseudobulkLSI <- cellEmbedding[pseudobulkCenters$barcode, 2]

# Within each age-group, re-order pseudobulks by LSI-2
reOrder <- c(order(pseudobulkLSI[pseudobulkAges == "Old"]), 
             order(pseudobulkLSI[pseudobulkAges == "Young"]) + sum(pseudobulkAges == "Old"))
substructureMat <- substructureMat[, reOrder]
groupATAC(project) <- groupATAC(project)[, reOrder]
groupRNA(project) <- groupRNA(project)[, reOrder]

# Re-label pseudobulks
colnames(substructureMat) <- pseudobulkNames
colnames(groupATAC(project)) <- pseudobulkNames
colnames(groupRNA(project)) <- pseudobulkNames

#######################################
# Old/Young-like old/Young assignment #
#######################################

seuratClusters <- as.character(scATACSeurat$seurat_clusters)
names(seuratClusters) <- colnames(scATACSeurat)
ageCluster <- seuratClusters[pseudobulkCenters$barcode]
names(ageCluster) <- pseudobulkCenters$group
ageCluster <- ageCluster[pseudobulkNames]
ageCluster <- ageCluster[reOrder]
names(ageCluster) <- pseudobulkNames
ageCluster[ageCluster == "0"] <- "Old"
ageCluster[ageCluster == "6"] <- "Young-like old"
ageCluster[ageCluster == "10"] <- "Young"

################################
# Motif scoring of pseudobulks #
################################

# Filter out sites with low signal
substructureFilter <- rowMaxs(substructureMat) > 0.3
substructureRanges <- substructureRanges[substructureFilter]
substructureMat <- substructureMat[substructureFilter,]

# Create SummarizedExperiment object for differential testing
mode <- "substructure"
if(mode == "substructure"){
  diffInds <- (abs(diffCompare$substructureDiff) > 1) & (abs(diffCompare$CREDiff) < 1)
  diffInds <- diffCompare$substructureInd[diffInds]
  diffSE <- SummarizedExperiment(list(counts = substructureMat[diffInds, ]),
                                 rowRanges = substructureRanges[diffInds])
}else if(mode == "CRE"){
  diffInds <- (abs(diffCompare$substructureDiff) < 1) & (abs(diffCompare$CREDiff) > 1)
  diffInds <- diffCompare$CREInd[diffInds]
  diffSE <- SummarizedExperiment(list(counts = groupATAC(project)[diffInds, ]),
                                 rowRanges = regions[diffInds])
}

# Filter features with zero signal
diffSE <- diffSE[rowSums(assay(diffSE)) > 0, ]

# Add GC bias
diffSE <- chromVAR::addGCBias(diffSE, genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)

# Remove NA values
bias <- rowRanges(diffSE)$bias
rowRanges(diffSE)$bias[is.na(rowRanges(diffSE)$bias)] <- mean(bias[!is.na(bias)])

# Get background substructures or CREs with matched GC bias and average signal
bgSites <- chromVAR::getBackgroundPeaks(diffSE,niterations = 100) 

# Get pseudobulk-by-TF motif score matrix
motifScores <- pbmcapply::pbmcmapply(
  function(TF){
    
    # Find motif matches
    motifMatches <- motifmatchr::matchMotifs(cisBPMotifs[TF], 
                                             diffSE, 
                                             genome = "mm10")
    
    # Compute motif scores (deviation score from background)
    deviation <- chromVAR::computeDeviations(object = diffSE, 
                                             annotations = motifMatches, 
                                             background_peaks = bgSites)
    deviation <- assay(deviation)
  },
  names(cisBPMotifs),
  mc.cores = 8
)
rownames(motifScores) <- colnames(diffSE)
saveRDS(motifScores, paste0("../../data/mHSCAging10xMultiome/", mode, "MotifScores.rds"))

# Differential testing comparing motifs scores of groupA and groupB
groupsA <- which(stringr::str_detect(rownames(motifScores), "Young"))
groupsB <- which(stringr::str_detect(rownames(motifScores), "Old"))
diffMotifScore <- pbmcapply::pbmclapply(
  names(cisBPMotifs),
  function(TF){
    test <- t.test(motifScores[groupsA, TF], motifScores[groupsB, TF])
    data.frame(pval = test$p.value, diffMean = mean(motifScores[groupsB, TF]) - mean(motifScores[groupsA, TF]))
  }
)
diffMotifScore <- as.data.frame(data.table::rbindlist(diffMotifScore))
rownames(diffMotifScore) <- names(cisBPMotifs)
diffMotifScore$TF <- names(cisBPMotifs)
diffMotifScore$logP <- -log10(diffMotifScore$pval)
diffMotifScore$fdr <- p.adjust(diffMotifScore$pval, method = "fdr")
diffMotifScore$signedLogFDR <- -log10(diffMotifScore$fdr) * sign(diffMotifScore$diffMean)
diffMotifScore <- diffMotifScore[order(diffMotifScore$signedLogFDR, decreasing = T),]
saveRDS(diffMotifScore, paste0("../../data/mHSCAging10xMultiome/",
                               mode, "DiffMotifScores.rds"))

######################################################
# Comapre diff motif scores with diff TF RNA results #
######################################################

ovTFs <- intersect(rownames(diffMotifScore), rownames(diffRNA))
diffMotifScore <- diffMotifScore[ovTFs, ]
diffTFRNA <- diffRNA[ovTFs, ]
diffMotifScore$FDR <- p.adjust(diffMotifScore$pval, method = "fdr")
diffTFRNA$FDR <- p.adjust(diffTFRNA$pvalue, method = "fdr")
diffMotifScore$signedLogFDR <- -log10(diffMotifScore$FDR) * sign(diffMotifScore$diffMean)
diffTFRNA$signedLogFDR <- -log10(diffTFRNA$FDR) * sign(diffTFRNA$log2FoldChange)

plotData <- data.frame(
  diffMotifScore = diffMotifScore$signedLogFDR,
  diffTFRNA = diffTFRNA$signedLogFDR,
  TF = ovTFs
)
selected <- ((plotData$diffMotifScore > 1) & (plotData$diffTFRNA > 10)) |
  ((plotData$diffMotifScore < -1) & (plotData$diffTFRNA < -10))
rownames(plotData) <- ovTFs
labels <- rep("", dim(plotData)[1])
labels[selected] <- plotData$TF[selected]
pdf(paste0("../../data/mHSCAging10xMultiome/plots/Old_HSC_vs_Young_HSC_diff_", 
           mode, "_motif_scores_RNA.pdf"),
    width = 7, height = 6)
ggplot(plotData) +
  geom_point(aes(x = diffMotifScore, y = diffTFRNA), color = "grey") +
  geom_point(data = plotData[selected,],
             aes(x = diffMotifScore, y = diffTFRNA), color = "red") +
  geom_text_repel(x = plotData$diffMotifScore, y = plotData$diffTFRNA, 
                  label = labels, max.overlaps = 1e5,
                  force = 1) +
  xlab("Differential motif score (Old - Young)\nsigned -log10(FDR)") +
  ylab("Differential RNA (Old - Young)\nsigned -log10(FDR)") +
  theme_classic()
dev.off()

######################################################
# Comapre diff motif scores with diff TF RNA results #
######################################################

diffCREMotifScores <- readRDS("../../data/mHSCAging10xMultiome/CREDiffMotifScores.rds")
diffSubstructureMotifScores <- readRDS("../../data/mHSCAging10xMultiome/substructureDiffMotifScores.rds")
ovTFs <- intersect(rownames(diffCREMotifScores), intersect(rownames(diffCREMotifScores), rownames(diffRNA)))
diffTFRNA <- diffRNA[ovTFs, ]
ovTFs <- ovTFs[diffTFRNA$padj < 0.01]
diffCREMotifScores <- diffCREMotifScores[ovTFs, ]
diffSubstructureMotifScores <- diffSubstructureMotifScores[ovTFs, ]

diffRNASigndLogFDR <- sign(diffTFRNA[ovTFs, ]$log2FoldChange) * -log10(diffTFRNA[ovTFs, ]$padj)
diffRNASigndLogFDR <- pmin(pmax(diffRNASigndLogFDR, -5), 5)
plotData <- data.frame(diffSubstructure = diffSubstructureMotifScores$signedLogFDR,
                       diffCRE = diffCREMotifScores$signedLogFDR,
                       TF = ovTFs,
                       diffRNA = diffRNASigndLogFDR)
rownames(plotData) <- ovTFs

# Label substructure-specific TFs
selected <- ((abs(diffSubstructureMotifScores$signedLogFDR) > 2) & (abs(diffCREMotifScores$signedLogFDR) < 1)) 
labelsA <- rep("", dim(plotData)[1])
labelsA[selected] <- plotData$TF[selected]

# Label shared TFs
selected <- ((abs(diffSubstructureMotifScores$signedLogFDR) > 2) & 
             (abs(diffCREMotifScores$signedLogFDR) > 2) &
               (diffCREMotifScores$signedLogFDR * diffSubstructureMotifScores$signedLogFDR > 0)) 
rownames(plotData) <- ovTFs
labelsB <- rep("", dim(plotData)[1])
labelsB[selected] <- plotData$TF[selected]

pdf("../../data/mHSCAging10xMultiome/plots/CRE_subtructure_diff_motif_score_compare.pdf",
    width = 7, height = 6)
ggplot(plotData) +
  geom_point(aes(x = diffCRE, y = diffSubstructure, color = diffRNA)) +
  xlim(-10,10) + ylim(-10,10) + 
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  xlab("Diff CRE motif scores signed log10(FDR)") +
  ylab("Diff substructure motif scores signed log10(FDR)") +
  geom_text_repel(x = plotData$diffCRE, y = plotData$diffSubstructure, 
                  label = labelsA, max.overlaps = 1e5,
                  force = 1, color = "Red") +
  geom_text_repel(x = plotData$diffCRE, y = plotData$diffSubstructure, 
                  label = labelsB, max.overlaps = 1e5,
                  force = 1, color = "Black") +
  geom_vline(xintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_vline(xintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  theme_classic()
dev.off()

###########################################
# Plot heatmap of TF motif scores and RNA #
###########################################

plotTFs <- "Specific"

if(plotTFs == "Shared"){
  # Select shared TFs
  selected <- ((abs(diffSubstructureMotifScores$signedLogFDR) > 2) & 
                 (abs(diffCREMotifScores$signedLogFDR) > 2) &
                 (diffCREMotifScores$signedLogFDR * diffSubstructureMotifScores$signedLogFDR > 0) &
                 (diffTFRNA[ovTFs, ]$log2FoldChange * diffSubstructureMotifScores$signedLogFDR > 0)) 
  
}else if(plotTFs == "Specific"){
  # Select substructure-specific TFs
  selected <- ((abs(diffSubstructureMotifScores$signedLogFDR) > 2) & (abs(diffCREMotifScores$signedLogFDR) < 1)) 
  
}

# Visualize motif scores per pseudobulk
pltMatrix <- motifScores[,plotData$TF[selected]]
pltMatrix <- t(pltMatrix)
pltMatrix <- t((pltMatrix - rowMeans(pltMatrix)) / rowSds(pltMatrix))
TFOrder <- rev(hclust(distCosine(t(motifScores[,plotData$TF[selected]])))$order)
pdf(paste0("../../data/mHSCAging10xMultiome/plots/", mode, plotTFs, "PseudobulkMotifScores.pdf"),
    width = 7, height = 6)
Heatmap(pltMatrix[, TFOrder],
        cluster_rows = F,
        cluster_columns = F,
        name = "Motif\nscore",
        row_split = factor(ageCluster, levels = c("Old", "Young-like old", "Young")),
        col = jdb_palette(n = 9, "solar_extra")[1:9])
dev.off()

# Visualize RNA of the corresponding TFs per pseudobulk
pltMatrix <- t(groupRNA(project)[plotData$TF[selected], ])
pltMatrix <- t(pltMatrix)
pltMatrix <- t((pltMatrix - rowMeans(pltMatrix)) / rowSds(pltMatrix))
pdf("../../data/mHSCAging10xMultiome/plots/pseudobulkTFRNA.pdf",
    width = 7, height = 6)
Heatmap(pltMatrix[, TFOrder],
        cluster_rows = F,
        cluster_columns = F,
        row_split = factor(ageCluster, levels = c("Old", "Young-like old", "Young")),
        name = "RNA",
        col = jdb_palette(n = 9, "solar_extra")[1:9])
dev.off()
