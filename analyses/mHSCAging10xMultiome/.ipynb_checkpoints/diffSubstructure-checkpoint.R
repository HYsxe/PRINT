# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getFootprints.R")
source("../../code/getSubstructures.R")
source("../../code/getTFBS.R")
source("../../code/visualization.R")
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

###################
# Load input data #
###################

# Load footprinting project
project <- readRDS("../../data/mHSCAging10xMultiome/project.rds")

regions <- regionRanges(project)

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

# Load LSI embedding of single cells
cellEmbedding <- readRDS("../../data/mHSCAging10xMultiome/LSIEmbedding.rds")

# Load pseudobulk centers
pseudobulkCenters <- read.table("../../data/mHSCAging10xMultiome/pseudobulkCenters.txt",
                                header = T)
rownames(pseudobulkCenters) <- pseudobulkCenters$group

# Load results from differential RNA testing
diffRNA <- read.table("../../data/mHSCAging10xMultiome/diffRNA.tsv")

# Load single cell ATAC data
scATACSeurat <- readRDS("../../data/mHSCAging10xMultiome/scATACSeurat.rds")

barcodeGroups <- read.table("../../data/mHSCAging10xMultiome/barcodeGrouping.txt",
                            header = T)

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

########################
# Differential testing #
########################

# Filter out sites with low signal
substructureFilter <- rowMaxs(substructureMat) > 0.3
substructureRanges <- substructureRanges[substructureFilter]
substructureMat <- substructureMat[substructureFilter,]

# Differential testing for CREs and candidate TF binding sites
CREMat <- groupATAC(project)
featureMatList <- list("CRE" = CREMat, "substructure" = substructureMat)
diffResults <- list("CRE" = list(), "substructure" = list())
for(feature in c("substructure", "CRE")){
  
  featureMat <- featureMatList[[feature]]
  
  pbulkAge <- unname(sapply(colnames(featureMat), function(x){strsplit(x, "_")[[1]][1]}))
  
  # Calculate mean feature scores per site
  featureMeans <- rowMeans(featureMat)
  
  # Calculate the difference between young and old
  if(feature == "CRE"){
    featureDiff <- log2(rowMeans(featureMat[, pbulkAge == "Old"] + 0.01) / rowMeans(featureMat[, pbulkAge == "Young"] + 0.01))
  }else if(feature == "substructure"){
    featureDiff <- rowMeans(featureMat[, pbulkAge == "Old"]) - rowMeans(featureMat[, pbulkAge == "Young"])
  }
  
  featurePvals <- twoSampTTest(featureMat[, pbulkAge == "Old"], featureMat[, pbulkAge == "Young"])
  featureFDRs <- p.adjust(featurePvals, method = "fdr")
  featureFDRs[is.na(featureFDRs)] <- 1
  diffResults[[feature]][["pvals"]] <- featurePvals
  diffResults[[feature]][["diff"]] <- featureDiff
  diffResults[[feature]][["FDRs"]] <- featureFDRs
  diffResults[[feature]][["average"]] <- featureMeans
  
}

# Match results for substructure and the corresponding CREs
overlap <- findOverlaps(regions, substructureRanges)
diffCompare <- data.frame(
  CREDiff = -log10(diffResults[["CRE"]][["FDRs"]][overlap@from]) * sign(diffResults[["CRE"]][["diff"]][overlap@from]),
  CREInd = overlap@from,
  substructureDiff = -log10(diffResults[["substructure"]][["FDRs"]][overlap@to]) * sign(diffResults[["substructure"]][["diff"]][overlap@to]),
  substructureInd = overlap@to
)

# For each CRE, keep the most differential substructure
diffCompare <- as.data.frame(diffCompare %>% group_by(CREInd) %>% filter(abs(substructureDiff) == max(abs(substructureDiff))))
write.table(diffCompare, "../../data/mHSCAging10xMultiome/diff_CRE_substructure_compare.tsv",
            sep = "\t", row.names = F, quote = F)

###################################################################
# Compare differential target gene RNA with differential CRE/TFBS #
###################################################################

# Find TSSs of differential genes
TSS <- BuenRTools::mm10TSSRanges
diffTSS <- TSS[TSS$gene_name %in% rownames(diffRNA)[diffRNA$padj < 0.1]]
diffTSS$log2FC  <- diffRNA$log2FoldChange[diffTSS$gene_name]

# Keep differential testing results around these TSS sites
diffGeneData <- diffCompare[diffCompare$CREInd %in% findOverlaps(diffTSS, regions)@to, ]

# Retrieve the corresponding gene name and signed log10(fdr)
diffGeneDataRNA <- data.table::rbindlist(pbmcapply::pbmclapply(
  diffGeneData$CREInd,
  function(CREInd){
    diffGene <- as.character(subsetByOverlaps(diffTSS, regions[CREInd])$gene_name)
    RNADiff <- -log10(diffRNA[diffGene, ]$padj) * sign(diffRNA[diffGene, ]$log2FoldChange)
    filter <- abs(RNADiff) == max(abs(RNADiff))
    data.frame(gene = diffGene[filter], RNADiff = RNADiff[filter])
  },
  mc.cores = 16
))
diffGeneData <- as.data.frame(cbind(diffGeneData, diffGeneDataRNA))

# Visualize the comparison
plotX <- diffGeneData$CREDiff
plotY <- diffGeneData$substructureDiff
plotX[!is.finite(plotX)] <- 0
plotY[!is.finite(plotY)] <- 0
density <- get_density(plotX, plotY)
plotData <- data.frame(CREDiff = plotX, substructureDiff = plotY, density = density,
                       RNADiff = diffGeneData)
pdf("../../data/mHSCAging10xMultiome/plots/diffCompareCRESubstructure.pdf",
    width = 5.5, height = 5)
ggplot(plotData) + 
  geom_point(aes(x = CREDiff, y = substructureDiff, color = density ^ 0.8), size = 0.5) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  xlab("Differential CRE signed log10(FDR)") + ylab("Differential substructure signed log10(FDR)") +
  ylim(-12, 12) + xlim(-12, 12) +
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

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Mus musculus","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")
filter <- (abs(diffGeneData$CREDiff) < 1) & (abs(diffGeneData$substructureDiff) > 1) 
enrichment <- pathwayEnrichment(
  fgGenes = diffGeneData$gene[filter], 
  bgGenes = as.character(diffTSS$gene_name), 
  geneSets = c5GO
)
as.data.frame(enrichment)
write.table(enrichment, "../../data/mHSCAging10xMultiome/diffSubstructurePathways.tsv",
            sep = "\t", quote = F)

# Visualize results
plotData <- as.data.frame(enrichment[10:1, ])
plotData$pathway <- stringr::str_replace_all(plotData$pathway, "_", " ")
plotData$pathway <- sapply(plotData$pathway, function(s){gsub('(.{1,25})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$pathway <- factor(plotData$pathway, levels = plotData$pathway) # This keeps the entries in the original order when plotting
plotData$logP <- -log10(plotData$pval)
pdf("../../data/mHSCAging10xMultiome/plots/diffSubstructurePathways.pdf",
    width = 5.5, height = 5)
ggplot(plotData) +
  geom_bar(aes(x = pathway, y = logP), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Enriched pathways")  + 
  ylab("Log10(p-value)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12))    
dev.off()

############################################
# Cluster pseudobulks and TF binding sites #
############################################

signedLogFDRs <- diffCompare$substructureDiff
selectedSites <- diffCompare$substructureInd[(abs(diffCompare$substructureDiff) > 1) & (abs(diffCompare$CREDiff) < 1)]
substructureMatClust <- substructureMat[selectedSites, ]

# Row normalize
substructureMatClust <- (substructureMatClust - rowMeans(substructureMatClust)) / rowSds(substructureMatClust)

# Categorize substructures
ageClusters <- unique(ageCluster)
substructureCluster <- pbmcapply::pbmcmapply(
  function(i){
    dt <- substructureMatClust[i, ]
    clusterDt <- sapply(ageClusters, function(x){mean(dt[ageCluster == x])})
    cluster <- names(clusterDt[clusterDt == max(clusterDt)])[1]
    if(cluster == "Old"){cluster <- "I"
    }else if(cluster == "Young"){cluster <- "IV"
    }else if(cluster == "Young-like old"){
      if(clusterDt["Old"] > clusterDt["Young"]){cluster <- "II"}
      else{cluster <- "III"}
    }
    cluster
  },
  1:dim(substructureMatClust)[1],
  mc.cores = 16
)

# Visualize clustering results
colors <- colorRamp2(seq(0.2, 1,length.out=9),
                     colors = jdb_palette("brewer_blue"))
pdf("../../data/mHSCAging10xMultiome/plots/diffSubstructureClusterMap.pdf",
    width = 6, height = 5)
Heatmap(matrix = t(substructureMatClust),
        col = colors,
        name = "Scaled\nsubstructure\nscores",
        row_split = factor(ageCluster, levels = c("Old", "Young-like old", "Young")),
        column_split = factor(substructureCluster, levels = c("I", "II", "III", "IV")),
        cluster_rows = F, cluster_columns = F)
dev.off()
