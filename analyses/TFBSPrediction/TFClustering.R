# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getCounts.R")
source("../../code/getBias.R")
source("../../code/getFootprints.R")
source("../../code/visualization.R")
source("../../code/getAggregateFootprint.R")

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)
library(ggpubr)

####################################################
# Retrieve multi-scale aggregate footprints of TFs #
####################################################

# Load multi-scale footprints from each dataset
datasets <- c("HepG2", "K562", "GM12878")
TFMultiScaleFootprints <- lapply(
  datasets,
  function(dataset){
    readRDS(paste0("../../data/", dataset, "/TFMultiScale.rds"))
  }
)
names(TFMultiScaleFootprints) <- datasets

# Then select representative kernel sizes and compare footprints in a local neighborhood
kernelSizes <- 2:100
selectedSizes <- c(10, 20, 30, 50, 80, 100)
radius <- 300
TFSelectedScaleFootprints <- Reduce(rbind, lapply(
  datasets,
  function(dataset){
    datasetFootprints <- t(sapply(
      TFMultiScaleFootprints[[dataset]],
      function(x){
        colnames(x) <- kernelSizes
        Reduce(c, lapply(
          selectedSizes,
          function(footprintRadius){
            x[(radius - 99):(radius + 100), as.character(footprintRadius)]
          }
        ))
      }
    ))
    rownames(datasetFootprints) <- paste0(names(TFMultiScaleFootprints[[dataset]]), "_", dataset)
    datasetFootprints
  }
))
TFNames <- Reduce(c, lapply(datasets, function(dataset){names(TFMultiScaleFootprints[[dataset]])}))
uniqueTFs <- unique(TFNames)

# If a TF appears in more than one cell line, keep the one with the highest signal
TFAvgFootprints <- t(
  sapply(
    uniqueTFs,
    function(TF){
      maxSignal <- rowMaxs(TFSelectedScaleFootprints[TFNames %in% TF, , drop = F])
      keptSample <- which(maxSignal == max(maxSignal))
      TFSelectedScaleFootprints[TFNames %in% TF, , drop = F][keptSample,]
    }
  )
)
rownames(TFAvgFootprints) <- uniqueTFs

###############################
# Assign TF family membership #
###############################

# Manually curated using information here: http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/family
TFFamilyMapping <- list(
  "Basic domains" = c("AP-2", "bHLH", "Nrf1", "RFX", "TF_bZIP", "TSC22"),
  "Beta-scaffold" = c("CP2", "STAT", "CBF", "CSD", "CSL", "GCM", "P53", "RHD",
                      "Runt"),
  "Helix-turn-helix" = c("ARID", "COE", "CUT", "E2F", "ETS", "Fork_head", 
                         "Homeobox", "HPD", "HSF", "HTH", "IRF", "MYB", 
                         "SRF", "TEA", "PAX", "Pou", "TF_Otx"),
  "Other alpha-helix" = c("HMG", "NF-YA","NF-YB", "NF-YC", "SAND", "GTF2I"),
  "Zinc-Coordinating" = c("DM","ESR-like","GCNF-like","Miscellaneous",
                          "NGFIB-like","RXR-like","SF-like","THAP","THR-like","ZBTB","zf-BED","zf-C2H2",
                          "zf-C2HC","zf-CCCH","zf-GAGA","zf-GATA","zf-LITAF-like","zf-MIZ","zf-NF-X1"),
  "Unclassified" = c("AF-4","CG-1","CSRNP","CTF/NFI",
                     "DACH","GCFC","HMGA","LRRFIP","MBD","MH1",
                     'NCU-G1',"NDT80/PhoG","PC4","T-box","Tub","Others")
)

# Map smaller families to larger families
TFFamilyMapping <- Reduce(c, lapply(
  names(TFFamilyMapping),
  function(Family){
    x <- rep(Family, length(TFFamilyMapping[[Family]]))
    names(x) <- TFFamilyMapping[[Family]]
    x
  }
))

# Retrieve TF family annotation
TFFamilies <- read.csv("../../data/shared/TFTypes.txt", sep = "\t")
selectedTFs <- intersect(uniqueTFs, unique(TFFamilies$Symbol))
TFAvgFootprints <- TFAvgFootprints[selectedTFs, ]
TFFamilies <- TFFamilies[match(selectedTFs, TFFamilies$Symbol), ]$Family
TFFamilies <- unname(TFFamilyMapping[TFFamilies])

###############################################
# Cluster TFs based on multi-scale footprints #
###############################################

# Cluster and order TFs based on multi-scale footprints
TFClust <- hclust(distCosine(TFAvgFootprints))
nClusters <- 6
clusterLabels <- cutree(TFClust, k = nClusters)
TFOrder <- TFClust$order

#####################
# Visualize results #
#####################

clusterColors <- as.character(MetBrewer::met.brewer("Derain", nClusters))
names(clusterColors) <- 1:nClusters
leftAnno <- rowAnnotation(Cluster=factor(clusterLabels[TFOrder]),
                          border=TRUE,
                          col=list(Cluster=clusterColors),
                          show_annotation_name=FALSE)

familyColors <- get_palette("npg", length(unique(TFFamilies)))
names(familyColors) <- unique(TFFamilies)
rightAnno <- rowAnnotation(Family=factor(TFFamilies[TFOrder]),
                           col=list(Family=familyColors),
                           gap = unit(2, "points"),
                           show_annotation_name=FALSE)

# Split heatmap columns by kernel size
colGroups <- Reduce(c, lapply(
  selectedSizes,
  function(r){
    paste(rep(r, 200), "bp")
  }
))
colors <- colorRamp2(seq(0, 1, length.out=9),
                     colors = colorRampPalette(c(rep("white", 2),  
                                                 "#9ECAE1", "#08519C", "#08306B"))(9))
system("mkdir -p ../../data/TFBSPrediction/plots")
pdf("../../data/TFBSPrediction/plots/TFClusters.pdf", width = 20, height = 22)
Heatmap(TFAvgFootprints[TFOrder,],
        left_annotation = leftAnno,
        right_annotation = rightAnno,
        column_split = factor(colGroups, levels = paste(selectedSizes, "bp")),
        col = colors,
        border = TRUE,
        name = "Footprint score",
        cluster_columns = F,
        cluster_rows = F)
dev.off()

# Save clustering results to a file
clusterDf <- data.frame(
  TF = uniqueTFs,
  cluster = clusterLabels
)
clusterDf <- clusterDf[order(clusterDf$cluster),]
write.table(clusterDf, "../../data/TFBSPrediction/clusterLabels.txt",
            row.names = F, quote = F, sep = "\t")

#################################################################################
# Also assign cluster membership to TFs not included in the above ChIP datasets #
#################################################################################

# Load PWM data
cisBPMotifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))

# Using motif PWM similarity, match each TF to a TF with pre-existing cluster label
matchedTFs <- pbmcapply::pbmcmapply(
  function(TF){
    similarity <- TFBSTools::PWMSimilarity(cisBPMotifs[clusterDf$TF], cisBPMotifs[[TF]], method = "Pearson")
    matchedTF <- names(sort(similarity, decreasing = T)[1])
  },
  names(cisBPMotifs),
  mc.cores = 16
)

# Assign cluster membership to each TF
motifClustering <- data.frame(
  TF = names(matchedTFs),
  cluster = clusterDf$cluster[match(matchedTFs, clusterDf$TF)]
)

# Re-order results by cluster
motifClustering <- motifClustering[order(motifClustering$cluster),]

write.table(motifClustering, "../../data/TFBSPrediction/clusterLabelsAllTFs.txt",
            row.names = F, quote = F, sep = "\t")

#####################################################
# Quantify number of TFs in each family and cluster #
#####################################################

uniqueFamilies <- unique(TFFamilies)
plotData <- data.table::rbindlist(lapply(
  1:nClusters,
  function(clusterInd){
    freq <- table(factor(TFFamilies[clusterLabels == clusterInd], levels = uniqueFamilies))
    data.frame(
      Cluster = as.character(clusterInd),
      Family = uniqueFamilies,
      Frequency = as.integer(unname(freq))
    )
  }
))

pdf("../../data/TFBSPrediction/plots/TFClusterBarplot.pdf", width = 6, height = 6)
ggpubr::ggbarplot(plotData, x = "Cluster", y = "Frequency",
                  fill = "Family", color = "Family", width = 0.5,
                  palette = get_palette("npg", 7)) +
  labs(x = "TF Cluster", y = "Number of TFs") 
dev.off()
