# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/getPseudoBulks.R")
library(SummarizedExperiment)
library(hdf5r)
library(reticulate)
library(ggplot2)
library(cisTopic)

#####################
# 1.Load input data #
#####################

# Load single cell ATAC
scATAC <- readRDS("../../data/BMMC/atac.se.rds")

# Load cisTopic results
cisTopicObject <- readRDS("../../data/BMMC/cisTopicObject.rds")
cellEmbedding <- modelMatSelection(cisTopicObject, 'cell', 'Z-score')
cellEmbedding <- t(cellEmbedding)
cellEmbedding <- cellEmbedding[colnames(scATAC), ]
saveRDS(cellEmbedding, "../../data/BMMC/cellEmbedding.rds")
rm(cisTopicObject)

# Get region ranges
regions <- rowRanges(scATAC)
regions <- IRanges::resize(regions, width = 1000, fix = "center")
saveRDS(regions, "../../data/BMMC/regionRanges.rds")

###############################################
# 2. Only keep motifs enriched in CRE regions #
###############################################

# Load PWM data
motifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))
filtMotifs <- filterMotifs(regions = regions, 
                           motifs = motifs,
                           genome = "hg38")
saveRDS(motifs, "../../data/BMMC/filteredMotifs.rds")

###################################
# 3.Test pseudobulking parameters #
###################################

# Get cell type labels per pseudobulk
cellTypeLabels <- as.character(scATAC@colData$cistopic.assign.l2.rank)

# Load PWM data
motifs <- readRDS(paste0(projectMainDir, "/data/shared/cisBP_human_pwms_2021.rds"))

coverageDf <- lapply( 
  seq(100, 1000, 100),
  function(nGroups){
    
    # Pseudobulk cells with selected nGroups
    print(paste0("Testing nGroups = ", nGroups))
    pseudobulking <- getPseudobulks(cellEmbedding, nGroups = nGroups)
    groupMembers <- pseudobulking$groupMembers
    centerInds <- pseudobulking$centerInds
    
    # Calculate coverage for each cell type
    coveredInds <- unique(c(groupMembers))
    coveredCellTypes <- cellTypeLabels[coveredInds]
    cellTypeCoverage <- sapply(
      sort(unique(cellTypeLabels)),
      function(cellType){
        sum(coveredCellTypes %in% cellType) / sum(cellTypeLabels %in% cellType)
      }
    )
    cellTypeCoverage[["All"]] <- length(coveredInds) / length(cellTypeLabels)
    
    as.data.frame(cbind(names(cellTypeCoverage), 
                        cellTypeCoverage,
                        nGroups)) 
  }
)
coverageDf <- data.table::rbindlist(coverageDf)
colnames(coverageDf) <- c("Cell type", "Coverage", "Pseudobulk number")
coverageDf$Coverage <- as.numeric(coverageDf$Coverage)
coverageDf$`Pseudobulk number` <- as.integer(coverageDf$`Pseudobulk number`)

# Assign a color to each cell type
cellTypeColors <- c("#0E4421","#00AF99","#9AD9E9","#A896C8","#FAA31B","#F28238",
                    "#F26625","#D13C27","#46A147","#EF3741","#CC1F39","#901838",
                    "#A896C8","#B280B9","#C757A1","#854199","#005F9F","#1481C4","#492264",
                    "#67B545", "Black")
names(cellTypeColors) <- c("HSC/MPP", "LMPP","CLP", "pro/pre-B","GMP", "CD14mono", "CD16mono", "DC",
                           "CMP", "MEP", "early-Ery", "late-Ery", "NaiveB", "MemoryB", "plasmaB", 
                           "pDC", "CD4", "CD8", "NK", "Baso", "All")

if(!dir.exists("../../data/BMMC/plots/preprocessing")){
  system("mkdir -p ../../data/BMMC/plots/preprocessing")
}
pdf("../../data/BMMC/plots/preprocessing/nPseudobulks.pdf",height = 6,width = 7)
ggpubr::ggline(coverageDf, 
               x = "Pseudobulk number", 
               y = "Coverage",
               color = "Cell type") +
  scale_colour_manual(values = cellTypeColors)
dev.off()

##########################
# 4.Pseudobulk the cells #
##########################

# Pseudobulk the cells. Get the barcodes of pseudobulk members and the index of pseudobulk center cell.
nGroups <- 1000
pseudobulking <- getPseudobulks(cellEmbedding, nGroups = nGroups)
groupMembers <- pseudobulking$groupMembers
centerInds <- pseudobulking$centerInds

# Get cell type labels per pseudobulk
cellTypeLabels <- as.character(scATAC@colData$cistopic.assign.l2.rank)
groupCellTypes <- cellTypeLabels[centerInds]

# Check depth of pseudobulks
pseudobulkDepths <- sapply(
  1:nGroups,
  function(i){
    sum(scATAC$depth[groupMembers[i,]])
  }
)

# Check average pseudobulk depth per cell type
cellTypeDepths <- sapply(
  unique(cellTypeLabels),
  function(celltype){
    mean(pseudobulkDepths[groupCellTypes %in% celltype])
  }
)
names(cellTypeDepths) <- unique(cellTypeLabels) 
ggpubr::gghistogram(pseudobulkDepths) + 
  xlab("Log10(depth)") +
  ylab("Number of pseudobulks")

# Visualize pseudobulk members
pbulkInd <- sample(which(cellTypeLabels[centerInds] == "early-Ery"), 1)
plotDf <- data.frame(
  UMAP1 = scATAC@colData$cisTopic.umap1,
  UMAP2 = scATAC@colData$cisTopic.umap2,
  celltype = cellTypeLabels
)
ggplot(plotDf) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = celltype), size = 0.1) +
  geom_point(data = plotDf[groupMembers[pbulkInd,],],
             aes(x = UMAP1, y = UMAP2), color = "black") +
  ggtitle(paste0("Pseudobulk No.", pbulkInd, 
                 " cell type: ", cellTypeLabels[centerInds][pbulkInd])) +
  theme_classic()

# Reformat cell barcodes
stopifnot(all(rownames(cellEmbedding) == colnames(scATAC)))
scBarcodes <- rownames(cellEmbedding)

# Save barcodes of each pseudobulk
barcodeGroups <- pbmcapply::pbmclapply(
  1:nGroups,
  function(groupInd){
    data.frame(barcode = scBarcodes[groupMembers[groupInd,]],
               group = groupInd)
  }
)
barcodeGroups <- data.table::rbindlist(barcodeGroups)
write.table(barcodeGroups, "../../data/BMMC/barcodeGrouping.txt",
            row.names = F, sep = "\t", quote = F)

###################################
# 5.Get pseudotime of pseudobulks #
###################################

set.seed(42)
nScaffold <- 100000
palantirPath <-"../../data/BMMC/palantir.h5"
if(file.exists(palantirPath)){
  system(paste0("rm ", palantirPath))
}

# Use scaffold cells to increase resolution of pseudotime
scaffoldInds <- sample(setdiff(1:dim(cellEmbedding)[1], centerInds), nScaffold)

# Write input data of Palantir to a file
h5File <- H5File$new(palantirPath, mode="w")
h5File[["cellCisTopic"]] <- cellEmbedding[c(centerInds, scaffoldInds),]
h5File[["barcodes"]] <- scBarcodes[c(centerInds, scaffoldInds)]
h5File[["rootCellBarcode"]] <- "T128.R1.125.R2.039.R3.149.P1.04" # Specify a root cell
h5File$close_all()

# Pass argumet to python script
write.table("../../data/BMMC/palantir.h5", 
            "palantir_args.txt",
            quote = F, col.names = F, row.names = F)
py_run_file("../../code/Palantir.py")
system("rm palantir_args.txt")

# Read in the pseudotime computed using Palantir
h5File <- H5File$new(palantirPath, mode="r")
pseudotime <- h5File[["pseudotime"]]
pseudotime <- pseudotime[1:pseudotime$dims]
pseudotime <- pseudotime[1:nGroups]

# Read lineage probability predicted by Palantir
lineageProbs <- h5File[["lineageProbs"]]
lineageProbs <- lineageProbs[1:lineageProbs$dims[1],
                             1:lineageProbs$dims[2]]
lineageProbs <- lineageProbs[, 1:nGroups]
lineageProbs <- t(lineageProbs)
lineageIdentities <- sapply(1:dim(lineageProbs)[1], 
                            function(x){which(lineageProbs[x,] == max(lineageProbs[x,]))})

###################################
# 6. Visualize pseudotime results #
###################################

# Some specifications for plotting
mytheme2 <- theme(text = element_text(size=24*0.2),
                  axis.text = element_text(size=20*0.2), 
                  axis.line = element_line(colour = 'black', size = 1*0.2),
                  axis.ticks = element_line(colour = 'black', size = 1*0.2),
                  axis.ticks.length=unit(.2*0.2, "cm"),
                  legend.text=element_text(size=16*0.2),
                  legend.title=element_text(size=16*0.2),
                  legend.key.height = unit(.15, "cm"),
                  legend.key.width = unit(.15, "cm"))

# Get pseudotime for each single cell by smoothing
smoothKNN <- FNN::get.knnx(data = cellEmbedding[centerInds,],
                           query = cellEmbedding)$nn.index
smoothedPseudoTime <- pbmcapply::pbmcmapply(
  function(cellInd){
    mean(pseudotime[smoothKNN[cellInd, ]])
  },
  1:dim(smoothKNN)[1],
  mc.cores = 16
)

# Visualize pseudotime on the UMAP
pltDf <- data.frame(UMAP1 = scATAC$cisTopic.umap1,
                    UMAP2 = scATAC$cisTopic.umap2,
                    Pseudotime = smoothedPseudoTime)
ptimePlot <- ggplot(pltDf) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = Pseudotime), 
             size = 0.02, alpha = 0.5, stroke = 0) +
  scale_color_gradientn(colors = BuenColors::jdb_palette("solar_extra")) +
  BuenColors::pretty_plot() + BuenColors::L_border() +
  mytheme2
ggsave("../../data/BMMC/plots/preprocessing/pseudoTimeErythroid.png", 
       plot=ptimePlot, width=4.5, height=3, units = "cm", dpi=1200)

# Visualize lineage probabilities
lineageInd <- 1
smoothedLineageProbs <- pbmcapply::pbmcmapply(
  function(cellInd){
    mean(lineageProbs[smoothKNN[cellInd, ], lineageInd])
  },
  1:dim(smoothKNN)[1],
  mc.cores = 16
)

# Visualize lineage probabilities
pltDf <- data.frame(UMAP1 = scATAC$cisTopic.umap1,
                    UMAP2 = scATAC$cisTopic.umap2,
                    lineageP = smoothedLineageProbs)
linPlot <- ggplot(pltDf) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = lineageP), 
             size = 0.02, stroke = 0) +
  BuenColors::pretty_plot() + BuenColors::L_border() +
  mytheme2
ggsave(paste0("../../data/BMMC/plots/preprocessing/lineage",
              lineageInd, "Prob.png"), 
       plot=linPlot, width=4.5, height=3, units = "cm", dpi=1200)

###########################
# 7. Save results to file #
###########################

# Annotate the lineages
colnames(lineageProbs) <- c("MyeloidProbs", "ErythroidProbs", "BLymphoidProbs", "TLymphoidProbs")

# Save metadata of pseudobulks
groupInfo <- data.frame(
  groupID = 1:nGroups,
  cellType = cellTypeLabels[centerInds],
  UMAP1 = scATAC$cisTopic.umap1[centerInds],
  UMAP2 = scATAC$cisTopic.umap2[centerInds],
  Pseudotime = pseudotime,
  color = cellTypeColors[as.character(scATAC@colData$cistopic.assign.l2.rank)[centerInds]]
)
groupInfo <- cbind(groupInfo, lineageProbs)

write.table(groupInfo, "../../data/BMMC/groupInfo.txt",
            row.names = F, sep = "\t", quote = F)
