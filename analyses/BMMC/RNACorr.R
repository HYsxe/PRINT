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
source("../../code/getTFBS.R")
source("../../code/getGroupData.R")
source("../../code/getSubstructures.R")
source("../../code/getGeneCorr.R")

set.seed(42)
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

###################
# Load input data #
###################

# Initialize a footprintingProject object
projectName <- "BMMC"
projectMainDir <- "../../"
projectDataDir <- paste0(projectMainDir, "data/", projectName, "/")

# Load the footprintingProject object
project <- readRDS(paste0(projectDataDir, "footprintingProject.rds"))

# Get color mapping for cell types
groupInfo <- read.table(paste0(projectDataDir, "groupInfo.txt"), header = T,
                        comment.char = "")
cellTypes <- unique(groupInfo$cellType)
cellTypeColors <- groupInfo$color[match(cellTypes, groupInfo$cellType)]
names(cellTypeColors) <- cellTypes

# Get substructure-by-pseudobulk SummarizedExperiment object
substructurePath <- "../../data/BMMC/substructureSE.rds"
if(file.exists(substructurePath)){
  substructureSE <- readRDS("../../data/BMMC/substructureSE.rds")
}else{
  substructureSE <- getSubstructureSE(project)
  saveRDS(substructureSE, substructurePath)
}

# Filter low signal substructures
substructureSignal <- rowMaxs(substructureSE@assays@data$substructures)
substructureFilter <- which(substructureSignal > 0.3)
substructureSE <- substructureSE[substructureFilter, ]

# Get CRE-by-pseudobulk SummarizedExperiment object
CRESE <- SummarizedExperiment(assays = list("CREs" = groupATAC(project)),
                              rowRanges = regionRanges(project))

##################################################################
# Calculate substructure-RNA correlation and CRE-RNA correlation #
##################################################################

# Compute substructure-gene correlation
substructureGeneCorr <- geneCorr(ATAC = substructureSE, 
                                 RNA = groupRNA(project),
                                 genome = refGenome(project),
                                 nCores = 16)
colnames(substructureGeneCorr) <- c("Gene", "substructureInd", "rObs", "pvalZ")
saveRDS(substructureGeneCorr, "../../data/BMMC/substructureGeneCorr.rds")

# Compute CRE-gene correlation
CREGeneCorr <- geneCorr(ATAC = CRESE, 
                        RNA = groupRNA(project),
                        genome = refGenome(project),
                        nCores = 16)
colnames(CREGeneCorr) <- c("Gene", "CREInd", "rObs", "pvalZ")
saveRDS(CREGeneCorr, "../../data/BMMC/CREGeneCorr.rds")

############################################################
# Comapre substructure-gene correlation with CRE-gene correlation #
############################################################

substructureGeneCorr <- readRDS("../../data/BMMC/substructureGeneCorr.rds")
CREGeneCorr <- readRDS("../../data/BMMC/CREGeneCorr.rds")

ovList <- findOverlaps(rowRanges(substructureSE), rowRanges(CRESE))

rownames(CREGeneCorr) <- paste0(CREGeneCorr$Gene,
                                "_", CREGeneCorr$CREInd)

selectedCREs <- 1:length(regionRanges(project))
selectedGene <- NULL
corCompare <- pbmcapply::pbmclapply(
  selectedCREs,
  function(CREInd){
    
    # Find all substructures mapped to the current CRE
    substructureInds <- ovList@from[ovList@to == CREInd]
    
    if(sum(substructureGeneCorr$substructureInd %in% substructureInds) > 0){
      
      # Retreive footprint-gene correlation for the selected footprints
      if(!is.null(selectedGene)){
        substructureCorrs <- substructureGeneCorr[(substructureGeneCorr$substructureInd %in% substructureInds) &
                                                    (substructureGeneCorr$Gene == selectedGene),]
      } else{
        substructureCorrs <- substructureGeneCorr[(substructureGeneCorr$substructureInd %in% substructureInds),]
      }
      
      # Retrieve the corresponding region-gene correlation
      CRECorrs <- CREGeneCorr[paste(substructureCorrs$Gene, CREInd, sep = "_"),]
      
      # Combine the results
      comparison <- data.frame(CREObsR = CRECorrs$rObs, 
                               CREPvalZ = CRECorrs$pvalZ, 
                               substructureObsR = substructureCorrs$rObs, 
                               substructurePvalZ = substructureCorrs$pvalZ, 
                               CREInd = CREInd, 
                               substructureInd = substructureCorrs$substructureInd)
      comparison$Gene <- substructureCorrs$Gene
      comparison <- comparison %>% group_by(Gene) %>% filter(abs(substructurePvalZ) == min(abs(substructurePvalZ)))
      
      comparison
    }
  },
  mc.cores = 16
)
corCompare <- corCompare[sapply(corCompare, function(x){!is.null(x)})]
corCompare <- data.table::rbindlist(corCompare)
corCompare <- as.data.frame(corCompare)
corCompare[, 1:6] <- apply(corCompare[, 1:6], 2,  as.numeric)
corCompare <- corCompare[!is.na(rowSums(corCompare[,1:6])),]
colnames(corCompare) <- c("CREGeneCorr", "CREGenePval",
                          "substructureGeneCorr", "substructureGenePval",
                          "CREInd", "substructureInd", "Gene")

# Multiple testing correction
CREGeneFDR <- p.adjust(corCompare$CREGenePval, method = "fdr")
substructureGeneFDR <- p.adjust(corCompare$substructureGenePval, method = "fdr")
CRESignedLogFDR <- -log10(CREGeneFDR) * sign(corCompare$CREGeneCorr)
substructureSignedLogFDR <- -log10(substructureGeneFDR) * sign(corCompare$substructureGeneCorr)

# Get density for plotting
plotX <- CRESignedLogFDR
plotY <- substructureSignedLogFDR
plotX[!is.finite(plotX)] <- 0
plotY[!is.finite(plotY)] <- 0
density <- get_density(plotX, plotY)

# Visualize results
plotData <- data.frame(
  CRESignedLogFDR = CRESignedLogFDR,
  substructureSignedLogFDR = substructureSignedLogFDR,
  Significant = (abs(CRESignedLogFDR) < 1) & (abs(substructureSignedLogFDR) > 1),
  density = density ^ 0.2
) 
png("../../data/BMMC/plots/RNACorr/RNACorr.png",
    width = 5.5, height = 5)
ggplot(plotData) +
  geom_point(aes(x = CRESignedLogFDR, y = substructureSignedLogFDR, color = density),
             size = 0.5) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  geom_vline(xintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_vline(xintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = 1, linetype ="dashed", 
             color = "black", size = 0.5) +
  geom_hline(yintercept = -1, linetype ="dashed", 
             color = "black", size = 0.5) +
  xlab("CRE-gene association signed FDR") +
  ylab("Substructure-gene association signed FDR") +
  theme_classic()
dev.off()

###############################
# Pathway enrichment analysis #
###############################

# Load pathway gene sets from hypeR
library(hypeR)
c5GO <- msigdb_gsets(species = "Homo sapiens","C5","BP",clean = TRUE)$genesets
names(c5GO) <- stringr::str_replace_all(names(c5GO), " ",  "_")

# Get foreground and background genes
substructureSpecificGenes <- unique(corCompare$Gene[(abs(CRESignedLogFDR) < 1) & (abs(substructureSignedLogFDR) > 1)])
bgGenes <- unique(corCompare$Gene[(abs(CRESignedLogFDR) > 1) & (abs(substructureSignedLogFDR) < 1)])
bgGenes <- unique(c(bgGenes, substructureSpecificGenes))

# Calculate pathway enrichment
enrichment <- pathwayEnrichment(substructureSpecificGenes, bgGenes, c5GO)

# Visualize results
plotData <- as.data.frame(enrichment[10:1, ])
plotData$pathway <- stringr::str_replace_all(plotData$pathway, "_", " ")
plotData$pathway <- sapply(plotData$pathway, function(s){gsub('(.{1,25})(\\s|$)', '\\1\n', s)}) # Add newline to long strings
plotData$pathway <- factor(plotData$pathway, levels = plotData$pathway) # This keeps the entries in the original order when plotting
plotData$logP <- -log10(plotData$pval)
pdf("../../data/BMMC/plots/RNACorr/pathwayEnrichment.pdf",
    width = 5.5, height = 5)
ggplot(plotData) +
  geom_bar(aes(x = pathway, y = logP), stat = "identity", width = 0.5, fill = "#CA9B80") +
  xlab("Enriched pathways")  + 
  ylab("Log10(p-value)") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 12))  
dev.off()

# Save results to a file
write.table(enrichment, "../../data/BMMC/substructureSpecificGenePathwayEnrichment.tsv",
            quote = F, sep = "\t", col.names = T, row.names = F)
