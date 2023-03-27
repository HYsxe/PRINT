# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
library(ggplot2)
library(ggrepel)

# Load input data
dataset <- "GM12878"
TFFootprintContribution <- read.table(paste0("../../data/TFBSPrediction/", dataset, "_TF_fp_contribution.tsv"))
rownames(TFFootprintContribution) <- TFFootprintContribution$TF
performance <- read.table(paste0("../../data/TFBSPrediction/", dataset, "Performance.txt"))
rownames(performance) <- performance$TF

# Only keep overlapping TFs
keptTFs <- intersect(performance$TF, TFFootprintContribution$TF)
performance <- performance[keptTFs,]
TFFootprintContribution <- TFFootprintContribution[keptTFs,]

# Compare 20 bp footprint contribution with model performnace
plotData <- data.frame(
  TF = keptTFs,
  precision = performance$precision,
  contribution = TFFootprintContribution$contribution
)

# Label top TFs
labels <- character(dim(plotData)[1])
labelInds <- (plotData$precision > 0.85)
labels[labelInds] <- keptTFs[labelInds]

# Visualize results
pdf(paste0("../../data/TFBSPrediction/plots/", dataset, "_TF_fp_contribution_vs_precision.pdf"), 
    height = 5, width = 6)
ggplot(data = plotData) +
  geom_point(aes(x = contribution, y = precision)) +
  geom_text_repel(x = plotData$contribution, 
                  y = plotData$precision, 
                  label = labels, size = 20) +
  xlab("TF footprint contribution score\n(40 bp scale)") +
  ylab("Model precision") +
  ggtitle(dataset) +
  theme_classic()
dev.off()