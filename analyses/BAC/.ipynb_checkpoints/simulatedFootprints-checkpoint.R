# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source("../../code/utils.R")
source("../../code/getFootprints.R")
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(RColorBrewer)

#########################################################
# Test multi-scale footprinting on simulated footprints #
#########################################################

# Load dispersion models
dispModels <- list()
for(kernelSize in 2:100){
  dispModels[[as.character(kernelSize)]] <- readRDS(paste0("../../data/shared/dispModel/dispersionModel", kernelSize, "bp.rds"))
}

# Within a 1kb window, simulate footprints of different sizes, depths and shapes
Tn5Insertion <- rnorm(n = 1000, sd = 2, mean = 10)
scaling <- c(0.75, 0.9, 0.9, 0.9, 0.9, 0.9)
position <- c(150, 300, 450, 490, 650, 850)
width <- c(5, 10, 10, 10, 20, 40)
for(i in 1:length(scaling)){
  footprint <- dnorm(1:1000, mean = position[i], sd = width[i])
  footprint <- 1 - footprint / max(footprint) * scaling
  Tn5Insertion <- Tn5Insertion * footprint
}

# Plot Simulated Tn5 insertion
positions <- 1:1000
barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                      y1 = 0, y2 = Tn5Insertion)
ggplot(barData) + 
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
            fill = "#006838") +
  scale_y_continuous(expand = c(0, 0)) + xlab("") + ylab("Simulated Tn5 insertion") +
  theme_classic() 

# Calculate multi-scale footprints using the simulated data, assuming Tn5 bias = 1 for all positions
multiScaleFootprints <- pbmcapply::pbmcmapply(
  function(footprintRadius){
    footprintPvals <- footprintScoring(
      Tn5Insertion = Tn5Insertion,
      Tn5Bias = rep(1, length(Tn5Insertion)),
      dispersionModel = dispModels[[as.character(footprintRadius)]],
      footprintRadius = footprintRadius,
      flankRadius = footprintRadius
    )
    # Remove regions affected by edge effect
    # For edges, the sum of bias and counts on the left / right flank might
    # be much lower than what the model has seen, since we're adding many zeros.
    # As a result, the model prediction of ratioSD will not be accurate
    footprintPvals[is.na(footprintPvals)] <- 1
    footprintScores <- -log10(footprintPvals)
  },
  100:4,
  mc.cores = 16
)

# Visualize multi-scale footprints
footprintColors <- colorRamp2(seq(-log10(0.05),3,length.out=9),
                              colors = jdb_palette("brewer_blue"))
Heatmap(t(multiScaleFootprints), 
        cluster_rows = F, 
        cluster_columns = F,
        col = footprintColors)
