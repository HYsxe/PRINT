# If running in Rstudio, set the working directory to current path
if (Sys.getenv("RSTUDIO") == "1"){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(ggplot2)
library(Matrix)

# Randomly sample local neighborhoods in the cell state space and group them into pseudobulks
getPseudobulks <- function(cellEmbedding, # Cell-by-feature dimensionality reduction matrix (e.g., cell-by-LSI, cell-by-PC)
                           nGroups = 1000, # Number of pseudobulks to generate
                           nScaffoldCells = 10000, # Number of scaffold cells to sample for calculation of pseudotime
                           pseudobulkSize = 5000, # Number of cells per pseudobulk
                           nNbrs = 10 # To maximize coverage of cell states, we construct a KNN graph and try to sample
                           # the cells with lower in-degrees in the graph as pseudobulk centers. This is the k used in the
                           # construction of the KNN graph
                           ){
  
  set.seed(42)
  
  nCells <-  dim(cellEmbedding)[1]
  
  # Randomly select scaffold cells
  scaffoldInds <- sample(1:nCells, nScaffoldCells, replace = F)
  scaffoldEmbedding <- cellEmbedding[scaffoldInds,]
  
  # Construct KNN graph among the scaffold cells
  cat("Constructing KNN graph among scaffold cells\n")
  scaffoldKNN <- FNN::get.knn(scaffoldEmbedding, k = nNbrs)$nn.index
  KNNGraphIJV <- lapply(1:dim(scaffoldKNN)[1], function(i){
    as.data.frame(cbind(rep(i, nNbrs),
                        scaffoldKNN[i,],
                        rep(1, nNbrs)))
  })
  KNNGraphIJV <- as.data.frame(data.table::rbindlist(KNNGraphIJV))
  KNNGraph <- sparseMatrix(i = KNNGraphIJV[, 1],
                           j = KNNGraphIJV[, 2],
                           x = KNNGraphIJV[, 3],
                           dims = c(nScaffoldCells, nScaffoldCells))
  
  # For every cell, count its in-degree (i.e, it is in the KNN of how many cells) in the KNN graph.
  inDegrees <- table(factor(c(scaffoldKNN), levels = 1:nScaffoldCells))
  
  # Select pseudobulk center cells 
  cat("Designating scaffold cells with lowest in-degrees as pseudobulk centers\n")
  centerInds <- scaffoldInds[as.integer(names(sort(inDegrees))[1:nGroups])]
  
  # For each pseudobulk center, find its k nearest neighbor cells to form the pseudobulk
  cat("Finding member cells for each pseudobulk center\n")
  groupMembers <- FNN::get.knnx(cellEmbedding, 
                                cellEmbedding[centerInds,],
                                k = pseudobulkSize)$nn.index
  
  list("centerInds" = centerInds, 
       "groupMembers" = groupMembers)
}
