# Get position-by-pseudobulk matrix of the feature of interest
getFeatureMatrix <- function(project, # footprintingProject object
                             regionInd, # Integer. Index of the region to plot
                             feature, # Character. Feature to plot. Can be one of "footprint", "insertion", and "TFBS"
                             groupIDs, # IDs of pseudobulks to plot
                             footprintRadius = NULL, # Radius of footprint. 
                             smoothRadius = NULL # Radius used for smoothing plot
                             ){
  
  # Retrieve data for plotting
  if(feature == "footprint"){
    
    if(is.null(footprintRadius)){
      stop("When plotting footprints as features, footprintRadius must be provided!")
    }
    
    # Get position-by-pseudobulk matrix of footprint scores
    featureMatrix <- regionFootprintMatrix(
      project = project, 
      regionInd = regionInd, 
      footprintRadius = footprintRadius,
      lineageGroups = groupIDs)
    
  }else if (feature == "insertion"){
    
    ret <- getCountData(project, regionInd)
    countData <- ret[["countData"]]
    adjustedRegionInd <- ret[["regionInd"]]
    
    # Get position-by-pseudobulk ATAC insertion matrix
    regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, regionWidth(project)[regionInd])
    
    # Smoothe data across genomic positions
    featureMatrix <- sapply(
      1:dim(regionATAC)[2],
      function(i){
        conv(regionATAC[,i], smoothRadius) / (2 * smoothRadius)
      }
    )
    
  }else if (feature == "TFBS"){
    
    if(length(TFBindingModel(project)) == 0){
      stop("No TFBS model loaded!")
    }
    
    # Get position-by-pseudobulk matrix of TF binding scores
    featureMatrix <- getRegionTFBSMatrix(project = project, 
                                       regionInd = regionInd, 
                                       groupIDs = groupIDs)
  }
  
  featureMatrix
}

# Plot position-by-pseudobulk heatmap of a selected feature for a specific CRE region
setGeneric("plotFeatureHeatmap", function(project, # footprintingProject object
                                          regionInd, # Integer. Index of the region to plot
                                          feature = "footprint",  # Character. Feature to plot. Can be one of "footprint", "insertion", and "TFBS"
                                          gene = NULL, # Whether to also plot the RNA levels of a specific gene as a side track. 
                                                       # If this is not null, you must load groupRNA data in project in advance
                                          cellTypeColors = NULL, # Named character vector. Names are cell type names while values are colors.
                                          cellTypeOrder = NULL, # Character vector. Specifies the order of cell types in the legend
                                          rowOrder = NULL, # Integer vector specifying the order of rows (pseudobulks).
                                          lineageGroups = NULL, # IDs of the pseudobulks to plot
                                          signalThreshold = 0.1, # Mask signal below this value
                                          footprintRadius = 20, # Radius of footprinting
                                          smoothRadius = 5, # Radius for smoothing the plot
                                          heatmapPalette = NULL, # Color palette for heatmap
                                          vmax = NULL, # Values above this value are capped to this value
                                          vmin = NULL, # Values below this value are floored to this value
                                          ...) 
           standardGeneric("plotFeatureHeatmap"))

setMethod("plotFeatureHeatmap", "footprintingProject", 
          function(project, 
                   regionInd, 
                   feature = "footprint",
                   gene = NULL,
                   cellTypeColors = NULL,
                   cellTypeOrder = NULL,
                   rowOrder = NULL,
                   lineageGroups = NULL,
                   signalThreshold = 0.1,
                   footprintRadius = 20,
                   smoothRadius = NULL,
                   heatmapPalette = NULL,
                   vmax = NULL,
                   vmin = NULL){
            
            set.seed(42)
            library(ComplexHeatmap)
            library(BuenColors)
            library(circlize)
            library(RColorBrewer)
            
            width <- regionWidth(project)
            if(is.null(lineageGroups)){
              groupIDs <- groups(project)
            }else{
              groupIDs <- lineageGroups
            }
            
            if(is.null(smoothRadius)){
              smoothRadius <- as.integer(footprintRadius / 2)
            }
            
            # Get data matrix for plotting
            featureMatrix <- getFeatureMatrix(project = project,
                                              regionInd = regionInd,
                                              feature = feature,
                                              groupIDs = groupIDs,
                                              footprintRadius = footprintRadius,
                                              smoothRadius = smoothRadius)
            
            # Load total accessibility in this region for each pseudobulk
            ATACTrack <- groupATAC(project)[regionInd,]
            
            # Get cell type labels of pseudobulks
            cellTypes <- groupCellType(project)
            if(length(cellTypes) == 0){cellTypes <- groups(project)}
            names(cellTypes) <- groups(project)
            
            # Assign colors to cell types if no palettes are provided
            if(is.null(cellTypeColors)){
              cellTypeColors <- as.character(MetBrewer::met.brewer("Austria", length(cellTypes)))
              names(cellTypeColors) <- cellTypes
            }
            
            # Filter out low signal footprints
            featureMatrix[featureMatrix < signalThreshold] <- 0
            
            # Get row order
            if(is.null(rowOrder)){
              rowOrder <- order(cellTypes)
            }
            
            # Get cell type order
            if(is.null(cellTypeOrder)){
              cellTypeOrder <- gtools::mixedsort(unique(cellTypes))
            }
            
            # If a lineage is specified, retrieve data for the lineage and order by pseudotime
            if(!is.null(lineageGroups)){
              ATACTrack <- ATACTrack[lineageGroups]
              cellTypes <- cellTypes[lineageGroups]
              cellTypeColors <- cellTypeColors[unique(cellTypes)]
              if(length(groupPseudoTime(project)) > 0){
                lineagePseudoTime <- groupPseudoTime(project)[lineageGroups]
                rowOrder <- order(lineagePseudoTime)
              }else{
                rowOrder <- 1:length(lineageGroups)
              }
            }
            
            # If gene is not null, retrieve RNA data for the gene
            if(!is.null(gene)){
              # Get RNA of pseudobulks
              RNATrack <- groupRNA(project)[gene, ]
              if(!is.null(lineageGroups)){
                RNATrack <- RNATrack[lineageGroups]
              }
            }else{
              RNATrack <- NULL
            }
            
            plotFeatureHeatmapHelper(project, 
                                     featureMatrix = featureMatrix,
                                     feature = feature,
                                     ATACTrack = ATACTrack,
                                     RNATrack = RNATrack,
                                     gene = gene,
                                     cellTypes = cellTypes,
                                     cellTypeColors = cellTypeColors,
                                     cellTypeOrder = cellTypeOrder,
                                     heatmapPalette = heatmapPalette,
                                     rowOrder = rowOrder,
                                     vmax = vmax,
                                     vmin = vmin)
            
          })

# Helper function for plotFeatureHeatmap
plotFeatureHeatmapHelper <- function(project, 
                                     featureMatrix,
                                     feature = "footprint",
                                     ATACTrack = NULL,
                                     RNATrack = NULL,
                                     gene = NULL,
                                     cellTypes = NULL,
                                     cellTypeColors = NULL,
                                     cellTypeOrder = NULL,
                                     heatmapPalette = NULL, 
                                     rowOrder = NULL,
                                     vmax = NULL,
                                     vmin = NULL){
  
  # Get max color value
  if(is.null(vmax)){
    if(quantile(featureMatrix, 0.95) > 0){
      vmax <- quantile(featureMatrix, 0.95)
    }else{
      vmax <- max(featureMatrix)
    }
    if(feature == "TFBS"){
      vmax <- 0.8
    }
  }
  
  # Get min color value
  if(is.null(vmin)){
    
    vmin <- quantile(featureMatrix, 0.5)
    
    if(feature == "TFBS"){
      vmin <- 0
    }
  }
  
  # Get colors for the heatmap
  if(is.null(heatmapPalette)){
    heatmapPalette <- colorRampPalette(c(rep("white", 2),  "#9ECAE1", "#08519C", "#08306B"))(9)
  }
  heatmapColors <- colorRamp2(seq(vmin, vmax, length.out=9),
                              colors = heatmapPalette)

  # Get colors for the ATAC side bar
  if(quantile(ATACTrack, 0.99) > 0){
    ATACColors <- colorRamp2(seq(0, quantile(ATACTrack, 0.99),length.out=9),
                             colors = jdb_palette("brewer_green"))
  }else{
    ATACColors <- colorRamp2(seq(0, max(ATACTrack),length.out=9),
                             colors = jdb_palette("brewer_green"))
  }
  
  # Use cell type labels as left side bar
  leftAnno <- rowAnnotation(cellType=factor(cellTypes[rowOrder], levels = cellTypeOrder),
                            col=list(cellType=cellTypeColors),
                            border=TRUE,show_annotation_name=FALSE)
  
  # Get colors for feature-RNA correlation top bar
  corColors <- colorRamp2(seq(-1, 1, length.out=9),colors = colorRampPalette(c("blue", "white", "red"))(9))
  
  # Prepare RNA data and calculate RNA-feature correlation
  if(!is.null(RNATrack)){
    
    # Ger colors for the RNA side bar
    if(quantile(RNATrack, 0.99) > 0){
      RNAColors <- colorRamp2(seq(quantile(RNATrack, 0.01),
                                  quantile(RNATrack, 0.95),length.out=9),
                              colors = jdb_palette("brewer_purple"))
    }else{
      RNAColors <- colorRamp2(seq(quantile(RNATrack, 0.01),
                                  max(RNATrack),length.out=9),
                              colors = jdb_palette("brewer_purple"))
    }
    
    # Get the  feature-RNA correlation top bar
    RNACorr <- cor(t(featureMatrix), RNATrack, method = "spearman")
    RNACorr[is.na(RNACorr)] <- 0
    topAnno <- HeatmapAnnotation("Correlation" = RNACorr, 
                                 col = list(Correlation=corColors))
    
    # Get ATAC and RNA side bars
    rightAnno <- rowAnnotation("ATAC" = ATACTrack[rowOrder],
                               "RNA" = RNATrack[rowOrder],border=TRUE,
                               gap = unit(2, "points"),
                               col=list(ATAC=ATACColors,RNA=RNAColors))
  }else{
    topAnno <- NULL
    rightAnno <- rowAnnotation("ATAC" = ATACTrack[rowOrder],
                               gap = unit(2, "points"),
                               col=list(ATAC=ATACColors))
  }
  
  featureName <- list(
    "insertion" = "Tn5 insertion",
    "footprint" = "Footprint\nscore",
    "TFBS" = "TF binding\nscore"
  )[[feature]]
  
  # If rowOrder is not provided, order the rows using hierarchical clustering
  if(is.null(rowOrder)){
    rowOrder <- hclust(dist(t(featureMatrix)))$order
  }
  
  Heatmap(t(featureMatrix)[rowOrder,],
          use_raster = TRUE,
          col = heatmapColors,
          cluster_rows = F,
          clustering_distance_rows = "pearson",
          show_row_dend = FALSE,
          cluster_columns = FALSE,
          name = featureName,
          right_annotation = rightAnno,
          left_annotation = leftAnno,
          border = TRUE,
          top_annotation = topAnno,
          column_names_gp = gpar(fontsize = 7),
          column_title = "Position (bp)",
          column_title_side = "bottom")
}

setGeneric("plotPseudotimeHeatmap", function(project, # footprintingProject object
                                             regionInd, # Integer. Index of the region to plot
                                             feature = "footprint", # Character. Feature to plot. Can be one of "footprint", "insertion", and "TFBS"
                                             lineageGroups, # IDs of the pseudobulks to plot
                                             gene = NULL, # Whether to also plot the RNA levels of a specific gene as a side track. 
                                             # If this is not null, you must load groupRNA data in project in advance
                                             cellTypeColors = NULL, # Named character vector. Names are cell type names while values are colors.
                                             footprintRadius = 20, # Radius of footprinting
                                             heatmapPalette = NULL, # Color palette for heatmap
                                             pseudoTimeWindowSize = 10, # We aggregate data using a running window across pseudotime. This is the size of the window
                                             vmax = NULL, # Max color value. Values above this value are capped to this value
                                             vmin = NULL # Min calor value.
                                             ) 
  standardGeneric("plotPseudotimeHeatmap"))

setMethod("plotPseudotimeHeatmap", "footprintingProject", 
          function(project, 
                   regionInd, 
                   feature = "footprint",
                   lineageGroups,
                   gene = NULL,
                   cellTypeColors = NULL,
                   footprintRadius = 20,
                   heatmapPalette = NULL,
                   pseudoTimeWindowSize = 10,
                   vmax = NULL,
                   vmin = NULL){
            
            set.seed(42)
            library(ComplexHeatmap)
            library(BuenColors)
            library(circlize)
            library(RColorBrewer)
            
            width <- regionWidth(project)
            smoothRadius <- as.integer(footprintRadius / 2)
            
            # Load Tn5 insertion count tensor for the corresponding chunk
            ret <- getCountData(project, regionInd)
            countData <- ret[["countData"]]
            adjustedRegionInd <- ret[["regionInd"]]
            
            # Get sequence bias
            seqBias <- regionBias(project)
            rownames(seqBias) <- 1:length(regionRanges(project))
            
            # Load background dispersion model
            dispersionModel <- dispModel(project, as.character(footprintRadius))
            
            # Get cell type labels
            cellTypes <- groupCellType(project)
            
            # Remove ties from pseudotime values and re-order pseudobulks by pseudo time
            lineagePseudoTime <- groupPseudoTime(project)[lineageGroups]
            lineagePseudoTime <- rank(lineagePseudoTime, ties.method = "random")
            lineageGroups <- lineageGroups[order(lineagePseudoTime)]
            
            if(is.null(lineageGroups)){
              groupIDs <- groups(project)
            }else{
              groupIDs <- lineageGroups
            }
            
            # Get position-by-pseudobulk ATAC insertion matrix
            regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, width[regionInd])
            
            # Get Tn5 bias for every bp in the current region
            Tn5Bias <- regionBias(project)[regionInd, ]
            
            # Create a sliding window along pseudotime
            pseudoTimeWindows <- lapply(1:(length(lineageGroups) - pseudoTimeWindowSize + 1),
                                        function(i){i:(i + pseudoTimeWindowSize - 1)})
            
            if(feature == "TFBS"){
              
              # Calculate TF binding patterns for each pseudobulk
              regionTFBSMatrix <- getRegionTFBSMatrix(project = project, 
                                                  regionInd = regionInd, 
                                                  groupIDs = groupIDs)
              
              # Smooth using running pseudo-time window
              featureMatrix <- sapply(
                pseudoTimeWindows,
                function(pseudoTimeWindow){
                  rowMeans(regionTFBSMatrix[,pseudoTimeWindow])
                }
              )
              
            }else if(feature == "footprint"){
              
              # For every pseudotime window, we pool the data of mega cells in the same window
              # Then we calculate the position-by-pseudoTimeWindow matrix of footprint scores
              featureMatrix <- sapply(
                pseudoTimeWindows,
                function(pseudoTimeWindow){
                  regionWindowATAC <- rowSums(regionATAC[,pseudoTimeWindow])
                  pvals <- footprintScoring(regionWindowATAC, 
                                            Tn5Bias,
                                            dispersionModel,
                                            footprintRadius = footprintRadius,
                                            flankRadius = footprintRadius)
                  scores <- -log10(pvals)
                  scores <- caTools::runmax(scores, 2 * smoothRadius)
                  scores <- conv(scores, smoothRadius) / (2 * smoothRadius)
                }
              )
            }else if(feature == "insertion"){
              featureMatrix <- sapply(
                1:dim(regionATAC)[2],
                function(i){
                  conv(regionATAC[,i], smoothRadius) / (2 * smoothRadius)
                }
              )
            }
            
            # Get total region accessibility for each pseudotime window
            ATACTrack <- sapply(
              pseudoTimeWindows,
              function(pseudoTimeWindow){
                sum(groupATAC(project)[regionInd, lineageGroups][pseudoTimeWindow])
              }
            )
            
            # If one pseudotime window contains pseudobulks of multiple cell types, choose the cell type
            # with the highest frequency as the cell type label for this window
            cellTypes <- sapply(
              pseudoTimeWindows,
              function(pseudoTimeWindow){
                names(sort(table(cellTypes[lineageGroups][pseudoTimeWindow]), decreasing = T)[1])
              }
            )
            
            if(!is.null(gene)){
              # Get RNA of pseudobulks
              RNATrack <- sapply(
                pseudoTimeWindows,
                function(pseudoTimeWindow){
                  sum(groupRNA(project)[gene, lineageGroups][pseudoTimeWindow])
                }
              )
            }else{
              RNATrack <- NULL
            }
            
            plotFeatureHeatmapHelper(project, 
                                     featureMatrix = featureMatrix,
                                     feature = feature,
                                     ATACTrack = ATACTrack,
                                     RNATrack = RNATrack,
                                     gene = gene,
                                     cellTypes = cellTypes,
                                     cellTypeColors = cellTypeColors,
                                     cellTypeOrder = sort(unique(cellTypes)),
                                     heatmapPalette = heatmapPalette,
                                     rowOrder = 1:dim(featureMatrix)[2],
                                     vmax = vmax,
                                     vmin = vmin)
            
          })

setGeneric("plotSegmentHeatmap", function(project, # footprintingProject object
                                          regionInd, # Integer. Index of the region to plot
                                          lineageGroups = NULL, # IDs of the pseudo-bulks to plot
                                          footprintThreshold = 0.1, # Signal below this threshold is set to zero
                                          smoothRadius = NULL, # Radius for smoothing
                                          windowSize = 10, # Window size for segmentation
                                          contextRadius = 100, # Input radius of the TFBS model
                                          ...) 
  standardGeneric("plotSegmentHeatmap"))

setMethod("plotSegmentHeatmap", "footprintingProject", 
          function(project, 
                   regionInd, 
                   lineageGroups = NULL, 
                   footprintThreshold = 0.1, 
                   smoothRadius = NULL, 
                   windowSize = 10, 
                   contextRadius = 100 
                   ){
            
            ######################## First calculate TF binding scores for the region across pseudo-bulks ########################
            
            if(is.null(lineageGroups)){
              groupIDs <- groups(project)
            }else{
              groupIDs <- lineageGroups
            }
            
            # Retrieve Tn5 insertion count tensor for the corresponding chunk
            ret <- getCountData(project, regionInd)
            countData <- ret[["countData"]]
            adjustedRegionInd <- ret[["regionInd"]]
            
            # Get position-by-pseudobulk ATAC insertion matrix
            regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, regionWidth(project)[regionInd])
            
            # Load TF binding prediction model
            TFBSModel <- TFBindingModel(project)
            scales <- TFBSModel$scales # Scales are footprint radii we use for multi-scale footprinting
            
            # Load the corresponding background dispersion models
            dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
            names(dispModels) <- scales
            
            # Retrieve the region of inetrest
            region <- regionRanges(project)[regionInd]
              
            # Calculate footprint scores of the current region
            regionTFBS <- getRegionTFBS(regionATAC = regionATAC,
                                        Tn5Bias = regionBias(project)[regionInd,],
                                        region = region,
                                        dispModels = dispModels,
                                        TFBSModel = TFBindingModel(project))

            ######################## Segmentation of CRE into SUCROs ########################
            
            # Footprint correlation map
            fpCor <- cor(t(regionTFBS$TFBSScores))
            fpCor[is.na(fpCor)] <- 0
            
            # Segment the enhancer region
            segmentation <- CRESegmentation(fpCor, windowSize = windowSize)
            boundaries <- segmentation[["Boundaries"]]
            boundaries <- sapply(boundaries, function(i){regionTFBS$position[i]})
            BI <- segmentation[["boundaryIndex"]]
            
            # Get categorical color maps for coloring segments
            library(RColorBrewer)
            colorPalette = brewer.pal.info[brewer.pal.info$category == 'qual',]
            colorVector = unique(unlist(mapply(brewer.pal, colorPalette$maxcolors, rownames(colorPalette))))
            
            # Get track of segment ID
            segmentIDs <- rep(0, dim(regionATAC)[1])
            for(i in 1:(length(boundaries) - 1)){
              segmentIDs[boundaries[i]:boundaries[i+1]] <- i 
            }
            
            # Get colors of each segment
            segmentColors <- colorVector[unique(segmentIDs) + 1]
            segmentIDs[segmentIDs == 0] <- ""
            names(segmentColors) <- unique(segmentIDs)
            
            ######################## Visualize segmentation results ########################
            
            # Go through each tile and fill in the TF binding score values into an empty matrix
            TFBSMat <- array(0, dim = dim(regionATAC))
            rownames(TFBSMat) <- rep("", dim(TFBSMat)[1])
            for(ind in 1:length(regionTFBS$position)){
              siteStart <- start(regionTFBS$sites[ind]) - start(region) + 1
              siteEnd <- end(regionTFBS$sites[ind]) - start(region) + 1
              for(pos in siteStart:siteEnd){
                TFBSMat[pos, ] <- pmax(TFBSMat[pos, ],regionTFBS$TFBSScores[ind,])
              }
            }
            
            # Set palette for plotting
            Colors=c("blue", "white", "red")
            Colors=colorRampPalette(Colors)(9)
            
            # Visualize results
            mat <- cor(t(TFBSMat))
            mat[is.na(mat)] <- 0
            
            # Remove edge (positions too close to the edge will not have valid prediction scores)
            start <- contextRadius + 1
            end <- dim(mat)[1] - contextRadius
            mat <- mat[start : end, start : end]
            segmentIDs <- segmentIDs[start : end]
            
            heat.cols <- colorRamp2(seq(-1, 1, length.out=9),
                                    colors = c("#133e7e", "#225bb2", "#447fdd", "#7db0ea",
                                               "#ffffff", "#ee956a", "#da6c42", "#973d21", "#6b200c"))
            topAnno <- HeatmapAnnotation(Segment = segmentIDs,
                                       col=list(Segment=segmentColors),
                                       border=TRUE,show_annotation_name=FALSE)
            Heatmap(t(mat),
                    use_raster = TRUE,
                    name = "Correlation",
                    col = heat.cols,
                    cluster_rows = F,
                    cluster_columns = F,
                    top_annotation = topAnno,
                    border = TRUE,
                    column_title = paste("Position (bp) in region", regionInd),
                    column_title_side = "bottom",
                    row_names_side = "left",
                    row_title = paste("Position (bp) in region", regionInd))
        
          })

# For a specific region, compute and visualize multi-scale footprints
setGeneric("plotMultiScale", function(project,# footprintingProject object
                                          regionInd, # Integer. Index of the region to plot
                                          lineageGroups = NULL, # IDs of the pseudobulks to plot
                                          kernelSizes = 2:100, # footprint radii to use for multi-scale footprinting
                                          smoothRadius = 5 # Radius for smoothing
                                          ) 
  standardGeneric("plotMultiScale"))


setMethod("plotMultiScale", "footprintingProject", 
          function(project, 
                   regionInd, 
                   lineageGroups = NULL,
                   kernelSizes = 2:100,
                   smoothRadius = 5){
            
            library(ComplexHeatmap)
            library(BuenColors)
            library(circlize)
            library(RColorBrewer)
            
            if(is.null(lineageGroups)){
              groupIDs <- groups(project)
            }else{
              groupIDs <- lineageGroups
            }
            
            # Load predicted Tn5 bias track
            Tn5Bias <- regionBias(project)[regionInd, ]
            width <- regionWidth(project)
            
            # Retrieve Tn5 insertion count tensor for the corresponding chunk
            ret <- getCountData(project, regionInd)
            countData <- ret[["countData"]]
            adjustedRegionInd <- ret[["regionInd"]]
            
            # Get position-by-pseudobulk ATAC insertion matrix
            regionATAC <- getRegionATAC(countData, adjustedRegionInd, groupIDs, width[regionInd])
            
            # Go through each footprint radius and compute multi-scale footprints
            kernelScanning <- pbmcapply::pbmcmapply(
              function(footprintRadius){
                footprintPvals <- footprintScoring(
                  Tn5Insertion = rowSums(regionATAC),
                  Tn5Bias = Tn5Bias,
                  dispersionModel = dispModel(project, as.character(footprintRadius)),
                  footprintRadius = footprintRadius,
                  flankRadius = footprintRadius
                )
                footprintScores <- -log10(footprintPvals)
                footprintScores <- caTools::runmax(footprintScores, smoothRadius)
                footprintScores <- conv(footprintScores, smoothRadius) / (2 * smoothRadius)
                footprintScores[1:footprintRadius] <- 0 # Mask left edge due to incomplete left flank 
                footprintScores[(width[regionInd] - footprintRadius):width[regionInd]] <- 0 # Mask right edge due to incomplete right flank
                footprintScores
              },
              kernelSizes,
              mc.cores = 16
            )
            
            footprintColors <- colorRamp2(seq(0.5, 2,length.out=9),
                                          colors = jdb_palette("brewer_blue"))
            
            colnames(kernelScanning) <- rep("", length(kernelSizes))
            colnames(kernelScanning)[kernelSizes %% 10 == 0] <- kernelSizes[kernelSizes %% 10 == 0] * 2
            
            positions <- 1:width[regionInd]
            rownames(kernelScanning) <- rep("", width[regionInd])
            rownames(kernelScanning)[positions %% 100 == 0] <- positions[positions %% 100 == 0] - width[regionInd] / 2
            
            plotPositions <- max(kernelSizes):(width[regionInd] - max(kernelSizes))
            Heatmap(t(kernelScanning[plotPositions, rev(1:length(kernelSizes))]),
                    use_raster = TRUE,
                    col = footprintColors,
                    cluster_rows = F,
                    show_row_dend = F,
                    cluster_columns = F,
                    name = "Footprint\nScore",
                    border = TRUE,
                    column_title = "Position (bp)",
                    column_title_side = "bottom",
                    column_names_rot = 0)
            
          })

# Plot position-by-pseudobulk heatmap of a selected feature for a specific CRE region
setGeneric("plotTracks", function(project, # footprintingProject object
                                  regionInd, # Integer. Index of the region to plot
                                  lineageGroups = NULL, # IDs of the pseudobulks to plot
                                  footprintRadius = 20, # Radius of footprinting
                                  contextRadius = 100,# Input radius of the TFBS model
                                  ...) 
  standardGeneric("plotTracks"))

setMethod("plotTracks", "footprintingProject", 
          function(project, 
                   regionInd,
                   lineageGroups = NULL,
                   footprintRadius = 50,
                   contextRadius = 100){
            
            # Get position-by-pseudobulk matrix of ATAC signal for the current region
            ret <- getCountData(project, regionInd)
            countData <- ret[["countData"]]
            chunkRegionInd <- ret[["regionInd"]]
            region <- regionRanges(project)[regionInd]
            
            if(!is.null(lineageGroups)){
              groupIDs <- lineageGroups
            }else{
              groupIDs <- groups(project)
            }
            
            # Get position-by-pseudobulk ATAC insertion matrix
            width <- width(region)
            regionATAC <- getRegionATAC(countData, chunkRegionInd, groupIDs, width)
            
            # Calculate footprint scores
            footprintPvals <- footprintScoring(
              Tn5Insertion = rowSums(regionATAC),
              Tn5Bias = regionBias(project)[regionInd, ],
              dispersionModel = dispModel(project, footprintRadius),
              footprintRadius = footprintRadius,
              flankRadius = footprintRadius
            )
            smoothRadius <- as.integer(footprintRadius / 2)
            footprintScores <- -log10(footprintPvals)
            footprintScores <- caTools::runmax(footprintScores, smoothRadius)
            footprintScores <- conv(footprintScores, smoothRadius) / (2 * smoothRadius)
            footprintScores[1:footprintRadius] <- 0
            footprintScores[(width - footprintRadius):width] <- 0
            
            # Calculate TFBS scores
            TFBSModel <- TFBindingModel(project) 
            scales <- as.character(TFBSModel$scales)
            dispModels <- lapply(scales, function(scale){dispModel(project, as.character(scale))})
            names(dispModels) <- scales
            sites <- GenomicRanges::slidingWindows(region, width = 1, step = 1)[[1]]
            sites$score <- 1
            TFBS <- getRegionTFBS(regionATAC = as.matrix(rowSums(regionATAC), ncol = 1),
                                  Tn5Bias = regionBias(project)[regionInd, ],
                                  region = region,
                                  dispModels = dispModels,
                                  TFBSModel = TFBSModel,
                                  sites = sites)
            TFBSScores <- rep(0, width)
            for(i in 1:length(TFBS$TFBSScores)){
              TFBSScores[start(TFBS$sites[i]):end(TFBS$sites[i]) - start(region)] <- TFBS$TFBSScores[i, 1]
            }
            
            # Determine plotting region
            plotRadius <- as.integer(width / 2)
            positions <- (-plotRadius + contextRadius + 1):(plotRadius - contextRadius)
            pltRange <- (contextRadius + 1):(width - contextRadius)
            
            # Plot Tn5 insertion counts
            insertion <- rowSums(regionATAC)[pltRange]
            barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                                  y1 = 0, y2 = insertion)
            p1 <- ggplot(barData) + 
              geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
                        fill = "#7F3F98") +
              scale_y_continuous(expand = c(0, 0)) + ylab("Tn5 insertion") +
              theme_classic() 
            
            # Plot ATAC coverage
            coverage <- conv(insertion, 200) / 400
            barData <- data.frame(x1 = positions - 2, x2 = positions + 2, 
                                  y1 = 0, y2 = coverage)
            p2 <- ggplot(barData) + 
              geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), 
                        fill = "#7F3F98") +
              scale_y_continuous(expand = c(0, 0)) + ylab("ATAC coverage") +
              theme_classic() 
            
            # Plot footprint scores
            plotData <- data.frame(
              "Position" = positions,
              "footprints" = footprintScores[pltRange],
              "TFBS" = TFBSScores[pltRange]
            )
            p3 <- ggplot(plotData) + 
              geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "footprints"), 
                          fill = "#69ACD5") + xlab("") + ylab("Footprint scores") +
              theme_classic()
            
            # Plot TFBS scores
            p4 <- ggplot(plotData) + 
              geom_ribbon(aes_string(x = "Position", ymin = 0, ymax = "TFBS"), 
                          fill = "#08519C") + xlab("") + ylab("TF binding scores") +
              theme_classic()
            
            library(patchwork)
            p1 / p2 / p3 / p4
          })

# Visualize marker feature (e.g., genes) signal across cell groups
dotPlot <- function(dataMtx, # Input data matrix of k features (such as genes) by n cells
                    # Row names should be feature names and column names should be cell IDs
                    markers, # Character array. Names of features we want to plot
                    groupLabels, # Character array. Group labels of each single cell
                    scaleFeature = T, # Whether to scale the values per feature across cells
                    palette = NULL, # Custom palette. 
                    vmax = NULL, # Value cap
                    vmin = NULL # Value floor
                    ){
  
  # Only keep data for selected marker features
  dataMtx <- dataMtx[markers, ]
  
  # Calculate marker-by-group matrix of detection percentage
  groups <- unique(groupLabels)
  pctMtx <- pbmcapply::pbmcmapply(
    function(group){
      groupData <- dataMtx[markers, groupLabels == group]
      percentage <- rowMeans(groupData > 0)
    },
    groups,
    mc.cores = 8
  ) * 100
  colnames(pctMtx) <- groups
  
  # Calculate marker-by-group matrix of average signal
  expMtx <- pbmcapply::pbmcmapply(
    function(group){
      groupData <- dataMtx[markers, groupLabels == group]
      exp <- rowMeans(groupData)
    },
    groups,
    mc.cores = 8
  )
  colnames(expMtx) <- groups
  
  # Scale each feature
  if(scaleFeature){
    expMtx <- expMtx / rowMeans(expMtx)
  }
  
  # Visualize results
  plotData <- as.data.frame(expand.grid(markers, groups))
  colnames(plotData) <- c("markers", "groups")
  plotData$exp <- sapply(1:dim(plotData)[1], function(i){expMtx[as.character(plotData$markers)[i], 
                                                                as.character(plotData$groups)[i]]})
  plotData$percent <- sapply(1:dim(plotData)[1], function(i){pctMtx[as.character(plotData$markers)[i], 
                                                                    as.character(plotData$groups)[i]]})
  if(!is.null(vmax)){ plotData$exp <- pmin(plotData$exp, vmax) }
  if(!is.null(vmin)){ plotData$exp <- pax(plotData$exp, vmin) }
  if(is.null(palette)){palette <- jdb_palette("brewer_blue")}
  
  p <- ggplot(plotData) +
    geom_point(aes(x = markers, y = groups, color = exp, size = percent)) +
    scale_color_gradientn(colors = palette) +
    theme_classic() + 
    RotatedAxis()
  
  p
}