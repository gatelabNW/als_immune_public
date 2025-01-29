
run_doubletfinder <- function(s) {

  # TODO: Find out why this value was chosen
  # The csf_immunity_age project uses 0.010
  # I think I remember Natalie found some documentation about the procedure used to generate my data.
  # IIRC, those docs said to use 0.044
  # TODO: Link to this doc in a comment
  # TODO: Also do this for all other params in this file
  doublet_formation_rate <- 0.044
  print(paste0("Using doublet formation rate of ", doublet_formation_rate))

  # Process normally
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Run TSNE clustering
  s <- RunPCA(s, features = VariableFeatures(object = s))
  s <- FindNeighbors(s, dims = 1:12)
  s <- FindClusters(s, resolution = 0.3)
  s <- RunTSNE(s, dims = 1:12)

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep_v3(s, PCs = 1:12, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats) # mean-variance-normalized bimodality coefficient, find maxBCmvn per pK
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_formation_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  print(s)

  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  s <- doubletFinder_v3(s, PCs = 1:12, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  ## Rename column name for consistency
  colnames(s@meta.data)[ grep("DF.classifications*", colnames(s@meta.data)) ] <- "DF.classifications"
  print(head(s@meta.data))

  return(s)
}