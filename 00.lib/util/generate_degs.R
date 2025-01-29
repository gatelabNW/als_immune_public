find_degs <- function(seurat, group_1, group_2, genes_to_test, logfc.threshold = 0.585, save_path = NULL, latent.vars = c(),
                      mode = "reg") {
  print(glue("Finding DEGs for {group_1} vs. {group_2}"))
  # join layers
  seurat <- JoinLayers(seurat)
  if(mode == "SCT"){
    print("Working with SCT counts.")
    # no need as we have only one SCT model
    seurat <- PrepSCTFindMarkers(seurat)
    # seurat[["LOL"]] <- CreateAssayObject(counts = seurat[["SCT"]]$counts)
    # #"Copy "SCT scale.data" to "LOL scale.data"
    # seurat <- SetAssayData(seurat, assay = "LOL", slot = "scale.data", new.data =  seurat[["SCT"]]$scale.data); gc()
    # seurat <- SetAssayData(seurat, assay = "LOL", slot = "data", new.data =  seurat[["SCT"]]$data); gc()
    # DefaultAssay(seurat) <- "LOL"
    # seurat <- RenameAssays(seurat, LOL = 'SCT2')
    degs<-FindMarkers(
      object = seurat,
      ident.1 = group_1,
      ident.2 = group_2,
      test.use = "MAST",
      logfc.threshold = logfc.threshold,
      features = genes_to_test,
      # min.pct = 0.1,
      # min.cells.group = 1,
      # min.cells.feature = 1,
      latent.vars = latent.vars,
      assay = "SCT",
      slot="data",
      recorrect_umi = FALSE
      # using SCT counts per this post https://github.com/satijalab/seurat/issues/3923
      # using recorrect_umi since we subset to ident.1 and ident.2 https://github.com/satijalab/seurat/issues/6427
    )
  }else if(mode == "reg"){
    print("Working with reg counts.")
    degs<-FindMarkers(
      object = seurat,
      ident.1 = group_1,
      ident.2 = group_2,
      test.use = "MAST",
      logfc.threshold = logfc.threshold,
      features = genes_to_test,
      # min.pct = 0.1,
      # min.cells.group = 1,
      # min.cells.feature = 1,
      latent.vars = latent.vars,
      assay = "RNA",
      slot="data"
    )
  }else{
    print("No mode selected. No degs generated!!!")
  }

  # remove ribo and mt gene
  degs[["gene_id"]] <- rownames(degs)
  degs$BH <- p.adjust(degs$p_val, method = "BH")

  if (!is.null(save_path)) {
    print(glue("Savings DEGs degs to {save_path}"))
    fwrite(degs, file = save_path, row.names = TRUE)
    print("Success")
  }
  degs
}