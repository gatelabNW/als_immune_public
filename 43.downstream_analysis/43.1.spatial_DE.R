source("../00.ref/config/spatial_config.R")
options(echo = TRUE)
source(glue("{repo_root}/00.lib/plot/ggplot_formatting.R"))
source(glue("{repo_root}/00.lib/plot/volcano_plot.R"))
source(glue("{repo_root}/00.lib/util/generate_degs.R"))

DEG_out_dir <- glue("{out_dir_root}/04.downstream_analysis/DE")
input_data_dir <- glue("{out_dir_root}/03.seurat_process")
out_data_dir <- glue("{DEG_out_dir}/data")
out_plot_dir <- glue("{DEG_out_dir}/plot")

dir.create(DEG_out_dir, recursive = T, showWarnings = F)
BH_thres <- 0.001
lfc_thres <- 0.585

s <- readRDS(glue("{input_data_dir}/data/all_samples_seurat_01.rds"))
DefaultAssay(s) <- "Spatial"
cluster_col <- "integrated_snn_res.0.25"
method <- "SCT"
clusters <- s@meta.data[[cluster_col]] |> unique() |> sort()


s[["condition_general"]] <- sapply(s$condition, function(x){
  if(x == "Control"){
    out <- x
  }else if(x == "sALS"){
    out <- "ALS"
  }else if(x == "C9orf72"){
    out <- "ALS"
  }else{
    out <- "NA"
  }
  out
})
condition_col <- "condition_general"
sample_conditions <- s@meta.data[[condition_col]] |> unique() |> sort()
comparison <- c("ALS", "Control")
s@meta.data$diagnosis_general_cluster <- glue('{s@meta.data[[cluster_col]]} - {s@meta.data[[condition_col]]}')
Idents(s) <- 'diagnosis_general_cluster'

for(cluster in clusters) {

  degs_path<-glue("{out_data_dir}/{comparison[1]}_vs._{comparison[2]}")
  plots_path<-glue("{out_plot_dir}/{comparison[1]}_vs._{comparison[2]}")
  dir.create(degs_path, recursive = T, showWarnings = F)
  dir.create(plots_path, recursive = T, showWarnings = F)

  title <- glue("cluster{cluster}: {comparison[1]} vs. {comparison[2]}")|> str_replace_all("_", " ")
  filename  <-  gsub(" ", "_", glue("cluster{cluster}__{comparison[1]}_vs._{comparison[2]}"))
  compA_cellnum<-sum(s$diagnosis_general_cluster==glue("{cluster} - {comparison[1]}"))
  compB_cellnum<-sum(s$diagnosis_general_cluster==glue("{cluster} - {comparison[2]}"))
  if(compA_cellnum<10 | compB_cellnum<10){
    print(glue("INFO: Not enough cells for {title}. Skip!"))
  }else{
    cur_deg_out_file <- glue("{degs_path}/{condition_col}___{filename}.csv")
    if(file.exists(cur_deg_out_file)) {
      if(overwrite){
        degs <- find_degs(s,
                          save_path = cur_deg_out_file,
                          group_1 = glue("{cluster} - {comparison[1]}"),
                          group_2 = glue("{cluster} - {comparison[2]}"),
                          logfc.threshold = -Inf,
                          # latent.vars = NULL,
                          mode = method
        )
      }else{
        print("INFO: Overwrite false. Only remaking volcano plot!")
        degs <- read.csv(cur_deg_out_file, row.names = 1)
      }
    }else{
      degs <- find_degs(s,
                        save_path = cur_deg_out_file,
                        group_1 = glue("{cluster} - {comparison[1]}"),
                        group_2 = glue("{cluster} - {comparison[2]}"),
                        logfc.threshold = -Inf,
                        # latent.vars = cur_cov,
                        mode = method
      )
    }
    degs[["gene"]] <- rownames(degs)
    volcano_plot(degs,
                 file = glue("{plots_path}/{condition_col}___{filename}.pdf"),
                 title = title,
                 padj.thresh = BH_thres,
                 lfc.threshold = lfc_thres
    )
  }
}

condition_col <- "condition"
s@meta.data$diagnosis_cluster <- glue('{s@meta.data[[cluster_col]]} - {s@meta.data[[condition_col]]}')
# generate all comparison permutations
clusters <- s@meta.data[[cluster_col]] |> unique() |> sort()
sample_conditions <- s@meta.data[[condition_col]] |> unique() |> sort()
Idents(s) <- 'diagnosis_cluster'

comparisons<-list(
  "1"=c("C9orf72", "Control"),
  "2"=c("sALS", "Control")
)
all_cov <- "no_cov"
for(cluster in clusters) {
  for(comparison in comparisons) {

    degs_path<-glue("{out_data_dir}/{comparison[1]}_vs._{comparison[2]}")
    plots_path<-glue("{out_plot_dir}/{comparison[1]}_vs._{comparison[2]}")
    dir.create(degs_path, recursive = T, showWarnings = F)
    dir.create(plots_path, recursive = T, showWarnings = F)

    title <- glue("cluster{cluster}: {comparison[1]} vs. {comparison[2]}")|> str_replace_all("_", " ")
    filename  <-  gsub(" ", "_", glue("cluster{cluster}__{comparison[1]}_vs._{comparison[2]}"))
    compA_cellnum<-sum(s$diagnosis_cluster==glue("{cluster} - {comparison[1]}"))
    compB_cellnum<-sum(s$diagnosis_cluster==glue("{cluster} - {comparison[2]}"))
    if(compA_cellnum<10 | compB_cellnum<10){
      print(glue("INFO: Not enough cells for {title}. Skip!"))
    }else{
      cur_deg_out_file <- glue("{degs_path}/{condition_col}___{filename}.csv")
      if(file.exists(cur_deg_out_file)) {
        if(overwrite){
          degs <- find_degs(s,
                            save_path = cur_deg_out_file,
                            group_1 = glue("{cluster} - {comparison[1]}"),
                            group_2 = glue("{cluster} - {comparison[2]}"),
                            logfc.threshold = -Inf,
                            # latent.vars = NULL,
                            mode = method
          )
        }else{
          print("INFO: Overwrite false. Only remaking volcano plot!")
          degs <- read.csv(cur_deg_out_file, row.names = 1)
        }
      }else{
        degs <- find_degs(s,
                          save_path = cur_deg_out_file,
                          group_1 = glue("{cluster} - {comparison[1]}"),
                          group_2 = glue("{cluster} - {comparison[2]}"),
                          logfc.threshold = -Inf,
                          # latent.vars = cur_cov,
                          mode = method
        )
      }
      degs[["gene"]] <- rownames(degs)
      volcano_plot(degs,
                   file = glue("{plots_path}/{condition_col}___{filename}.pdf"),
                   title = title,
                   padj.thresh = BH_thres,
                   lfc.threshold = lfc_thres
      )
    }
  }
}
