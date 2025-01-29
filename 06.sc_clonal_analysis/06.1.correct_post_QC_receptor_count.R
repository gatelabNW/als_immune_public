####################################################################
# Charles Zhang
# Gate Lab
# Northwestern University
# 10/12/2023
####################################################################
# Perform correction on seurat object's receptor counts
# Expected output:
#       seurat_06_01.rds
####################################################################
# See README.md for additional links

source("../00.ref/config/immune_profiling_config.R")
seurat_output_dir <- glue("{out_dir_root}/04.seurat")
seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_03.rds")|>read_rds()
receptor_seurat_output_dir <- glue("{out_dir_root}/06.seurat")
out_file <- glue("{receptor_seurat_output_dir}/seurat_06_01.rds")
dir.create(receptor_seurat_output_dir, showWarnings = FALSE, recursive = TRUE)

# build count dictionary
mt <- seurat@meta.data

b_freq_df<-dplyr::select(seurat@meta.data, bcr_clonotype_id)
is.na(b_freq_df)|>sum()|>print()
b_freq_df<-b_freq_df[!is.na(b_freq_df$bcr_clonotype_id),]
t_freq_df<-dplyr::select(seurat@meta.data, tcr_clonotype_id)
is.na(t_freq_df)|>sum()|>print()
t_freq_df<-t_freq_df[!is.na(t_freq_df$tcr_clonotype_id),]

b_freq_dict<-table(b_freq_df)
t_freq_dict<-table(t_freq_df)

# update meta data
updated_mt<-seurat@meta.data
updated_mt <- updated_mt[,-grep("pANN", names(updated_mt))]
updated_mt[["bcr_frequency"]]<-sapply(updated_mt$bcr_clonotype_id, function(x){
  if(x %in% names(b_freq_dict)){
    out<-b_freq_dict[[x]]
  }else{
    out<-NA
  }
  out
})

updated_mt[["tcr_frequency"]]<-sapply(updated_mt$tcr_clonotype_id, function(x){
  if(x %in% names(t_freq_dict)){
    out<-t_freq_dict[[x]]
  }else{
    out<-NA
  }
  out
})

# check bcr
clonotypes <- sample(b_freq_dict,1000)|>names()
for(cur_clonotype in clonotypes){
  cur_mt <- dplyr::filter(updated_mt, bcr_clonotype_id == cur_clonotype)
  cur_count <- cur_mt$bcr_frequency
  if(nrow(cur_mt)!=unique(cur_count)){
    break
  }
}

# check tcr
clonotypes <- sample(t_freq_dict,1000)|>names()
for(cur_clonotype in clonotypes){
  cur_mt <- dplyr::filter(updated_mt, tcr_clonotype_id == cur_clonotype)
  cur_count <- cur_mt$tcr_frequency
  if(nrow(cur_mt)!=unique(cur_count)){
    break
  }
}

# replace old meta data
seurat@meta.data <- updated_mt

# save new object
saveRDS(seurat, out_file)