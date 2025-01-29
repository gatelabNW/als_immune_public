

source("../../00.ref/config/CRISPR_clean_config.R")
set_size_min <- 3

DEG_out_dir <- glue("{out_dir_root}/05.DEG_430_filter")
color_df<-read.csv(color_path)
celltype2color<-color_df$color
names(celltype2color)<-color_df$predicted.celltype.l2


###
# deg_threshold<-0.01
# fc_threshold<-0.585
# diagnosis_general_out_root <-glue("{DEG_out_dir}/degs/diagnosis_general/SCT")
# all_DEG_dirs <- list.dirs(diagnosis_general_out_root, recursive=FALSE)
# for(degs_path_root in all_DEG_dirs){
#   degs_path<-paste(degs_path_root,list.files(degs_path_root),sep="/")
#   gene_sets<-list()
#   for(cur_file in degs_path){
#     x<-tail(unlist(strsplit(cur_file,"/")),n=1)%>%
#       strsplit("_vs._")%>%
#       unlist()
#     if(!grepl("__", cur_file, fixed = TRUE)){
#       next
#     }
#     cur_cell_type<-strsplit(unlist(strsplit(x[1],"___"))[2], "__")[[1]][1]
#     df<-read.csv(cur_file)
#     df<-df[df["BH"]<deg_threshold,]
#     df<-df[df["avg_log2FC"]>fc_threshold | df["avg_log2FC"]< (-1)*fc_threshold,]
#
#     # no deg
#     if(nrow(df)!=0){
#       gene_sets[[cur_cell_type]]<-df[["X"]]
#     }
#
#   }
#   non_zero_count<-lapply(gene_sets,function(x){length(x)>0})%>%
#     unlist()%>%
#     sum()
#   gene_sets<-gene_sets[lengths(gene_sets)>set_size_min]
#
#   t<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq")
#   set_names<-t$labels
#   new_col<-sapply(set_names,function(x) celltype2color[[str_replace(x,"_"," ")]])
#   if(non_zero_count>1){
#     out_dir<-degs_path_root|>
#       str_replace("degs", "plots")
#     out_file <- glue("{out_dir}/upset_plot.pdf")
#     upset_plt<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq",sets.bar.color = new_col)
#     pdf(out_file, width = 10, height = 10)
#     print(upset_plt)
#     dev.off()
#     # grid.text(str_replace_all(condition_title,"_"," "),x = 0.65, y=0.95, gp=gpar(fontsize=20))
#   }
# }

###
deg_threshold<-0.01
fc_threshold<-0.585
diagnosis_general_out_root <-glue("{DEG_out_dir}/degs/sALS_hc/SCT")
all_DEG_dirs <- list.dirs(diagnosis_general_out_root, recursive=FALSE)
for(degs_path_root in all_DEG_dirs){
  degs_path<-paste(degs_path_root,list.files(degs_path_root),sep="/")
  gene_sets<-list()
  for(cur_file in degs_path){
    x<-tail(unlist(strsplit(cur_file,"/")),n=1)%>%
      strsplit("_vs._")%>%
      unlist()
    if(!grepl("__", cur_file, fixed = TRUE)){
      next
    }
    cur_cell_type<-strsplit(unlist(strsplit(x[1],"___"))[2], "__")[[1]][1]
    df<-read.csv(cur_file)
    df<-df[df["BH"]<deg_threshold,]
    df<-df[df["avg_log2FC"]>fc_threshold | df["avg_log2FC"]< (-1)*fc_threshold,]

    # no deg
    if(nrow(df)!=0){
      gene_sets[[cur_cell_type]]<-df[["X"]]
    }

  }
  non_zero_count<-lapply(gene_sets,function(x){length(x)>0})%>%
    unlist()%>%
    sum()
  gene_sets<-gene_sets[lengths(gene_sets)>set_size_min]

  t<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq")
  set_names<-t$labels
  new_col<-sapply(set_names,function(x) celltype2color[[str_replace(x,"_"," ")]])
  if(non_zero_count>1){
    out_dir<-degs_path_root|>
      str_replace("degs", "plots")
    out_file <- glue("{out_dir}/upset_plot.pdf")
    upset_plt<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq",sets.bar.color = new_col)
    pdf(out_file, width = 10, height = 10)
    print(upset_plt)
    dev.off()
    # grid.text(str_replace_all(condition_title,"_"," "),x = 0.65, y=0.95, gp=gpar(fontsize=20))
  }
}


###
deg_threshold<-0.01
fc_threshold<-0.585
c9_out_root <-glue("{DEG_out_dir}/degs/female_c9_hc/SCT")
all_DEG_dirs <- list.dirs(c9_out_root, recursive=FALSE)

for(degs_path_root in all_DEG_dirs){
  degs_path<-paste(degs_path_root,list.files(degs_path_root),sep="/")
  gene_sets<-list()
  for(cur_file in degs_path){
    x<-tail(unlist(strsplit(cur_file,"/")),n=1)%>%
      strsplit("_vs._")%>%
      unlist()
    if(!grepl("__", cur_file, fixed = TRUE)){
      next
    }
    cur_cell_type<-strsplit(unlist(strsplit(x[1],"___"))[2], "__")[[1]][1]
    df<-read.csv(cur_file)
    df<-df[df["BH"]<deg_threshold,]
    df<-df[df["avg_log2FC"]>fc_threshold | df["avg_log2FC"]< (-1)*fc_threshold,]

    # no deg
    if(nrow(df)!=0){
      gene_sets[[cur_cell_type]]<-df[["X"]]
    }

  }
  non_zero_count<-lapply(gene_sets,function(x){length(x)>0})%>%
    unlist()%>%
    sum()
  gene_sets<-gene_sets[lengths(gene_sets)>set_size_min]

  t<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq")
  set_names<-t$labels
  new_col<-sapply(set_names,function(x) celltype2color[[str_replace(x,"_"," ")]])
  if(non_zero_count>1){
    out_dir<-degs_path_root|>
      str_replace("degs", "plots")
    out_file <- glue("{out_dir}/upset_plot.pdf")
    upset_plt<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq",sets.bar.color = new_col)
    pdf(out_file, width = 10, height = 10)
    print(upset_plt)
    dev.off()
    # grid.text(str_replace_all(condition_title,"_"," "),x = 0.65, y=0.95, gp=gpar(fontsize=20))
  }
}

###
# deg_threshold<-0.01
# fc_threshold<-0.585
# diagnosis_out_root <-glue("{DEG_out_dir}/degs/diagnosis/SCT")
# all_DEG_cov_dirs <- list.dirs(diagnosis_out_root, recursive=FALSE)
# for(cur_cov_dir in all_DEG_cov_dirs){
#   all_DEG_files <- list.files(cur_cov_dir, recursive=FALSE, full.names = T)
#   for(cur_condition_string in c("als_fast_vs._healthy_control", "als_slow_vs._healthy_control", "als_fast_vs._als_slow")){
#     degs_path <- all_DEG_files[grepl(cur_condition_string, all_DEG_files)]
#     gene_sets<-list()
#     for(cur_file in degs_path){
#       x<-tail(unlist(strsplit(cur_file,"/")),n=1)%>%
#         strsplit("_vs._")%>%
#         unlist()
#       if(!grepl("__", cur_file, fixed = TRUE)){
#         next
#       }
#       cur_cell_type<-strsplit(unlist(strsplit(x[1],"___"))[2], "__")[[1]][1]
#       df<-read.csv(cur_file)
#       df<-df[df["BH"]<deg_threshold,]
#       df<-df[df["avg_log2FC"]>fc_threshold | df["avg_log2FC"]< (-1)*fc_threshold,]
#
#       # no deg
#       if(nrow(df)!=0){
#         gene_sets[[cur_cell_type]]<-df[["X"]]
#       }
#
#     }
#     non_zero_count<-lapply(gene_sets,function(x){length(x)>0})%>%
#       unlist()%>%
#       sum()
#     gene_sets<-gene_sets[lengths(gene_sets)>set_size_min]
#
#     t<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq")
#     set_names<-t$labels
#     new_col<-sapply(set_names,function(x) celltype2color[[str_replace(x,"_"," ")]])
#     if(non_zero_count>1){
#       out_dir<-cur_cov_dir|>
#         str_replace("degs", "plots")
#       out_file <- glue("{out_dir}/{cur_condition_string}___upset_plot.pdf")
#       upset_plt<-upset(fromList(gene_sets), nsets = length(gene_sets),nintersects=30, order.by = "freq",sets.bar.color = new_col)
#       pdf(out_file, width = 10, height = 10)
#       print(upset_plt)
#       dev.off()
#       # grid.text(str_replace_all(condition_title,"_"," "),x = 0.65, y=0.95, gp=gpar(fontsize=20))
#     }
#   }
# }