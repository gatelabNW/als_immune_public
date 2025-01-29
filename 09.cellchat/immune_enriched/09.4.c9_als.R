source("../../00.ref/config/immune_panel_config.R")
library(CellChat)


seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_SCT_mapped_04_02.rds")
cc_out_dir <- glue("{out_dir_root}/09.cellchat")
dir.create(cc_out_dir, showWarnings = FALSE, recursive = TRUE)
CellChatDB <- CellChatDB.human

s <- readRDS(input_seurat)
s@meta.data[["samples"]] <- s@meta.data$orig.ident

als_c9_s <- subset(s, diagnosis == "als_c9orf72")
Idents(als_c9_s) <- "predicted.celltype.l2"

cc_als_c9 <- createCellChat(object = als_c9_s)
cc_als_c9@DB <- CellChatDB
cc_als_c9 <- CellChat::subsetData(cc_als_c9) # This step is necessary even if using the whole database
cc_als_c9 <- identifyOverExpressedGenes(cc_als_c9)
cc_als_c9 <- identifyOverExpressedInteractions(cc_als_c9)
cc_als_c9 <- computeCommunProb(cc_als_c9, type = "triMean")
cc_als_c9 <- filterCommunication(cc_als_c9, min.cells = 10)
df.net <- subsetCommunication(cc_als_c9)
cur_out_file <- glue("{cc_out_dir}/als_c9.csv")
write.csv(df.net, file = cur_out_file, quote = F, row.names = F)
cc_als_c9 <- computeCommunProbPathway(cc_als_c9)
cc_als_c9 <- aggregateNet(cc_als_c9)
cc_als_c9 <- netAnalysis_computeCentrality(cc_als_c9, slot.name = "netP")
saveRDS(file = glue("{cc_out_dir}/als_c9.rds"), cc_als_c9)
groupSize <- as.numeric(table(cc_als_c9@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/als_c9_num_int.pdf",
    width = 12, height = 12)
p_num <- netVisual_circle(cc_als_c9@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/als_c9_weight_int.pdf",
    width = 12, height = 12)
p_weight <- netVisual_circle(cc_als_c9@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/als_c9_heatmap.pdf",
    width = 12, height = 12)
netVisual_heatmap(cc_als_c9, color.heatmap = "Reds")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/als_c9_bubble.pdf",
    width = 12, height = 12)
netAnalysis_signalingRole_scatter(cc_als_c9)
dev.off()