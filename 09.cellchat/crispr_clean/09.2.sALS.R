source("../00.ref/config/CRISPR_clean_config.R")
library(CellChat)

seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_04_05.rds")
cc_out_dir <- glue("{out_dir_root}/09.cellchat")
dir.create(cc_out_dir, showWarnings = FALSE, recursive = TRUE)
CellChatDB <- CellChatDB.human

s <- readRDS(input_seurat)
s@meta.data[["samples"]] <- s@meta.data$orig.ident

s <- subset(s, subset = diagnosis_general!="healthy_control")
als_s <- subset(s, subset = diagnosis!="als_c9orf72")
print(unique(als_s$diagnosis))
Idents(als_s) <- "predicted.celltype.l2"

cc_als <- createCellChat(object = als_s)
cc_als@DB <- CellChatDB
cc_als <- CellChat::subsetData(cc_als) # This step is necessary even if using the whole database
cc_als <- identifyOverExpressedGenes(cc_als)
cc_als <- identifyOverExpressedInteractions(cc_als)
cc_als <- computeCommunProb(cc_als, type = "triMean")
cc_als <- filterCommunication(cc_als, min.cells = 10)
df.net <- subsetCommunication(cc_als)
cur_out_file <- glue("{cc_out_dir}/sALS.csv")
write.csv(df.net, file = cur_out_file, quote = F, row.names = F)
cc_als <- computeCommunProbPathway(cc_als)
cc_als <- aggregateNet(cc_als)
cc_als <- netAnalysis_computeCentrality(cc_als, slot.name = "netP")
saveRDS(file = glue("{cc_out_dir}/sALS.rds"), cc_als)
groupSize <- as.numeric(table(cc_als@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/sALS_num_int.pdf",
    width = 12, height = 12)
p_num <- netVisual_circle(cc_als@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/sALS_weight_int.pdf",
    width = 12, height = 12)
p_weight <- netVisual_circle(cc_als@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/sALS_heatmap.pdf",
    width = 12, height = 12)
netVisual_heatmap(cc_als, color.heatmap = "Reds")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/sALS_bubble.pdf",
    width = 12, height = 12)
netAnalysis_signalingRole_scatter(cc_als)
dev.off()