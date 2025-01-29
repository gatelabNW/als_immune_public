source("../00.ref/config/CRISPR_clean_config.R")
library(CellChat)

seurat_output_dir <- glue("{out_dir_root}/04.seurat")
input_seurat <- glue("{seurat_output_dir}/seurat_04_05.rds")
cc_out_dir <- glue("{out_dir_root}/09.cellchat")
dir.create(cc_out_dir, showWarnings = FALSE, recursive = TRUE)
CellChatDB <- CellChatDB.human

s <- readRDS(input_seurat)
s@meta.data[["samples"]] <- s@meta.data$orig.ident

female_hc_s <- subset(s, subset = sex == "f"  & diagnosis == "healthy_control")
Idents(female_hc_s) <- "predicted.celltype.l2"

cc_female_hc <- createCellChat(object = female_hc_s)
cc_female_hc@DB <- CellChatDB
cc_female_hc <- CellChat::subsetData(cc_female_hc) # This step is necessary even if using the whole database
cc_female_hc <- identifyOverExpressedGenes(cc_female_hc)
cc_female_hc <- identifyOverExpressedInteractions(cc_female_hc)
cc_female_hc <- computeCommunProb(cc_female_hc, type = "triMean")
cc_female_hc <- filterCommunication(cc_female_hc, min.cells = 10)
df.net <- subsetCommunication(cc_female_hc)
cur_out_file <- glue("{cc_out_dir}/female_hc.csv")
write.csv(df.net, file = cur_out_file, quote = F, row.names = F)
cc_female_hc <- computeCommunProbPathway(cc_female_hc)
cc_female_hc <- aggregateNet(cc_female_hc)
cc_female_hc <- netAnalysis_computeCentrality(cc_female_hc, slot.name = "netP")
saveRDS(file = glue("{cc_out_dir}/female_hc.rds"), cc_female_hc)
groupSize <- as.numeric(table(cc_female_hc@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/female_hc_num_int.pdf",
    width = 12, height = 12)
p_num <- netVisual_circle(cc_female_hc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/female_hc_weight_int.pdf",
    width = 12, height = 12)
p_weight <- netVisual_circle(cc_female_hc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/female_hc_heatmap.pdf",
    width = 12, height = 12)
netVisual_heatmap(cc_female_hc, color.heatmap = "Reds")
dev.off()

pdf("/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/09.cellchat/female_hc_bubble.pdf",
    width = 12, height = 12)
netAnalysis_signalingRole_scatter(cc_female_hc)
dev.off()