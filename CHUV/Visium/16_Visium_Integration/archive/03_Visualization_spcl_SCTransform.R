library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
# disease = "breast"
disease = "lung"
# disease = "dlbcl"

# # Before spotclean ------------------------------------------------------
# datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
# savepath.filtered = paste0(datapath.filtered, "/Results")
# UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
# save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTnoSpotClean_filtered.rds")
# savepath = savepath.filtered

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean

# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath, save_rds_name))

pdf(paste0(savepath, paste0(UMAP_name, "Level4_decon_max", ".pdf")), width = 11, height = 12)
DimPlot(seu, label=TRUE, group.by="Level4_decon_max", pt.size = 1) + labs(y= "UMAP_2", x = "UMAP_1") + theme(legend.position = "bottom")
dev.off()

pdf(paste0(savepath, paste0(UMAP_name, "Region", ".pdf")), width = 8.5, height = 6)
DimPlot(seu, label=TRUE, group.by="Region") + labs(y= "UMAP_2", x = "UMAP_1")
dev.off()

pdf(paste0(savepath, paste0(UMAP_name, "Cluster", ".pdf")), width = 7, height = 6)
print(DimPlot(seu, label=TRUE) + labs(y= "UMAP_2", x = "UMAP_1"))
dev.off()

pdf(paste0(savepath, paste0(UMAP_name, "Sample", ".pdf")), width = 7, height = 6)
print(DimPlot(seu, label=TRUE, group.by="sample_id") + labs(y= "UMAP_2", x = "UMAP_1"))
dev.off()


# # DLBCL post spotclean DE ----------------------------------------------------------------
# seu_de <- PrepSCTFindMarkers(seu)
# all.markers <- FindAllMarkers(object = seu_de, only.pos = TRUE)
# all.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1) %>%
#   slice_head(n = 20) %>%
#   ungroup() -> top20
# 
# p <- DoHeatmap(seu_de, features = top20$gene) + NoLegend()
# 
# pdf(paste0(savepath, paste0("/DLBCL_DE_Heatmap.pdf")), width = 30, height = 30)
# print(p)
# dev.off()
# 
# write.csv(all.markers, paste0(savepath, paste0("/DLBCL_Integrated_Visium_DE.csv")))
