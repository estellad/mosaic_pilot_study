library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"


# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath.spotclean, save_rds_name))

pdf(paste0(figpath, paste0(UMAP_name, "Cluster", ".pdf")), width = 7, height = 6)
print(DimPlot(seu, label=TRUE) + labs(y= "UMAP_2", x = "UMAP_1"))
dev.off()


# DLBCL post spotclean DE ----------------------------------------------------------------
seu_de <- SCTransform(seu, verbose = FALSE, assay = "Spatial") # try re-running SCTransform

# Or try DE on log counts

# seu_de <- PrepSCTFindMarkers(seu)
all.markers <- FindAllMarkers(object = seu_de, only.pos = TRUE)
all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 40) %>%
  ungroup() -> top20

p <- DoHeatmap(seu_de, features = top20$gene, size = 15, angle = 90) + NoLegend() + theme(axis.text = element_text(size = 30))

pdf(paste0(figpath, paste0("/DLBCL_Vis_DE_Heatmap_40.pdf")), width = 40, height = 30)
print(p)
dev.off()

write.csv(all.markers, paste0(savepath, paste0("/DLBCL_Integrated_Visium_DE.csv")))