library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"


# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath.spotclean, save_rds_name))
# seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

seu_de <- PrepSCTFindMarkers(seu)
all.markers <- FindAllMarkers(object = seu_de, only.pos = TRUE)

all.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 40) %>%
  ungroup() -> top20

p <- DoHeatmap(seu_de, features = top20$gene, size = 15, angle = 0) + 
  NoLegend() + 
  theme(axis.text = element_text(size = 30)) #+ 
p



