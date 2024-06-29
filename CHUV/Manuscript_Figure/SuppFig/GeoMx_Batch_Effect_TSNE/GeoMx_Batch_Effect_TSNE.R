library(ggspavis)
library(patchwork)
library(SpatialExperiment)
geopath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/GeoMx_Batch_Effect_TSNE"

disease = "breast"
disease = "lung"
disease = "dlbcl"
spe <- readRDS(file.path(geopath, paste0(disease, "_spe.rds")))
spe_ruv <- readRDS(file.path(geopath, paste0(disease, "_spe_ruv.rds")))

## In total 8 slides across three indications 
library(RColorBrewer)
n <- 8
color_slides <- brewer.pal(n, "Set2")
pie(rep(1,n), col=color_slides)
palette <- list(
  breast = color_slides[1:2],
  lung = color_slides[3:5],
  dlbcl = color_slides[6:8]
)

## Cell fraction color
if(disease %in% c("breast", "lung")){
  cell_fraction_order = c("Macro", "Malignant", "Other", "T cells", "PanCK-")
  cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
  names(cell_fraction_color) <- cell_fraction_order
}else{
  cell_fraction_order = c("B cells", "Macro", "Other", "T cells")
  cell_fraction_color = c("#FFD700", "#9A32CD", "#388E8E", "#4169E1")
  names(cell_fraction_color) <- cell_fraction_order
}


# -------------------------------------------------------------------------
plotDimRed_publication <- function(spe, reddim = "PCA", annotate = "section_id", legen.lab = "Section"){
  p <- plotDimRed(spe, type = reddim, annotate = annotate, pt.size = 3) + 
    theme(legend.position = "bottom",
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks = element_blank(),
          legend.title= element_text(size = 10),
          legend.key.size = unit(0.8, "lines")) +
    guides(colour = guide_legend(nrow = 1)) + 
    labs(color = legen.lab)
  p
}

plot_save_reddim <- function(reddim = "PCA", plot_title = "geo_lung_PCA.pdf"){
  sample_name <- ifelse(disease == "dlbcl", "patient", "section_id")
  
  p1 <- plotDimRed_publication(spe, reddim = reddim, annotate = "slide_name", legen.lab = "Slide") + 
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) + theme(legend.position = "none") +
    scale_color_manual(values = palette[[disease]]) |
    plotDimRed_publication(spe, reddim = reddim, annotate = sample_name, legen.lab = "Sample") + 
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) + theme(legend.position = "none") |
    plotDimRed_publication(spe, reddim = reddim, annotate = "cell_fraction", legen.lab = "Cell fraction") +
    theme(plot.margin = margin(5.5, 5.5, 20, 35)) + theme(legend.position = "none") +
    scale_color_manual(values = cell_fraction_color)
  
  p2 <- plotDimRed_publication(spe_ruv, reddim = reddim, annotate = "slide_name", legen.lab = "Slide") +
    theme(plot.margin = margin(20, 35, 5.5, 5.5)) +
    scale_color_manual(values = palette[[disease]]) |
    plotDimRed_publication(spe_ruv, reddim = reddim, annotate = sample_name, legen.lab = "Sample") +
    theme(plot.margin = margin(5.5, 35, 5.5, 5.5)) |
    plotDimRed_publication(spe_ruv, reddim = reddim, annotate = "cell_fraction", legen.lab = "Cell fraction") +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 35)) +
    scale_color_manual(values = cell_fraction_color)
  
  pdf(file = file.path(figpath, plot_title),
      width = 20,
      height = 10)
  print(p1 / p2)
  dev.off()
}

plot_save_reddim(reddim = "PCA", plot_title = paste0("geo_", disease, "_PCA.pdf"))
plot_save_reddim(reddim = "UMAP", plot_title = paste0("geo_", disease, "_UMAP.pdf"))
plot_save_reddim(reddim = "TSNE", plot_title = paste0("geo_", disease, "_TSNE.pdf"))





# # -------------------------------------------------------------------------
# p_combined <- (p_breast_1 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_1 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_1 + theme(plot.margin = margin(0, 0, 50, 50))) /
#   (p_breast_2 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_2 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_2 + theme(plot.margin = margin(0, 0, 50, 50))) /
#   (p_breast_3 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_3 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_3 + theme(plot.margin = margin(0, 0, 50, 50))) /
#   (p_breast_4 + theme(plot.margin = margin(0, 50, 0, 0)) | p_lung_4 + theme(plot.margin = margin(0, 50, 0, 50)) | p_dlbcl_4 + theme(plot.margin = margin(0, 0, 0, 50))) 
# 
# figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"
# pdf(file.path(figpath, "Visium_UMAP_Samples.pdf"), width = 26, height = 30)
# print(p_combined)
# dev.off()
# 
# 
# # pdf(file = "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig1d/geo_lung_batch_PCA_cellfraction.pdf",
# #     width = 5,
# #     height = 5)
# # print(p)
# # dev.off()
