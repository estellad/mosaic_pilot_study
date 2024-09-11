library(BayesSpace)
library(patchwork)
library(ggplot2)

disease = "Breast"
sample = "B3_2"
# resolution = "spots"
# resolution = "subspots"
read_path = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", disease, "/", sample)

# if(resolution == "spot"){ # spot
sce <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered.rds"))) # spot
assay(sce, "log1p") <- log1p(counts(sce))
# raw_assay_name = "counts"

# }else{ # subspot
sce_enhanced <- readRDS(file.path(read_path, paste0(sample, "_baye_clustered_all_enhanced_expr.rds"))) # subspot
# raw_assay_name = "decont_subspot"
# }

dim(rowData(sce_enhanced))

# -------------------------------------------------------------------------
# markers <- c("MS4A1", "CXCL13", "CXCR5", "CCL19")
# markers <- c("CD68", "CD163", "CSF1R", "ITGAM") # General macrophage
markers <- c("CD68", "TNFAIP1") # General macro, M1 (pro-inflam), M2 (pro-tumor)


plot_expression <- function(sce_obj = sce, marker = "CD4"){
  featurePlot(sce_obj, marker, assay.type = "log1p", color=NA) +
    scale_x_reverse() + 
    # viridis::scale_fill_viridis(option="A") +
    scale_color_gradientn(colors = colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, 
                                                                                              name = "Spectral")))(100)) + 
    labs(title=marker, fill="Log\nexpression") +
    theme(plot.title = element_blank())
  # theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}

# With legend
feature.plots <- purrr::map(markers, function(x) plot_expression(sce, x)) 
enhanced.feature.plots <- purrr::map(markers, function(x) plot_expression(sce_enhanced, x))

# Without legend
feature.plots2 <- purrr::map(markers, function(x) plot_expression(sce, x) + theme(legend.position = "none")) 
enhanced.feature.plots2 <- purrr::map(markers, function(x) plot_expression(sce_enhanced, x) + theme(legend.position = "none"))

# -------------------------------------------------------------------------
saveFeaturePlots <- function(feature.plots, enhanced.feature.plots, legend, width, height){
  p1 <- patchwork::wrap_plots(feature.plots, ncol = 4) + plot_layout(guides = "collect") 
  p2 <- patchwork::wrap_plots(enhanced.feature.plots, ncol = 4) + plot_layout(guides = "collect") 
  
  fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5"
  pdf(file = file.path(paste0(fig_path, "/L1/"),
                       paste0(sample, "_TLS_subspot_spot_", legend, ".pdf")),
      width = width,
      height = height)
  print(p2 / p1)
  dev.off()
}

saveFeaturePlots(feature.plots1, enhanced.feature.plots1, "legend", width = 10, height = 15)
saveFeaturePlots(feature.plots2, enhanced.feature.plots2, "nolegend", width = 10, height = 4.5)

# plotSpots(sce_enhanced, annotate = "CD68", in_tissue = NULL, assay = "log1p",
#           x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres")

