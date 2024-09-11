library(BayesSpace)
library(patchwork)
library(ggplot2)

disease = "Lung"
sample = "L4_2"
sample = "L1_4"
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
markers <- list()
# markers[["B-cell"]] <- c("CXCL13", "MS4A1")#, "CD23")
# markers[["T-cell"]] <- c("CD4", "PDCD1", "CD3E", "CD3D", "CD3G")
# markers[["Myeloid"]] <- c(# "CD11B", 
#                           "CD14", # "CD15", 
#                           "CD33", # "CD64",
#                           "CD68", # "CD117", 
#                           # "CD123", 
#                           "CD163" #, 
#                           # "CD169"
#                           )

# # B-cell
# markers[["B-cell"]] <- c("CXCL13", "MS4A1")
# # Follicular Dendritic Cell (FDC)
# markers[["FDC"]] <- c("FDCSP")
# # T follicular helper cells (TFH)
# markers[["TFH"]] <- c("PDCD1", "CD4", "CXCR5") # "CD3", 
# markers <- unlist(markers)

# c("CD19", "CD8", "CXCL9", "CXCL10", "CXCL11", 
#   "LAMP2", "CD209", "ITGAX", "CD74", "HHLA3", "HHLA2",
#   "CD40L", "IL21", "IL6",
#   "CCL19",
#   "PECAM1", "SELL", "CHST4",
#   "MKI67")[c("CD19", "CD8", "CXCL9", "CXCL10", "CXCL11", 
# "LAMP2", "CD209", "ITGAX", "CD74", "HHLA3", "HHLA2",
# "CD40L", "IL21", "IL6",
# "CCL19",
# "PECAM1", "SELL", "CHST4",
# "MKI67") %in% rownames(sce)]

markers <- c("CXCL9", "CXCL10", "CXCL11", "LAMP2", "CD209", "ITGAX", "CD74", "HHLA3", "HHLA2", "IL6", "CCL19", "PECAM1", "SELL", "CHST4", "MKI67")



plot_expression <- function(sce_obj = sce, marker = "CD4"){
  featurePlot(sce_obj, marker, assay.type = "log1p", color=NA) +
    viridis::scale_fill_viridis(option="A") +
    labs(title=marker, fill="Log\nexpression") +
    theme(plot.title = element_blank())
    # theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))
}

each_counts <- function(sce, feature) {
  assay(sce, "log1p")[feature, ]
}

feature.plots <- purrr::map(markers[8:15], function(x) plot_expression(sce, x)) 
# enhanced.feature.plots_no_legend <- purrr::map(markers[1:4], function(x) plot_expression(sce_enhanced, x) + theme(legend.position = "none"))
enhanced.feature.plots_w_legend <- purrr::map(markers[8:15], function(x) plot_expression(sce_enhanced, x))
enhanced.feature.plots <- enhanced.feature.plots_w_legend # append(enhanced.feature.plots_no_legend, enhanced.feature.plots_w_legend)

spot_feature_counts <- purrr::map(markers, function(x) each_counts(sce, x))
subspot_feature_counts <- purrr::map(markers, function(x) each_counts(sce_enhanced, x))
# -------------------------------------------------------------------------
p1 <- patchwork::wrap_plots(feature.plots, ncol = 8) + plot_layout(guides = "collect") & 
  scale_colour_continuous(limits = range(unlist(spot_feature_counts)))

p2 <- patchwork::wrap_plots(enhanced.feature.plots, ncol = 8) + plot_layout(guides = "collect") & 
  scale_colour_continuous(limits = range(unlist(subspot_feature_counts)))



# aggregate CD3 subtypes --------------------------------------------------
sce$CD3_spot <- colSums(assay(sce, "log1p")[c("CD3E", "CD3D", "CD3G"), ])
sce_enhanced$CD3_subspot <- colSums(assay(sce_enhanced, "log1p")[c("CD3E", "CD3D", "CD3G"), ])

library(ggspavis)
sce$pxl_col_in_fullres <- NULL
sce$pxl_row_in_fullres <- NULL
p_CD3 <- plotSpots(sce, annotate = "CD3_spot", x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", pt.size = 1) + scale_colour_viridis_c(option = "magma") + theme(plot.title = element_blank()
                                                                                                                                                                               , legend.position = "none")
psub_CD3 <- plotSpots(sce_enhanced, annotate = "CD3_subspot", x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", in_tissue = NULL) + scale_colour_viridis_c(option = "magma") + theme(plot.title = element_blank()
                                                                                                                                                                                                   , legend.position = "none")

psub_CD3 / p_CD3



fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5"
pdf(file = file.path(paste0(fig_path, "/L4/"),
                     paste0(sample, "_TLS_subspot_spot.pdf")),
    width = 10,
    height = 6)
print(p2 / p1)
dev.off()




