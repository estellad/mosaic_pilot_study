disease = "Lung"
sample = "L4_2"
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

# -------------------------------------------------------------------------
markers <- list()
markers[["B-cell"]] <- c("CXCL13", "MS4A1")
markers[["T-cell"]] <- c("CD4")

markers <- unlist(markers)

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

feature.plots <- purrr::map(markers, function(x) plot_expression(sce, x)) 
enhanced.feature.plots_no_legend <- purrr::map(markers[1:2], function(x) plot_expression(sce_enhanced, x) + theme(legend.position = "none"))
enhanced.feature.plots_w_legend <- purrr::map(markers[3], function(x) plot_expression(sce_enhanced, x))
# enhanced.feature.plots_w_legend <- list(plot_expression(sce_enhanced, markers[3]))
enhanced.feature.plots <- append(enhanced.feature.plots_no_legend, enhanced.feature.plots_w_legend)

spot_feature_counts <- purrr::map(markers, function(x) each_counts(sce, x))
subspot_feature_counts <- purrr::map(markers, function(x) each_counts(sce_enhanced, x))
# -------------------------------------------------------------------------
p1 <- patchwork::wrap_plots(feature.plots, ncol=3) + plot_layout(guides = "collect") & 
  scale_colour_continuous(limits = range(unlist(spot_feature_counts)))

# plot_title = "Spot_TLS_logexpr.pdf"
# pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig2e_TLS/", plot_title),
#     width = 8,
#     height = 6)
# print(p1)
# dev.off()

p2 <- patchwork::wrap_plots(enhanced.feature.plots, ncol=3) + plot_layout(guides = "collect") & 
  scale_colour_continuous(limits = range(unlist(subspot_feature_counts)))

# plot_title = "Subspot_TLS_logexpr.pdf"
# pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig2e_TLS/", plot_title),
#     width = 8,
#     height = 6)
# print(p2)
# dev.off()

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig5"
pdf(file = file.path(paste0(fig_path, "/L4/"),
                     paste0(sample, "_TLS_subspot_spot.pdf")),
    width = 10,
    height = 6)
print(p2 / p1)
dev.off()




