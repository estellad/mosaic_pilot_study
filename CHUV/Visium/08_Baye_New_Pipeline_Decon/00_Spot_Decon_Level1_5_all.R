# disease = "lung"
# disease = "breast"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
if(disease == "breast"){
  chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
}else if(disease == "lung"){
  chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
}else{ # DLBCL
  chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
}

save_bs_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/"

vis_rawmat <- "counts"

for(i in 1:nsamples){
  print(i)
  
  # chrom_pt_tumor <- chrom_classes[grepl(paste0("Tu_L", i), chrom_classes)]
  # chrom <- chrom[, chrom$Harmonised_Level4 %in% c(chrom_immune_include, chrom_pt_tumor)]
  
  ## Get Spot gene expr by sample -----------------------------------------
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  
  save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
  
  ## Decon -------------------------------------------------------------------
  # Ref: prep count matrix and cell anno
  cell_types <- data.frame(
    cellID = colnames(chrom),
    cellType = as.factor(chrom@meta.data[["level1_5"]]), 
    sampleInfo = "sample1")
  rownames(cell_types) <- colnames(chrom)
  # chrom_mat <- chrom@assays$SoupX@counts # decont chromium 
  chrom_mat <- chrom@assays$RNA@counts # raw chromium: use raw chrom for decont
  
  # sce: prep coords and count matrix
  coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
  counts_sce <- assay(sce, vis_rawmat)
  
  CARD_obj = createCARDObject(
    sc_count = chrom_mat,
    sc_meta = cell_types,
    spatial_count = counts_sce,
    spatial_location = coords,
    ct.varname = "cellType",
    ct.select = unique(cell_types$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5) 
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
  
  results = data.frame(CARD_obj@Proportion_CARD)
  write.csv(results, paste0(save_bs_path_sample, save_names[i], "_spot_level1_5_decon.csv"))
  
  
  # Visualize each cell type ------------------------------------------------
  # ## select the cell type that we are interested
  # ct.visualize = levels(cell_types$cellType)
  # ## visualize the spatial distribution of the cell type proportion
  # p2 <- CARD.visualize.prop(
  #   proportion = CARD_obj@Proportion_CARD,        
  #   spatial_location = CARD_obj@spatial_location, 
  #   ct.visualize = ct.visualize,                 ### selected cell types to visualize
  #   colors = c("darkblue","orange","red"),       ### if not provide, we will use the default colors
  #   NumCols = 4,                                 ### number of columns in the figure panel
  #   pointSize = 0.2) + ggtitle(save_names[i])    ### point size in ggplot2 scatterplot  
  # # print(p2)
  # 
  # each_plot_space <- 2
  # plot_ratio <- (max(sce$pxl_col_in_fullres) - min(sce$pxl_col_in_fullres)) /
  #   (max(sce$pxl_row_in_fullres) - min(sce$pxl_row_in_fullres))
  # plot_width <- plot_ratio * each_plot_space * 4 
  # plot_height <- each_plot_space * ceiling(length(colnames(results))/4)
  # 
  # pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_all_spot_level1_5_decon_bycelltype.pdf")),
  #     height = plot_height, width = plot_width)
  # print(p2)
  # dev.off()
  
  
  # # Decon label ---------------------------------------------------------
  # # 1st dominant cell type
  # results$max <- colnames(results)[max.col(results, ties.method="first")]
  # if(nrow(results) < ncol(sce)){
  #   sce$barcode <- colnames(sce)
  #   results$barcode <- rownames(results)
  #   
  #   CD <- as.data.frame(colData(sce))
  #   CD <- CD %>% 
  #     left_join(results %>% select(max, barcode)) %>%
  #     rename(level1_5_decon_max = max)
  #   rownames(CD) <- CD$barcode; CD$barcode <- NULL
  #   
  #   colData(sce) <- as(CD, "DFrame")
  # }else{
  #   sce[["level1_5_decon_max"]] <- as.factor(results$max)
  # }
  
  # # Save back the object ------------------------------------------------
  # saveRDS(sce, readpath)
}
  

# clusterPlot_by_cluster <- function(sce_obj, annotate = "spatial.cluster", title = " "){
#   spatialCoords(sce_obj) <- NULL
#   vec <- as.factor(sce_obj[[annotate]])
#   for (j in 1:nlevels(vec)){
#     sce_obj[[levels(vec)[j]]] <- ifelse(sce_obj[[annotate]] == levels(vec)[j], TRUE, FALSE)
#   }
#   clus <- levels(vec)
#   feat.plots <- purrr::map(clus, function(x) plotSpotQC(sce_obj, annotate = x, type = "spots",
#                                                         x_coord = "pxl_col_in_fullres", y_coord = "pxl_row_in_fullres", y_reverse = FALSE) +
#                              theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
#                              ggtitle(x))
#   clusters_plot <- patchwork::wrap_plots(feat.plots, ncol=4) +
#     plot_annotation(title = title) &
#     theme(plot.title = element_text(hjust = 0.5))
# 
#   return(clusters_plot)
# }
# 
# 
# clusterPlot_by_cluster(sce_obj = sce, annotate = "level1_5_decon_max", title = "CARD level 1.5 spot deconvolution")








