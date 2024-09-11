# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

# i = 2
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))

# -------------------------------------------------------------------------
vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
results <- read.csv(file.path(vis_decon_path, paste0(save_names[i], ".csv")), row.names = 1) 

# C2L all genes -----
if(save_names[i] == "B3_2"){
  results <- results %>%
    dplyr::rename(Fibroblast = Fibroblast_B3)
}else if(save_names[i] %in% c("B1_2", "B1_4", "B2_2")){
  results <- results %>%
    mutate(Fibroblast = Fibroblast + Fibroblast_B3) %>%
    select(-Fibroblast_B3)
}
# # -------------------
# 
# image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
# image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
# 
# # # if need to coord_flip an image
# # image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
# # image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
# 
# plot_width <- 10
# plot_height <- plot_width * (image_height/image_width) -2

# Decon label ---------------------------------------------------------
# 1st dominant cell type
results$max <- colnames(results)[max.col(results, ties.method="first")]

result_max <- NULL
for(s in 1:nrow(results)){
  results_s <- results[s, ]
  if(results_s[results_s$max] < 0.2 ){
    results_s_max = "Mix"
  }else{
    results_s_max = results_s$max
  }
  result_max <- c(result_max, results_s_max)
}
results2 <- results
results2$max <- result_max

if(nrow(results2) < ncol(sce)){
  sce$Level4_decon_max <- NULL
  sce$barcode <- colnames(sce)
  results2$barcode <- rownames(results2)
  
  CD <- as.data.frame(colData(sce))
  CD <- CD %>% 
    left_join(results2 %>% select(max, barcode)) %>%
    rename(Level4_decon_max = max) %>%
    mutate(Level4_decon_max = factor(Level4_decon_max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix")))
  rownames(CD) <- CD$barcode; CD$barcode <- NULL
  
  colData(sce) <- as(CD, "DFrame")
}else{
  results2 <- results2[colnames(sce), ]
  sce[["Level4_decon_max"]] <- factor(results2$max, levels = c(sort(unique(results2$max))[sort(unique(results2$max)) != "Mix"], "Mix"))
}



level4_decon_cols_max <- plyr::mapvalues(levels(sce[["Level4_decon_max"]]), from = level4_cellnames, to = level4_cellcolors)

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/08_Baye_New_Pipeline_Decon/helpers.R")

sce_plot <- sce[, !is.na(sce$Level4_decon_max)]
p2 <- plotSpots_deconmax(sce_plot, 
                         annotate = "Level4_decon_max",
                         pt.size = 2.5,
                         palette = level4_decon_cols_max) # + 
# coord_flip() # + 
# scale_y_reverse() +
# scale_x_reverse()
p2
saveRDS(sce, file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))

# # -------------------------------------------------------------------------
# save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
# pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_spot_level4_decon_bycelltype_majorityvote_.pdf")),
#     height = plot_height, width = plot_width)
# print(p2)
# dev.off()



