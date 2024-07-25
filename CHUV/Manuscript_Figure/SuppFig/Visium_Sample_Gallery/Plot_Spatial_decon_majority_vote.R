# disease = "breast"
disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

i = 4
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))


image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])

# # if need to coord_flip an image
# image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
# image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) -2

level4_decon_cols_max <- plyr::mapvalues(levels(sce[["Level4_decon_max"]]), from = level4_cellnames, to = level4_cellcolors)

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Manuscript_Figure/SuppFig/Plot_Spots_Gallery/helpers.R")

sce_plot <- sce[, !is.na(sce$Level4_decon_max)]
p2 <- plotSpots_deconmax(sce_plot, 
                         annotate = "Level4_decon_max",
                         pt.size = 7.5,
                         palette = level4_decon_cols_max) + 
# coord_flip() #+ 
scale_y_reverse() #+
# scale_x_reverse()

save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_spot_level4_decon_bycelltype_majorityvote_.pdf")),
    height = plot_height, width = plot_width)
print(p2)
dev.off()



