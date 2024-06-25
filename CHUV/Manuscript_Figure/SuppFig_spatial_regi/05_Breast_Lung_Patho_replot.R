disease = "lung"
disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

i = 5
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))


image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])

# # if need to coord_flip an image
# image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
# image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) -2
patho_color <- get(paste0(disease, "_patho_color"))

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/08_Baye_New_Pipeline_Decon/helpers.R")

sce_plot <- sce[, !is.na(sce$Region)]
p2 <- plotSpots_deconmax(sce_plot, 
                         annotate = "Region",
                         pt.size = 7.3,
                         palette = patho_color[names(table(sce$Region))]) + 
# coord_flip() # +
scale_y_reverse() # +
# scale_x_reverse()

save_bs_path_sample <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/SuppFig_regi_spatial/With_legend/Pathology"
pdf(file.path(save_bs_path_sample, paste0(save_names[i], "_patho.pdf")),
    height = plot_height, width = plot_width)
print(p2)
dev.off()



