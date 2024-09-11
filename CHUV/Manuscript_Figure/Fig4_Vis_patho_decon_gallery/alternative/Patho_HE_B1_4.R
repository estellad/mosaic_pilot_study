library(ggspavis)
library(ggplot2)

disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Visium/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/color_palette.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/Manuscript_Figure/Patho_hires_HE_helper.R")
sample = "B1_4"
if(disease == "lung"){
  patho_color <- lung_patho_color
}else if(disease == "breast"){
  patho_color <- breast_patho_color
}else{
  patho_color <- dlbcl_patho_color
}

sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(sce_path, paste0(sample, "_qcd.rds")))

datapath_rerun <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/VisiumRerun/", sample)
coord_new <- read.csv(file.path(datapath_rerun, "tissue_positions.csv"))

CD <- data.frame(colData(sce))

CD_new <- coord_new %>% 
  mutate(Barcode = barcode) %>% 
  select(Barcode, pxl_col_in_fullres, pxl_row_in_fullres)

CD <- CD %>%
  left_join(CD_new, by = "Barcode")
rownames(CD) <- CD$Barcode

colData(sce) <- as(CD, "DFrame")


sce_ <- sce
colors <- patho_color[names(table(sce_$Region))]

scalef <- rjson::fromJSON(file = file.path(datapath_rerun, "scalefactors_json.json"))

# imgData(sce_)

sce_ <- rmvImg(sce_)
sce_ <- addImg(sce_, 
               sample_id = "sample01", # sample, 
               image_id = "hires",
               imageSource = file.path(datapath_rerun, "tissue_hires_image.png"), 
               scaleFactor = scalef$tissue_hires_scalef, 
               load = TRUE)
# img <- imgRaster(sce_,
#                  sample_id = "sample01", # sample,
#                  image_id = "hires")
# plot(img)
# imgData(sce_)

sce_$Region <- as.factor(sce_$Region)

########################### Patho HE full slide #########################
p2 <- plotVisium_new(
  sce_, 
  annotate = "Region", 
  palette = patho_color[names(table(sce_$Region))], 
  zoom = TRUE, 
  guide.pt.size = 3,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2) + 
  theme(legend.key.size = unit(0.8, 'lines')) 


# Image size and save ---------
image_height <- (range(spatialCoords(sce_)[, 2])[2] - range(spatialCoords(sce_)[, 2])[1])
image_width <- (range(spatialCoords(sce_)[, 1])[2] - range(spatialCoords(sce_)[, 1])[1])

plot_width <- 10.5
plot_height <- plot_width * (image_height/image_width) - 3

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_Patho_HE"
pdf(file = file.path(fig_path, paste0(sample, "_patho_HE.pdf")),
    width = plot_width,
    height = plot_height)
print(p2)
dev.off()


########################################################################
# With coordinates to see where to zoom
p3 <- plotVisium_new(
  sce_,
  annotate = "Region",
  palette = patho_color[names(table(sce_$Region))],
  zoom = TRUE,
  spots = FALSE,
  show_axis = TRUE,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2
)

########################################################################
# Just HE but zoomed to B cell
# Too blurry, but as reference for python crop ------------------
p3_HE_Bcell <- plotVisium_new(
  sce_,
  annotate = "Region",
  palette = patho_color[names(table(sce_$Region))],
  zoom = TRUE,
  # show_axis = TRUE,
  spot = FALSE,
  image = TRUE,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2
) + xlim(c(500, 606)) + ylim(c(637, 790)) + theme(legend.position = "none")

plot_width <- 1.5
plot_height <- plot_width * (790-637)/(606-500)

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped/archive"
pdf(file = file.path(fig_path, paste0(sample, "_HE_Bcell_archive.pdf")),
    width = plot_width,
    height = plot_height)
print(p3_HE_Bcell)
dev.off()


# Just HE but zoomed to T cell
# Too blurry, but as reference for python crop ------------------
p3_HE_Tcell <- plotVisium_new(
  sce_,
  annotate = "Region",
  palette = patho_color[names(table(sce_$Region))],
  zoom = TRUE,
  # show_axis = TRUE,
  spot = FALSE,
  image = TRUE,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2
) + xlim(c(755, 920)) + ylim(c(1000, 1135)) + theme(legend.position = "none")

plot_width <- 2.5
plot_height <- plot_width * (1135-1000)/(920-755)

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped/archive"
pdf(file = file.path(fig_path, paste0(sample, "_HE_Tcell_archive.pdf")),
    width = plot_width,
    height = plot_height)
print(p3_HE_Tcell)
dev.off()

########################################################################
# Patho_HE - B cell -------------------------------------------
p3_Bcell <- plotVisium_new(
  sce_, 
  annotate = "Region", 
  palette = patho_color[names(table(sce_$Region))], 
  zoom = TRUE, 
  #show_axis = TRUE,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2
) + xlim(c(500, 606)) + ylim(c(637, 790)) + theme(legend.position = "none")

plot_width <- 1.5
plot_height <- plot_width * (790-637)/(606-500)

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
pdf(file = file.path(fig_path, paste0(sample, "_patho_HE_Bcell.pdf")),
    width = plot_width,
    height = plot_height)
print(p3_Bcell)
dev.off()


# Patho_HE - T cell -------------------------------------------
p3_Tcell <- plotVisium_new(
  sce_,
  annotate = "Region",
  palette = patho_color[names(table(sce_$Region))],
  zoom = TRUE,
  # show_axis = TRUE,
  pt.size = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2
) + xlim(c(755, 920)) + ylim(c(1000, 1135)) + theme(legend.position = "none")

plot_width <- 2.5
plot_height <- plot_width * (1135-1000)/(920-755)

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig4/B1_cropped"
pdf(file = file.path(fig_path, paste0(sample, "_patho_HE_Tcell.pdf")),
    width = plot_width,
    height = plot_height)
print(p3_Tcell)
dev.off()


# -------------------------------------------------------------------------
## Scatterpie


