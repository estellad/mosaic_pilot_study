library(magick)
library(ggplot2)
library(dplyr)
library(scatterpie)

disease = "breast"
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
sample = "B1_4"
if(disease == "lung"){patho_color <- lung_patho_color}else{patho_color <- breast_patho_color}

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

scalef <- rjson::fromJSON(file = file.path(datapath_rerun, "scalefactors_json.json"))
image <- image_read(file.path(datapath_rerun, "tissue_hires_image.png"))
info <- image_info(image)

sce$pxl_col_in_hires <- sce$pxl_col_in_fullres * scalef$tissue_hires_scalef
sce$pxl_row_in_hires <- sce$pxl_row_in_fullres * scalef$tissue_hires_scalef

sce$Region <- as.factor(sce$Region)


# -------------------------------------------------------------------------
p0 <- ggplot(
  data.frame(x = 0, y = 0),
  aes(x, y)
) +
  geom_blank() +
  coord_fixed(
    expand = FALSE,
    xlim = c(min(sce$pxl_col_in_hires) - info$width * 0.05, max(sce$pxl_col_in_hires) + info$width * 0.05),
    ylim = c(min(sce$pxl_row_in_hires) - info$height * 0.05, max(sce$pxl_row_in_hires) + info$height * 0.05)
  ) +
  annotation_raster(
    image_flip(image),
    0,
    info$width,
    0,
    info$height
  ) +
  theme_bw()

res_CARD <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/",
                            foldername, "/",
                            sample, "/",
                            sample, "_spot_Level4_decon.csv"), 
                     row.names = 1)
location <- data.frame(x = sce$pxl_col_in_hires,
                       y = sce$pxl_row_in_hires) 
rownames(location) <- colnames(sce)
colnames(location) <- c("x", "y")

overlap_barcodes <- intersect(rownames(res_CARD), rownames(location))
location <- location[overlap_barcodes, ]
res_CARD <- res_CARD[overlap_barcodes, ]
data = cbind(res_CARD, location)
ct.select = colnames(res_CARD)

colors = level4_cellcolors[names(level4_cellcolors) %in% sort(colnames(res_CARD))]
radius <- scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 1.7 

# -------------------------------------------------------------------------
data <- data %>%
  select(B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1, 
         T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg,
         x, y) %>%
  mutate(B_cells = B_cell + B_plasma_IGHA1 + B_plasma_IGHG1 + B_plasma_IGHG3 + B_plasma_IGHM + B_plasma_IGKC + B_plasma_IGLC1,
         T_cells = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg)

# -------------------------------------------------------------------------
p <- p0 + 
  ggforce::geom_circle(
    data = data,
    aes(
      x0 = x,
      y0 = y,
      r = scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 2,
      fill =  T_cells, # B_cells, 
      colour = T_cells  # B_cells
    ),
    inherit.aes = FALSE
  ) +
  labs(x = NULL, y = NULL) + 
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15)) +
  viridis::scale_fill_viridis(option="A") + 
  viridis::scale_color_viridis(option="A")
p


image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures"
pdf(file = file.path(paste0(fig_path, "/Fig4b_deconHE"), 
                     paste0(sample, "_decon_HE_T_cells.pdf")),
                     # paste0(sample, "_decon_HE_B_cells.pdf")),
    width = plot_width,
    height = plot_height)
print(p)
dev.off()


