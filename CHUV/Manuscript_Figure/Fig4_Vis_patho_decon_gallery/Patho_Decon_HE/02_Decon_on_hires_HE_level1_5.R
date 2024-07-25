library(magick)
library(ggplot2)
library(dplyr)
library(scatterpie)

disease = "breast"
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")
sample = "B1_2"
if(disease == "lung"){patho_color <- lung_patho_color}else{patho_color <- breast_patho_color}

save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

# Spot --------------------------------------------------------------------
sce_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
sce <- readRDS(file.path(sce_path, paste0(sample, "_qcd.rds")))

# # Subspot --------------------------------------------------------------------
# save_bs_path_sample <- paste0(save_bs_path, foldername, "/", sample, "/")
# readpath = file.path(save_bs_path_sample, paste0(sample, "_baye_clustered_all_enhanced_expr.rds")) # all genes
# sce <- readRDS(readpath) 

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

# CARD or RCTD
# -------------------------------------------------------------------------
res_CARD <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/",
                     foldername, "/",
                     sample, "/",
                     sample, 
                     "_spot_level1_5_decon.csv"), # CARD
                     # "_subspot_level1_5_decon.csv"), # CARD subspot
                     # "_spot_level1_5_decon_RCTD.csv"), # RCTD
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

colors = level1_5_cellcolors[names(level1_5_cellcolors) %in% sort(colnames(res_CARD))]
radius <- scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 1.7

p <- p0 + 
  geom_scatterpie(aes(x = x, y = y, r = radius), 
                  data = data, cols = ct.select, color = NA) +
  scale_fill_manual(values =  colors) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 16,face="bold"),
        legend.position="right")+
  guides(fill=guide_legend(title="Cell Type", ncol = 1))

if(sample == "B1_4"){
  image_height <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
  image_width <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
}else{
  image_height <- (range(spatialCoords(sce)[, 2])[2] - range(spatialCoords(sce)[, 2])[1])
  image_width <- (range(spatialCoords(sce)[, 1])[2] - range(spatialCoords(sce)[, 1])[1])
}

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures"
pdf(file = file.path(paste0(fig_path, "/Fig4b_deconHE"), 
                     paste0(sample, "_decon_HE_level1_5.pdf")), # CARD
                                                                # did not finish attempt of CARD subspot
                     # paste0(sample, "_decon_HE_level1_5_RCTD.pdf")), # RCTD
    width = plot_width,
    height = plot_height)
print(p)
dev.off()

