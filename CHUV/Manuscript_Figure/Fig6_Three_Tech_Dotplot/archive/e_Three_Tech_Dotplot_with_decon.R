library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)
library(patchwork)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig6/"

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

CT_order <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else",  "Macrophage", "T cells", "Tumor") # Fig 6

# Chromium ----------------------------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

## Separate tumor subtypes
Idents(chrom) <- chrom$level1_5_immune_tumor

# chrom <- chrom[rownames(chrom) %in% DLBCLnChromium_Marker_Gene_List, ]

p_chrom <- DotPlot(
  chrom,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "none"
  )


# Visium ------------------------------------------------------------------
vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis_small <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean_small_fig6e.rds"))

# vis_small <- subset(vis, features = DLBCLnChromium_Marker_Gene_List)
# vis[rownames(vis) %in% DLBCLnChromium_Marker_Gene_List, ]

p_vis <- DotPlot(
  vis_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank()
        ,legend.position = "none"
  )


# GeoMx ----------------------------------------------------------------
geo <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/dlbcl_seu_ruv.rds")
# head(Idents(geo)) # Levels: T cells Macrophage Other Tu_D1 Tu_D2 Tu_D3 Tu_D4 Tu_D5 Tu_D6
Idents(geo) <- factor(geo$clusters, levels = c("Other", "Macrophage", "T cells", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6"))

geo <- geo[rownames(geo) %in% DLBCLnChromium_Marker_Gene_List, ]
geo_small <- geo[, geo$consensus]

p_geo <- DotPlot(
  geo_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")

p_geo_ <- p_geo + scale_size(range = c(1, 8))

# Stitch plots ------------------------------------------------------------
p_final <- p_chrom + theme(plot.margin = margin(0, 4, 0, 0.5)) + 
  p_vis + theme(plot.margin = margin(0, 4, 0, 4)) +
  p_geo_ + theme(plot.margin = margin(0, 0.5, 0, 4)) +
  plot_layout(widths = c(1.20, 1.20, 0.9))


# -------------------------------------------------------------------------
plot_title = "e_Three_Tech_Cluster_Dotplot.pdf"
pdf(file = file.path(figpath, plot_title),
    width = 10,
    height = 6)
print(p_final)
dev.off()