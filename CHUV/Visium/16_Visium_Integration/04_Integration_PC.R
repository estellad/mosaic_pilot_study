library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(Seurat)
library(SingleCellExperiment)
library(standR)
library(cluster)
library(mclustcomp)


disease = "lung"
## LogNorm ##############################################################
## Before SpotClean -----------------------------------------------------
datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
savepath.filtered = paste0(datapath.filtered, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormnoSpotClean_filtered.rds")
savepath = savepath.filtered
pipeline = 1
seu <- readRDS(paste0(savepath, save_rds_name))
assign(paste0("seu_", pipeline), seu)


## SCT #################################################################
# After spotclean ------------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean
pipeline = 4

seu <- readRDS(paste0(savepath, save_rds_name))
seu$CellType <- case_when(seu$Level4_decon_max %in% c("B_cell", "B_plasma_IGHA1", "B_plasma_IGHG1", "B_plasma_IGHG3", "B_plasma_IGHM", "B_plasma_IGKC", "B_plasma_IGLC1") ~ "B cells",
                          seu$Level4_decon_max %in% c("Endothelia_vascular", "Fibroblast", "Muscle_smooth", "Pericyte") ~ "Stroma",
                          seu$Level4_decon_max %in% c("Epi_lung") ~ "Epithelia",
                          seu$Level4_decon_max %in% c("DC_1", "DC_2", "DC_activated", "DC_pc", "Granulocyte", "Mast_cell", "Macrophage", "Monocyte") ~ "Myeloid",
                          seu$Level4_decon_max %in% c("T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing", "NK") ~ "TNK",
                          grepl("Tu_L1", seu$Level4_decon_max) ~ "Tu_L1",
                          grepl("Tu_L2", seu$Level4_decon_max) ~ "Tu_L2",
                          grepl("Tu_L3", seu$Level4_decon_max) ~ "Tu_L3",
                          grepl("Tu_L4", seu$Level4_decon_max) ~ "Tu_L4",
                          TRUE ~ as.character(seu$Level4_decon_max))

DimPlot(seu, reduction = "pca", group.by = "CellType") | DimPlot(seu, reduction = "umap", group.by = "CellType")
seu_nomix <- seu[, seu$CellType != "Mix"]
DimPlot(seu_nomix, reduction = "pca", group.by = "CellType") | DimPlot(seu_nomix, reduction = "umap", group.by = "CellType")

DimPlot(seu, reduction = "pca", group.by = "sample_id", split.by = "sample_id")


