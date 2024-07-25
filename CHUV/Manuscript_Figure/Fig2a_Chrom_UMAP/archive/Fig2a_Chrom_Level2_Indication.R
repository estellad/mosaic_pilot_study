library(Seurat)
library(dplyr)

# Breast and lung -------------------------------------------------------
chrom_brlu_palette <- data.frame(
  rbind(
    c("B", "#EEEE00"), 
    c("Epithelia", "#BC8F8F"), 
    c("Fibro_muscle", "#388E8E"), 
    c("Granulocyte", "#FFC125"), 
    c("Mast_cell", "#FFD700"), 
    c("Myeloid", "#9A32CD"),
    c("T_NK", "#4169E1"),
    c("Vessel", "#FF7256"), 
    c("Tu_B1", "#00C78C"), 
    c("Tu_B3", "#66CDAA"), 
    c("Tu_B4", "#3CB371"), 
    c("Tu_L1", "#EEB4B4"), 
    c("Tu_L2", "#FFA54F"), 
    c("Tu_L3", "#CD853F"), 
    c("Tu_L4", "#8B4513")
    )
  )
color <- chrom_brlu_palette$X2
names(color) <- chrom_brlu_palette$X1

# DLBCL ----------------------------------------------------------------
chrom_dlbcl_palette <- data.frame(
  rbind(
    c("B", "#EEEE00"),
    c("Epithelia", "#BC8F8F"),
    c("Fibro_muscle", "#388E8E"),
    c("Myeloid", "#9A32CD"),
    c("T_NK", "#4169E1"),
    c("Vessel", "#FF7256"),
    c("Tu_D1", "#FF8C00"),
    c("Tu_D2", "#EEEE00"),
    c("Tu_D3", "#FFD700"),
    c("Tu_D4", "#A2CD5A"),
    c("Tu_D5", "#00EE76"),
    c("Tu_D6", "#ADFF2F")
    )
  )

color <- chrom_dlbcl_palette$X2
names(color) <- chrom_dlbcl_palette$X1

#########################################################################
# -------------------------------------------------------------------------
chrom_prev_brlu <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/final_owkin_annot.rds")
Idents(chrom_prev_brlu) <- factor(chrom_prev_brlu$Level2, levels = chrom_brlu_palette$X1)
p <- DimPlot(chrom_prev_brlu, cols = color)
  
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_BrLu_Level2.pdf"), width = 9, height = 6)
print(p)
dev.off()


# -------------------------------------------------------------------------
color <- c("#00c78cff", "#ffa54fff")
names(color) <- c("Breast", "Lung")
chrom_prev_brlu$indication <- 
  case_when(chrom_prev_brlu$sample_id %in% c("B3", "B2", "B1", "B4", "B3_rep", "B4_rep") ~ "Breast",
            chrom_prev_brlu$sample_id %in% c("L3", "L2", "L1", "L4") ~ "Lung")
Idents(chrom_prev_brlu) <- factor(chrom_prev_brlu$indication, levels = c("Breast", "Lung"))
p <- DimPlot(chrom_prev_brlu, cols = color)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_BrLu_Indication.pdf"), width = 9, height = 6)
print(p)
dev.off()


# -------------------------------------------------------------------------
chrom_dlbcl <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")) # 16066 39713
chrom_dlbcl$Level2 <- ifelse(chrom_dlbcl$Level2 == "Fibro_Muscle", "Fibro_muscle", chrom_dlbcl$Level2)
Idents(chrom_dlbcl) <- factor(chrom_dlbcl$Level2, levels = chrom_dlbcl_palette$X1)
p <- DimPlot(chrom_dlbcl, cols = color)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig2"
pdf(file.path(figpath, "/Chrom_UMAP_DLBCL_Level2.pdf"), width = 9, height = 6)
print(p)
dev.off()



















