library(Seurat)
library(dplyr)

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/"
chrom_breast <- readRDS(file.path(chrompath, "chrom_breast.rds")) # 17170 10689

Idents(chrom_breast) <- factor(chrom_breast$sample_id, levels = c("B1", "B2", "B3", "B3_rep", "B4", "B4_rep"))
p <- DimPlot(chrom_breast)

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/Fig1/Fig1d"
pdf(file.path(figpath, "/Chrom_UMAP_Breast_Sample.pdf"), width = 6, height = 4.5)
print(p)
dev.off()