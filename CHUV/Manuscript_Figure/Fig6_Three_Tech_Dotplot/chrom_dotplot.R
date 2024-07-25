library(Seurat)
library(ggplot2)
# chrom <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/DLBCL/dlbcl_final_owkin_annot_pub.rds")
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

## Level4 tumor subtypes
chrom$level1_5 <- chrom$Harmonised_Level4
chrom$level1_5 <- ifelse(chrom$Level1 == "B", "B", 
                         ifelse(chrom$Level1 == "Epithelia", "Epithelia", 
                                ifelse(chrom$Level1 == "Myeloid", "Myeloid", 
                                       ifelse(chrom$Level1 == "Stroma", "Stroma", 
                                              ifelse(chrom$Level1 == "T_NK", "T_NK", 
                                                     chrom$Harmonised_Level4)))))
chrom$level1_5 <- as.factor(chrom$level1_5)
Idents(chrom) <- chrom$level1_5

## Combined tumor subtypes
chrom$level1_5 <- chrom$Level2
chrom$level1_5 <- ifelse(chrom$level1_5 %in% c("Fibro_Muscle", "Vessel"), "Stroma", chrom$level1_5)
chrom$level1_5 <- as.factor(chrom$level1_5)
Idents(chrom) <- chrom$level1_5

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

# DLBCLnChromium_Marker_Gene_List <- DLBCLnChromium_Marker_Gene_List[order(DLBCLnChromium_Marker_Gene_List, decreasing = TRUE)]

chrom <- chrom[rownames(chrom) %in% DLBCLnChromium_Marker_Gene_List, ]

p <- DotPlot(
  chrom,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))

# plot_title = "Chrom_DLBLC_dotplot.pdf"
# pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig4_drug_target_dlbcl_geo/", plot_title),
#     width = 5,
#     height = 20)
# print(p)
# dev.off()