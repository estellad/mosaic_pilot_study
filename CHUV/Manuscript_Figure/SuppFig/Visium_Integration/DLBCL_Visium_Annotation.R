library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SCpubr)
library(patchwork)

disease = "dlbcl"

# After spotclean -----------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
savepath = savepath.spotclean

# Better plots for SCT before and after spotclean -----------------------
seu <- readRDS(paste0(savepath, save_rds_name))

DimPlot(seu, label = TRUE, group.by = "Level4_decon_max") + facet_wrap(~Level4_decon_max) + 
  theme(legend.position = "none")
DimPlot(seu, label = TRUE, group.by = "seurat_clusters") + facet_wrap(~seurat_clusters) + 
  theme(legend.position = "none")
p <- DimPlot(seu, label = TRUE, group.by = "Region") + facet_wrap(~Region) + 
  theme(legend.position = "none")
figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_UMAP_facet_patho.pdf"), width = 10, height = 10)
print(p)
dev.off()



# Plasma cells
FeaturePlot(seu, features = c("IGKC", "IGHG3", "IGHA", "JCHAIN", "IGHG1"))
# Epithelium
FeaturePlot(seu, features = c("MUCL3", "LCN2", "TFF1", "MUC5AC"))
# Stroma
FeaturePlot(seu, features = c("THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1"), ncol = 3)
# Vessels
FeaturePlot(seu, features = c("ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4"), ncol = 4)

# Marcophage 
FeaturePlot(seu, features = c("CD68", "CD163", "CD14", "ITGAM", "MRC1", "MSR1", "SIRPA", "ARG1"))
FeaturePlot(seu, features = rownames(seu)[grepl("^MMP", rownames(seu))])
FeaturePlot(seu, features = rownames(seu)[grepl("^TGFB", rownames(seu))])


table# Tu_D1_LMO2
FeaturePlot(seu, features = c("LMO2"))
# Tu_D1_SMIM14
FeaturePlot(seu, features = c("SMIM14"))

# # Tu_D2_mito
# FeaturePlot(seu, features = c("mito?"))
FeaturePlot(seu, features = rownames(seu)[grepl("^MT-", rownames(seu))])

# check marker gene from heatmap for clus 14 and 16
FeaturePlot(seu, features = c("SLC40A1"))

FeaturePlot(seu, features = c("MMP12"))

# Tu_D3_FAM3C
FeaturePlot(seu, features = c("FAM3C"))
# Tu_D4_BCL7A
FeaturePlot(seu, features = c("BCL7A"))
# Tu_D5_CCL22
FeaturePlot(seu, features = c("CCL22"))
# Tu_D6_BCL2
FeaturePlot(seu, features = c("BCL2"))

seu_sub <- seu[, seu$seurat_clusters == "18"]
table(seu_sub$sample_id, seu_sub$Level4_decon_max)

# SuppFig (cluster vs marker genes dotplot) -------------------------------
# Key tumor and annotation markers DLBCL
tu_dlbcl <- c(# "CD19", "CD79A", "BCL6", "IRF4", "MYC", "BCL2", "MKI67", "TP53", # Tumor
              "IGKC", "IGHG3", "IGHA", "JCHAIN", "IGHG1",                      # Plasma
              "MUCL3", "LCN2", "TFF1", "MUC5AC",                               # Epithelium
              "THBS1", "POSTN", "COL3A1", "COL1A2", "COL1A1",                  # Stroma
              "ACKR1", "MEOX2", "CCL14", "TSPAN7", "PLA1A", "VWF", "GIMAP4"    # Vessels
              )

p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
                        axis.text.x.angle = 0,
                        flip = TRUE# , legend.position = "right"
                        ) + 
  theme(axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_x_discrete(position = "top")
p

# SuppFig (tumor cluster vs pt-specific tumor marker genes dotplot) -------
p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
                        axis.text.x.angle = 0,
                        flip = TRUE# , legend.position = "right"
) + 
  theme(axis.ticks = element_blank(),
        axis.line = element_blank()) + 
  scale_x_discrete(position = "top")
p

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Pt_Clus_Heatmap.pdf"), width = 10, height = 5)
print(p)
dev.off()


# Heatmap of tumor clus versus sample id ----------------------------------
seu_sub <- seu 
seu_sub$sample_id <- paste0(substr(seu$sample_id, 1, 1), substr(seu$sample_id, 7, 7))
seu_sub$seurat_clusters <- droplevels(seu_sub$seurat_clusters)
Idents(seu_sub) <- seu_sub$seurat_clusters
df <- as.data.frame(table(seu_sub$seurat_clusters, seu_sub$sample_id))
df2 <- df %>%
  group_by(Var1) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(cluster = df2$Var1, sample = df2$Var2, per = df2$per)
data_frame$sample <- factor(data_frame$sample, levels = rev(levels(factor(data_frame$sample))))

# Create the heatmap pt per cluster
p <- ggplot(data_frame, aes(x = cluster, y = sample, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_fill_viridis_c(option = "mako", direction = -1) + 
  # scale_fill_gradient(low = "#faf9e4", high = "#12111d") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(fill = "Per-cluster Percentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Sample") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(23, "cm"), barheight = unit(0.5, "cm")))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Pt_Clus_Heatmap.pdf"), width = 10, height = 5)
print(p)
dev.off()


# Heatmap of tumor clus versus Region ----------------------------------
seu_sub <- seu 
df <- as.data.frame(table(seu_sub$seurat_clusters, seu_sub$Region))
df2 <- df %>%
  group_by(Var1) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(cluster = df2$Var1, Region = df2$Var2, per = df2$per)
data_frame$Region <- factor(data_frame$Region, levels = rev(levels(factor(data_frame$Region))))

# Create the heatmap pt per cluster
p <- ggplot(data_frame, aes(x = cluster, y = Region, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_fill_viridis_c(option = "rocket", direction = -1) + 
  # scale_fill_gradient(low = "#faf9e4", high = "#12111d") +
  theme_minimal() +
  theme(legend.position = "bottom") + 
  labs(fill = "Per-cluster Percentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Pathology Annotation") + 
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, 
                               barwidth = unit(21, "cm"), barheight = unit(0.5, "cm")))

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_DLBCL_Annotation"
pdf(file.path(figpath, "/DLBCL_Vis_Patho_Clus_Heatmap.pdf"), width = 10, height = 5)
print(p)
dev.off()


# Chrom pt spec marker in Visium  -------------------------------------------
plot_UMAP_agg_markers <- function(seu, genelist, genelistname = "Tu_D1"){
  seu <- AddModuleScore(object = seu, features = genelist, name = genelistname)
  p <- FeaturePlot(object = seu, features = paste0(genelistname, "1"), label = TRUE) + 
    scale_color_viridis_c(option = "magma", direction = -1) + 
    ggtitle(paste0(genelistname))
  return(p)
}

Tu_D1 <- c("NME2", "FCRL1", "LMO2", "IGKC", "MT-CO3", "TMSB4X", "GRHPR",
           "MPEG1", "CD1C", "CD52", "LTB", "SMIM14", "TNFRSF13C", "CD79A",
           "CD22", "FCRL2", "FCRL5")
Tu_D2 <- c("MT-CYB", "NIBAN3", "SPIB", "MT-ATP6", "HIST1H1C", "IGHM",
           "MT-ND4", "TCL1A", "MT-ND4L", "HIST1H1R") 
Tu_D3 <- c("CD83", "SWAP70", "SEL1L3", "FCRL3", "ACTB", "IL4I1", "FAM3C",
           "ACTG2", "CD74", "MS4A1", "MT-CYB", "CNN2", "LRMP", "LCP1")
Tu_D4 <- c("GGA2", "HELLS", "RASGRP2", "NKX6-3", "DTX1", "BCL7A", "MAT2A",
           "ARGLU1", "ALOX5", "AKNA", "ARHGEF1", "PNN", "TMC8")
Tu_D5 <- c("TSPAN33", "SPATC1", "NFKB2", "BCL2L1", "OGI", "CCL17",
           "CD40", "DUSP2", "ITPKB", "CCL22")
Tu_D6 <- c("POU2F2", "CCDC88A", "LENG8", "KLHL6", "TNFRSF13B", "BCL2",
           "NFATC1", "MT-ND6", "NIBAN3", "FCRL2")

markers <- lapply(list(Tu_D1, Tu_D2, Tu_D3, Tu_D4, Tu_D5, Tu_D6), function(x) x[x %in% rownames(seu)])
names(markers) <- c("Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
UMAP_plots <- purrr::map2(markers, names(markers), function(genelist, genelistname) plot_UMAP_agg_markers(seu, genelist, genelistname))
patchwork::wrap_plots(UMAP_plots, ncol=6) + plot_annotation(title = "Aggregated Markers")





FeaturePlot(seu, features = Tu_D1)

tu_dlbcl <- c("NME2", "FCRL1", "LMO2", "IGKC", "MT-CO3", "TMSB4X", "GRHPR",
              "MPEG1", "CD1C", "CD52", "LTB", "SMIM14", "TNFRSF13C", "CD79A",
              "CD22", "FCRL2", "FCRL5",                                    # D1 clus 2, 5, 13, 19:
              "MT-CYB", "NIBAN3", "SPIB", "MT-ATP6", "HIST1H1C", "IGHM",
              "MT-ND4", "TCL1A", "MT-ND4L", "HIST1H1R",                    # D2 clus 0:
              "CD83", "SWAP70", "SEL1L3", "FCRL3", "ACTB", "IL4I1", "FAM3C",
              "ACTG2", "CD74", "MS4A1", "MT-CYB", "CNN2", "LRMP", "LCP1",  # D3 clus 1, 6, 26:
              "GGA2", "HELLS", "RASGRP2", "NKX6-3", "DTX1", "BCL7A", "MAT2A",
              "ARGLU1", "ALOX5", "AKNA", "ARHGEF1", "PNN", "TMC8",         # D4 clus 9&22:
              "TSPAN33", "SPATC1", "NFKB2", "BCL2L1", "OGI", "CCL17",
              "CD40", "DUSP2", "ITPKB", "CCL22",                           # D5 clus 17:
              "POU2F2", "CCDC88A", "LENG8", "KLHL6", "TNFRSF13B", "BCL2",
              "NFATC1", "MT-ND6", "NIBAN3", "FCRL2"                        # D6 clus 12:
)


# p <- SCpubr::do_DotPlot(sample = seu, features = tu_dlbcl, 
#                         axis.text.x.angle = 45,
#                         flip = FALSE, legend.position = "right"
# ) + 
#   theme(axis.ticks = element_blank(),
#         axis.line = element_blank()) + 
#   scale_fill_viridis_c(option = "magma", direction = -1) + 
#   scale_x_discrete(position = "top")
# p

## C'est impossible

# -------------------------------------------------------------------------
# wide format
test <- table(seu$Level4_decon_max, seu$seurat_clusters)
# long format with counts
test2 <- as.data.frame(test)

# # Heatmap (decon x clus) --------------------------------------------------
# data_frame <- data.frame(cluster = test2$Var2, decon = test2$Var1, counts = test2$Freq)
# data_frame$decon <- factor(data_frame$decon, levels = rev(levels(factor(data_frame$decon))))
# 
# # Create the heatmap
# ggplot(data_frame, aes(x = cluster, y = decon, fill = counts)) + 
#   geom_tile() + 
#   scale_fill_gradient(low = "white", high = "red") +
#   theme_minimal() +
#   labs(title = "", x = "Level 4 Deconvolution Majority Vote", y = "Integrated DLBCL Visium Seurat Clusters")


# -------------------------------------------------------------------------
# long format with percentage per-clus
test3 <- test2 %>%
  group_by(Var2) %>%
  mutate(countT= sum(Freq)) %>%
  mutate(per=round(100*Freq/countT,0))

data_frame <- data.frame(cluster = test3$Var2, decon = test3$Var1, per = test3$per)
data_frame$decon <- factor(data_frame$decon, levels = rev(levels(factor(data_frame$decon))))

# Create the heatmap
ggplot(data_frame, aes(x = cluster, y = decon, fill = per)) + 
  geom_tile() + 
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(fill = "Per-cluster\nPercentage", title = "", 
       x = "Integrated DLBCL Visium Seurat Clusters", y = "Level 4 Deconvolution Majority Vote") 


# wide format
library(tidyverse)
test5 <- test3 %>%
  select(Var1, Var2, per) %>%
  spread(Var2, per)



