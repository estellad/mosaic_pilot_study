library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)
library(patchwork)

figpath_ridge <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Fig6_Ridge"

DLBCLnChromium_Marker_Gene_List <- c(
  "MS4A1", "TNFRSF13C", "CD79B", "CD37", "PSMB8", "CD19", "TYMS",
  "TUBB", "TOP2A", "POLD4", "CD47", "CD52", "BLK", "CD38", "MAP2K1",
  "CD40", "BCL2L1", "TNFRSF8", "SMO", "RARA", "TYK2", "TNFRSF10B"
)

vispath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/dlbcl/spotclean/Results/"
vis <- readRDS(file.path(vispath, "Dlbcl-merge-SCTpostSpotClean.rds"))

# Decon levels ------------------------------------------------------------
decon_result_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2_final_immune"
dlbcl_decon <- read.csv(file.path(decon_result_path, "vis_dlbcl_decon_immune_long.csv"))

dlbcl_decon_ <- dlbcl_decon %>%
  mutate(CellType = case_when(CellType == "Tumor" & Section == "D1" ~ "Tu_D1",
                              CellType == "Tumor" & Section == "D2" ~ "Tu_D2",
                              CellType == "Tumor" & Section == "D3" ~ "Tu_D3",
                              CellType == "Tumor" & Section == "D4" ~ "Tu_D4",
                              CellType == "Tumor" & Section == "D5" ~ "Tu_D5",
                              CellType == "Tumor" & Section == "D6" ~ "Tu_D6",
                              .default = as.character(CellType)
                              ))

palette_ridge <- c("#BC8F8F", "#388E8E", "#EEEE00", "#79CDCD", "#C71585", "#9A32CD", "#4169E1", "#FF8C00", "#EEEE00", "#FFD700", "#A2CD5A", "#00EE76", "#ADFF2F")
ridge_order <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else", "Macrophage", "T cells", "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6")
names(palette_ridge) <- ridge_order

# table(dlbcl_decon_$CellType)
# B cells    Epithelia   Macrophage Myeloid else           NK       Stroma      T cells        Tu_D1        Tu_D2        Tu_D3        Tu_D4        Tu_D5        Tu_D6 
# 18580        18580        18580        18580        18580        18580        18580         4704         4111         4951         1436         1729         1649 

dlbcl_decon_$CellType <- factor(dlbcl_decon_$CellType, levels = rev(ridge_order))


# Max cell type calculation -----------------------------------------------
dlbcl_decon_ <- dlbcl_decon_ %>%
  mutate(Barcode2 = paste0(Barcode, "_", substr(Section, 2, 2)))

dlbcl_decon_max <- dlbcl_decon_ %>%
  group_by(Barcode2) %>%
  arrange(desc(Fraction), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# table(dlbcl_decon_max$CellType)
# Tu_D6        Tu_D5        Tu_D4        Tu_D3        Tu_D2        Tu_D1      T cells   Macrophage Myeloid else           NK      B cells       Stroma    Epithelia 
# 158            9          383          274         3455         4179         1881          850           20           61           58         4824         2428

# A ridge plot with all spots, even max fraction of a spot is less than 50% ---- 
p <- ggplot(dlbcl_decon_max,
            aes(x = Fraction, y = CellType, fill = CellType)) +
  geom_density_ridges(scale = 4, alpha = 0.5) + 
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0)) +   
  coord_cartesian(clip = "off") + 
  xlim(c(0, 1)) +
  theme_ridges() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = palette_ridge)

sample_size_df <- data.frame(table(dlbcl_decon_max$CellType)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
                   aes(label = n, x= 1, y=ylabel, vjust = -2, hjust=-0.1),
                   position = position_stack(),
                   inherit.aes = FALSE) 

plot_title = "Vis_ridge_immune_decon_max.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
    width = 12.5,
    height = 8)
print(p)
dev.off()

# Now a ridge only if the max fraction is > 0.5 ---------------------------
dlbcl_decon_max_50 <- dlbcl_decon_max %>% filter(Fraction > 0.5)

p <- ggplot(dlbcl_decon_max_50,
            aes(x = Fraction, y = CellType, fill = CellType)) +
  geom_density_ridges(scale = 4, alpha = 0.5) + 
  scale_y_discrete(expand = c(0, 0)) +     
  scale_x_continuous(expand = c(0, 0)) +   
  coord_cartesian(clip = "off") + 
  xlim(c(0, 1)) +
  theme_ridges() + 
  theme(legend.position = "none") + 
  scale_fill_manual(values = palette_ridge)

sample_size_df <- data.frame(table(dlbcl_decon_max_50$CellType)) %>%
  dplyr::rename(n = Freq,
                ylabel = Var1) %>%
  mutate(n = paste0("n = ", n))

p <- p + geom_text(data=sample_size_df,
                   aes(label = n, x= 1, y=ylabel, vjust = -2, hjust=-0.1),
                   position = position_stack(),
                   inherit.aes = FALSE) 

plot_title = "Vis_ridge_immune_decon_max_50.pdf"
pdf(file = file.path(figpath_ridge, plot_title),
    width = 12.5,
    height = 8)
print(p)
dev.off()

# -------------------------------------------------------------------------
df_test <- dlbcl_decon_ %>% filter(Fraction > 0.5)

df_test <- df_test %>%
  mutate(Barcode2 = paste0(Barcode, "_", substr(Section, 2, 2)))

vis_small <- vis[, colnames(vis) %in% df_test$Barcode2]

CD <- data.frame(Barcode = colnames(vis_small))

decon_max_label <- data.frame(
  Barcode = df_test$Barcode2,
  decon_max = df_test$CellType)

CD <- CD %>%
  left_join(decon_max_label)

vis_small$new_annot <- CD$decon_max

Idents(vis_small) <- factor(vis_small$new_annot, levels = ridge_order)

saveRDS(vis_small, file.path(vispath, "Dlbcl-merge-SCTpostSpotClean_small_fig6e.rds"))

p <- DotPlot(
  vis_small,
  features = rev(DLBCLnChromium_Marker_Gene_List)) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5))

