
# Pathology ---------------------------------------------------------------
patho_anno <- c(
  "Breast_Carcinoma_in_situ",
  "Breast_Tumor_pure",
  "Breast_Most_likely_tumor",
  
  "Lung_Tumor_pure",
  "Lung_Most_likely_tumor",
  
  "Tumor_TIL",
  "Tumor_Immune_mix",
  "Lymphocytes",
  "Immune_Cell_mix",
  "Tumor_Stroma_mix",
  "Large_Vessel",
  "Intratumoral_Stroma",
  "Intratumoral_Vessel",
  'Normal_Lung_Epithelium',
  "Acellular mucin",
  "Necrosis_Debris"# ,
  # "Artefact_Fold_exclude"
)

patho_color <- c(
  "#00E5A1",
  "#03C78C",
  "#62e1bb",
  
  "#FFA54F",
  "#c47f3d",
  
  "#2da3f2",
  "#00e04e",
  "#4169E1",
  "#b2ff00",
  "#00c3c3",
  "#FF7256",
  "#388E8E",
  '#ba5440',
  "#BC8F8F",
  "#FFC1C1",
  "#666666"# ,
  # "#cccccc"
)

names(patho_color) <- patho_anno

df <- data.frame(
  x = sample(1:20, length(patho_color)),
  y = sample(1:20, length(patho_color)),
  annote = factor(patho_anno, levels = patho_anno)
)

p <- ggplot(df, aes(x = x, y = y, colour = annote),
       size = 5) + 
  geom_point() + 
  theme_bw()+
  scale_color_manual(name = "Pathology", values = patho_color[names(table(df$annote))]) + 
  guides(color = guide_legend(override.aes = list(size = 5)))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures"
pdf(file = file.path(paste0(fig_path, "/SuppFig_regi_spatial"), "SuppFig_spatial_regi_Patho_legend.pdf"),
    width = 8,
    height = 10)
print(p)
dev.off()



# Cell fraction -----------------------------------------------------------
cell_fraction_order = c("Macro", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
names(cell_fraction_color) <- cell_fraction_order

df <- data.frame(
  x = sample(1:20, length(cell_fraction_color)),
  y = sample(1:20, length(cell_fraction_color)),
  annote = factor(cell_fraction_order, levels = cell_fraction_order)
)

p <- ggplot(df, aes(x = x, y = y, colour = annote),
            size = 5) + 
  geom_point() + 
  theme_bw()+
  scale_color_manual(name = "Cell fraction", values = cell_fraction_color[names(table(df$annote))]) + 
  guides(color = guide_legend(override.aes = list(size = 5)))

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures"
pdf(file = file.path(paste0(fig_path, "/SuppFig_regi_spatial"), "SuppFig_spatial_regi_Cellfrac_legend.pdf"),
    width = 4,
    height = 5)
print(p)
dev.off()










