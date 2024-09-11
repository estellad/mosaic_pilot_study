# Reading mapping results -------------------------------------------------
mapped_all_spot <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Registration", "spot_mapped_cf_region.csv")) # read B1_2 annotated

fig_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/SuppFig_regi_subspot"
###########################################################################
# Spot  ----------------------------------------------
patho_order <- c("Tumor_pure", "Most_likely_tumor", "Tumor_Stroma_mix", 
                 "Intratumoral_Stroma", "Intratumoral_Vessel", "Acellular mucin", 
                 "Tumor_Immune_mix", "Immune_Cell_mix", "Lymphocytes")

# mapped_all_plt <- get_mapped_plt_df(mapped_all_df = mapped_all_spot, resolution = "Spot")
library(SeuratData)
# InstallData("pbmc3k")
pbmc3k <- LoadData("pbmc3k", type = "pbmc3k.final")

df <- mapped_all_spot %>% # 614
  filter(!(Region %in% c("NA", "FALSE"))) %>% # 614
  filter(cell_fraction != "PanCK-") %>% # 463
  dplyr::select(cell_fraction, Region) %>%
  mutate(cell_fraction = ifelse(cell_fraction == "T_cells", "T cells", cell_fraction))

test <- pbmc3k[1:500, 1:463]
test$cell_fraction <- factor(df$cell_fraction, levels = c("Malignant", "Other", "T cells"))
test$Region <- factor(df$Region, levels = patho_order)

patho_anno <- c(
  "Tumor_pure",
  "Most_likely_tumor",
  "Tumor_Immune_mix",
  "Lymphocytes",
  "Immune_Cell_mix",
  "Tumor_Stroma_mix",
  "Intratumoral_Stroma",
  "Intratumoral_Vessel",
  "Acellular mucin"
)

patho_color <- c(
  "#FFA54F",
  "#c47f3d",
  "#00e04e",
  "#4169E1",
  "#b2ff00",
  "#00c3c3",
  "#388E8E",
  "#ba5440",
  "#FFC1C1"
)

names(patho_color) <- patho_anno
patho_color <- patho_color[patho_order]

cell_fraction_order = c("Macro", "Malignant", "Other", "T cells", "PanCK-") # note with PanCK-
cell_fraction_color = c("#9A32CD", "#BC8F8F", "#388E8E", "#4169E1", "#09D0EF")
names(cell_fraction_color) <- cell_fraction_order


p <- SCpubr::do_ChordDiagramPlot(sample = test,
                                 from = "Region",
                                 to = "cell_fraction",
                                 colors.from = c("Tumor_pure" = "#FFA54F", "Most_likely_tumor" = "#c47f3d", "Tumor_Stroma_mix" = "#808000",
                                                 "Intratumoral_Stroma" = "#388E8E", "Intratumoral_Vessel" = "#388E8E", "Acellular mucin" = "#FFC1C1",
                                                 "Tumor_Immune_mix" = "#00e04e", "Immune_Cell_mix" = "#b2ff00", "Lymphocytes" = "#4169E1"),
                                 colors.to = c("Malignant" = "#65DBD7", "Other" = "#388E8E", "T cells" = "#4169E1"),
                                 big.gap = 40,
                                 directional = -1,
                                 alignment = "horizontal",
                                 padding_labels = 0)

p

pdf(file = file.path(fig_path, "Regi_spot_chord_diagram.pdf"),
    width = 13,
    height = 13)
print(p)
dev.off()




