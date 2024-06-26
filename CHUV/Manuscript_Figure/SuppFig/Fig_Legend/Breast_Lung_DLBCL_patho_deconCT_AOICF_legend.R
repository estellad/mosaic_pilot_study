disease = "dlbcl"
disease = "breast"
disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
save_path_legend <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Figure_Legend"


## Get level 4 decon majority vote legend ---------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

decon_ct <- NULL
for (i in 1:nsamples){
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  decon_ct_i <- levels(sce[["Level4_decon_max"]])
  
  decon_ct <- c(decon_ct, decon_ct_i)
}

decon_ct <- unique(decon_ct)

level4_decon_cols_max <- level4_cellcolors[decon_ct]

# Breast + Lung  ------------------------------------------------------------------
## Breast
decon_ct_breast <- decon_ct
level4_decon_cols_max_breast <- level4_decon_cols_max
# > decon_ct
# [1] "B_plasma_IGHM"        "Endothelia_vascular"  "Fibroblast"           "Macrophage"           "Pericyte"             "T_CD4"               
# [7] "T_CTL"                "Tu_B1_MUCL1_necrosis" "Mix"                  "Granulocyte"          "Muscle_smooth"        "T_reg"               
# [13] "TNK_dividing"         "Tu_B2"                "Fibroblast_B3"        "Tu_B3_CYP4F8"         "DC_1"                 "DC_2"                
# [19] "Monocyte"             "Tu_B4_RHOB" 
# > level4_decon_cols_max
# B_plasma_IGHM  Endothelia_vascular           Fibroblast           Macrophage             Pericyte                T_CD4                T_CTL 
# "#FF83FA"            "#FF7256"            "#388E8E"            "#FF9900"            "#A0522D"            "#4169E1"            "#009ACD" 
# Tu_B1_MUCL1_necrosis                  Mix          Granulocyte        Muscle_smooth                T_reg         TNK_dividing                Tu_B2 
# "#00EE76"            "#5A5A5A"            "#FFC125"            "#668B8B"            "#87CEFF"            "#97FFFF"            "#66FF66" 
# Fibroblast_B3         Tu_B3_CYP4F8                 DC_1                 DC_2             Monocyte           Tu_B4_RHOB 
# "#008080"            "#EEEE00"            "#C71585"            "#FFAEB9"            "#FF8000"            "#FFD700" 

## Lung
decon_ct_lung <- decon_ct
level4_decon_cols_max_lung <- level4_decon_cols_max

# Combine
decon_ct <- unique(sort(c(decon_ct_breast, decon_ct_lung)))
level4_decon_cols_max <- c(level4_decon_cols_max_breast, level4_decon_cols_max_lung)[decon_ct]

plt_df <- data.frame(
  CT = decon_ct,
  x = 1:length(decon_ct),
  y = 1:length(decon_ct)
)

p_ct <- ggplot(plt_df, aes(x = x, y = y, color = CT)) + 
  geom_point(size = 5) + 
  theme_bw() + 
  scale_color_manual(values = level4_decon_cols_max) + 
  theme(legend.position = "right") + 
  guides(color = guide_legend(ncol = 1, byrow = TRUE))

pdf(file.path(save_path_legend, "Breast_Lung_decon_majority_vote_CT.pdf"),
    height = 20, width = 5)
print(p_ct)
dev.off()

# DLBCL  ------------------------------------------------------------------
plt_df <- data.frame(
  CT = decon_ct,
  x = 1:length(decon_ct),
  y = 1:length(decon_ct)
)

p_ct <- ggplot(plt_df, aes(x = x, y = y, color = CT)) + 
  geom_point(size = 5) + 
  theme_bw() + 
  scale_color_manual(values = level4_decon_cols_max) + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

pdf(file.path(save_path_legend, "DLBCL_decon_majority_vote_CT.pdf"),
    height = 5, width = 20)
print(p_ct)
dev.off()

## Get pathology legend ---------------------------------
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

patho <- NULL
for (i in 1:6){
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  patho_i <- unique(sce[["Region"]])
  
  patho <- c(patho, patho_i)
}

patho <- unique(patho)

patho_color <- dlbcl_patho_color[patho]

plt_df <- data.frame(
  CT = patho,
  x = 1:length(patho),
  y = 1:length(patho)
)

p_patho <- ggplot(plt_df, aes(x = x, y = y, color = CT)) + 
  geom_point(size = 5) + 
  theme_bw() + 
  scale_color_manual(values = patho_color) + 
  theme(legend.position = "bottom") + 
  guides(color = guide_legend(nrow = 1, byrow = TRUE))


pdf(file.path(save_path_legend, "DLBCL_patho.pdf"),
    height = 5, width = 10)
print(p_patho)
dev.off()



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

pdf(file = file.path(save_path_legend, "Breast_Lung_patho.pdf"),
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


pdf(file = file.path(save_path_legend, "Breast_Lung_cell_fraction.pdf"),
    width = 4,
    height = 5)
print(p)
dev.off()
