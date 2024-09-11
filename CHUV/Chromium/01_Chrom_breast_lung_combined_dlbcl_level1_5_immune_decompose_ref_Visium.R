library(Seurat)
library(dplyr)
###########################################################################
#      Final Decon Matrix level1_5 (decompose T cell and Macrophage)      #
###########################################################################
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"

# DLBCL -------------------------------------------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

# dlbcl
chrom$level1_5_immune <- case_when(chrom$Level1 == "B" ~ "B cells",   
                                   chrom$Level1 == "Epithelia" ~ "Epithelia",
                                   chrom$Level1 == "Stroma" ~ "Stroma",
                                   chrom$Level1 == "Tumor_DLBCL" ~ "Tumor",
                                   chrom$Harmonised_Level4 %in% c("DC_1", "DC_2", "DC_pc") ~ "Myeloid else",
                                   chrom$Harmonised_Level4 == "NK" ~ "NK",
                                   chrom$Harmonised_Level4 == "Mono_Macro" ~ "Macrophage",
                                   chrom$Harmonised_Level4 %in% c("T_CD4", "T_CD4_reg", "T_CD8", "T_dividing") ~ "T cells") 

# DLBCL order
dlbcl_level1_5_immune_order <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else",  "Macrophage", "T cells", "Tumor") # -> Fig 6
chrom$level1_5_immune <- factor(chrom$level1_5_immune, levels = dlbcl_level1_5_immune_order)

table(chrom$level1_5_immune)
# Epithelia       Stroma      B cells           NK Myeloid else   Macrophage      T cells        Tumor 
#      2945         3388           77          349         1398         3075         3721        24760

saveRDS(chrom, file.path(chrompath, "chrom_dlbcl.rds"))


# Breast -------------------------------------------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))

# breast / lung
chrom$level1_5_immune <- case_when(chrom$Level1 == "B" ~ "B cells",   
                                   chrom$Level1 == "Stroma" ~ "Stroma",
                                   chrom$Level1 == "Tumor_Breast" ~ "Tumor",
                                   chrom$Harmonised_Level4 %in% c("DC_1", "DC_2", "DC_activated", "DC_pc", "Granulocyte", "Mast_cell") ~ "Myeloid else",
                                   chrom$Harmonised_Level4 == "NK" ~ "NK",
                                   chrom$Harmonised_Level4 %in% c("Macrophage", "Monocyte") ~ "Macrophage",
                                   chrom$Harmonised_Level4 %in% c("T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing") ~ "T cells") 

# Breast/lung order
breast_level1_5_immune_order <- c("Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") # Fig 2, 3
chrom$level1_5_immune <- factor(chrom$level1_5_immune, levels = breast_level1_5_immune_order)

table(chrom$level1_5_immune)
# Stroma        Tumor   Macrophage      T cells      B cells           NK Myeloid else 
#   7834         4691         6620        11977         8357          708         2880 
table(chrom$Harmonised_Level4[is.na(chrom$level1_5_immune)])

saveRDS(chrom, file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))


# Lung --------------------------------------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))

# breast / lung
chrom$level1_5_immune <- case_when(chrom$Level1 == "B" ~ "B cells",   
                                   chrom$Level1 == "Stroma" ~ "Stroma",
                                   chrom$Level1 == "Epithelia" ~ "Epithelia",
                                   chrom$Level1 == "Tumor_Lung" ~ "Tumor",
                                   chrom$Harmonised_Level4 %in% c("DC_1", "DC_2", "DC_activated", "DC_pc", "Granulocyte", "Mast_cell") ~ "Myeloid else",
                                   chrom$Harmonised_Level4 == "NK" ~ "NK",
                                   chrom$Harmonised_Level4 %in% c("Macrophage", "Monocyte") ~ "Macrophage",
                                   chrom$Harmonised_Level4 %in% c("T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg", "TNK_dividing") ~ "T cells") 


lung_level1_5_immune_order <- c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else") # Fig 2, 3
chrom$level1_5_immune <- factor(chrom$level1_5_immune, levels = lung_level1_5_immune_order)
table(chrom$level1_5_immune)
# Epithelia       Stroma        Tumor   Macrophage      T cells      B cells           NK Myeloid else 
#       648         7834         2928         6620        11977         8357          708         2880 
table(chrom$Harmonised_Level4[is.na(chrom$level1_5_immune)])

saveRDS(chrom, file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))


