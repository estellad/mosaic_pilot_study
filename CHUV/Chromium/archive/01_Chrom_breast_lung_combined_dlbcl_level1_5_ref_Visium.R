library(Seurat)

##############################################################################
#                               Breast and Lung                              #
##############################################################################
# breast
chrom <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast/breast_final_owkin_annot.rds")
chrom$level1_5 <- ifelse(chrom$Level2 %in% c("Tu_B1", "Tu_B3", "Tu_B4"), "Tumor", chrom$Level2)
chrom$level1_5 <- ifelse(chrom$level1_5 %in% c("Granulocyte", "Mast_cell"), "Granulocyte", chrom$level1_5)
chrom$level1_5 <- droplevels(as.factor(chrom$level1_5))

table(chrom$level1_5) #note: breast has no epithelia
#   B Fibro_muscle  Granulocyte      Myeloid         T_NK        Tumor       Vessel 
# 759         2198          131          758         1287         4691          865 

chrom_breast <- chrom
chrom_breast$sample_id <- as.character(droplevels(as.factor(chrom_breast$sample_id)))
saveRDS(chrom_breast, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_breast.rds")
rm(chrom)

# lung
chrom <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Lung/lung_final_owkin_annot.rds")
chrom$level1_5 <- ifelse(chrom$Level2 %in% c("Tu_L1", "Tu_L2", "Tu_L3", "Tu_L4"), "Tumor", chrom$Level2)
chrom$level1_5 <- ifelse(chrom$level1_5 %in% c("Granulocyte", "Mast_cell"), "Granulocyte", chrom$level1_5)
chrom$level1_5 <- droplevels(as.factor(chrom$level1_5))

table(chrom$level1_5)
#    B    Epithelia Fibro_muscle  Granulocyte      Myeloid         T_NK        Tumor       Vessel 
# 7598          648         3916          623         7988        11398         2928          855

chrom_lung <- chrom
chrom_lung$sample_id <- as.character(droplevels(as.factor(chrom_lung$sample_id)))
saveRDS(chrom_lung, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_lung.rds")
rm(chrom)


## Breast only add lung immune stroma ------------------------------------------
chrom_breast_add_lung_immune_stroma <- merge(chrom_breast, 
                                      chrom_lung[, chrom_lung$level1_5 %in% c("B", "Fibro_muscle", "Granulocyte", "Myeloid", "T_NK", "Vessel")])

saveRDS(chrom_breast_add_lung_immune_stroma, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_breast_add_lung_healthy.rds")



## Lung only add breast immune stroma ------------------------------------------
chrom_lung_add_breast_immune_stroma <- merge(chrom_lung, 
                                             chrom_breast[, chrom_breast$level1_5 %in% c("B", "Fibro_muscle", "Granulocyte", "Myeloid", "T_NK", "Vessel")])

saveRDS(chrom_lung_add_breast_immune_stroma, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_lung_add_breast_healthy.rds")



##############################################################################
#                                    DLBCL                                   #
##############################################################################
disease = "dlbcl"
chrompath <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/DLBCL")
filename <- ifelse(disease == "dlbcl", "DLBCL", disease)
chrom <- readRDS(file.path(chrompath, paste0(filename, "_final_owkin_annot.rds")))

chrom$level1_5 <- ifelse(chrom$Level2 %in% c("Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6"), "Tumor", chrom$Level2)

# table(chrom$level1_5)
# B    Epithelia Fibro_Muscle      Myeloid         T_NK        Tumor       Vessel 
# 77         2945         2954         4473         4070        24760          434 

saveRDS(chrom, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_dlbcl.rds")

# chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
# chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
# chrom$Level3 <- ifelse(chrom$Level3 == "pericytes", "Pericyte", chrom$Level3)
# chrom$Harmonised_Level4 <- ifelse(chrom$Harmonised_Level4 == "pericytes", "Pericyte", chrom$Harmonised_Level4)
# saveRDS(chrom, file.path(chrompath, "chrom_dlbcl.rds"))


###########################################################################
#                   Final Decon Matrix level1_5 (Prev)                    #
###########################################################################

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
table(chrom$level1_5)
# B       Epithelia   Fibro_Muscle                    Myeloid         T_NK        Tumor       Vessel          # DLBCL
# 77      2945        2954                            4473            4070        24760       434


# Definition of healthy cells that breast and lung share: "B", "Fibro_muscle", "Granulocyte", "Myeloid", "T_NK", "Vessel"
chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
chrom$level1_5 <- ifelse(chrom$level1_5 == "Fibro_muscle", "Fibro_Muscle", chrom$level1_5)
chrom$Level2 <- ifelse(chrom$Level2 == "Fibro_muscle", "Fibro_Muscle", chrom$Level2)
table(chrom$level1_5)
# B                   Fibro_Muscle   Granulocyte      Myeloid         T_NK        Tumor       Vessel          # Breast
# 8357                6114           754              8746            12685       4691        1720
saveRDS(chrom, file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))


chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
chrom$level1_5 <- ifelse(chrom$level1_5 == "Fibro_muscle", "Fibro_Muscle", chrom$level1_5)
chrom$Level2 <- ifelse(chrom$Level2 == "Fibro_muscle", "Fibro_Muscle", chrom$Level2)
table(chrom$level1_5)
# B       Epithelia   Fibro_Muscle   Granulocyte      Myeloid         T_NK        Tumor       Vessel          # Lung
# 8357    648         6114           754              8746            12685       2928        1720
saveRDS(chrom, file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))


