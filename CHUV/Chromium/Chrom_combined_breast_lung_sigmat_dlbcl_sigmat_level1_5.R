library(SpatialDecon)

# breast
chrompath <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/", foldername)
filename <- ifelse(disease == "dlbcl", "DLBCL", disease)
chrom <- readRDS(file.path(chrompath, paste0(filename, "_final_owkin_annot.rds")))
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

# table(chrom$level1_5)
# #    B    Epithelia Fibro_muscle  Granulocyte      Myeloid         T_NK        Tumor       Vessel 
# # 7596          648         3915          919         7910        11398         2873          855 

chrom_lung <- chrom
chrom_lung$sample_id <- as.character(droplevels(as.factor(chrom_lung$sample_id)))
saveRDS(chrom_lung, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_lung.rds")
rm(chrom)

chrom_breast_lung <- merge(chrom_breast, chrom_lung)

table(chrom_breast_lung$level1_5)

matrix <- as(chrom_breast_lung@assays$RNA@counts, "dgCMatrix")

annot_df <- data.frame(
  CellID = colnames(chrom_breast_lung),
  LabeledCellType = chrom_breast_lung$level1_5)

# ## Breast and lung combined ---------------------------
# custom_mtx <- create_profile_matrix(mtx = chrom_breast_lung@assays$RNA@counts,            # cell x gene count matrix
#                                     cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
#                                     cellTypeCol = "LabeledCellType",  # column containing cell type
#                                     cellNameCol = "CellID",           # column containing cell ID/name
#                                     matrixName = "custom_mini_colon", # name of final profile matrix
#                                     outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
#                                     normalize = FALSE,                # Should data be normalized? 
#                                     minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
#                                     minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
#                                     scalingFactor = 5,                # what should all values be multiplied by for final matrix
#                                     discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
# 
# 
# heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
#         labRow = NA, margins = c(10, 5), cexCol = 0.7)

## Breast only ------------------------------------------
annot_df <- data.frame(
  CellID = colnames(chrom_breast),
  LabeledCellType = chrom_breast$level1_5)

custom_mtx <- create_profile_matrix(mtx = chrom_breast@assays$RNA@counts,            # cell x gene count matrix
                                    cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

saveRDS(custom_mtx, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_breast.rds")

heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)


## Lung only ------------------------------------------
annot_df <- data.frame(
  CellID = colnames(chrom_lung),
  LabeledCellType = chrom_lung$level1_5)

custom_mtx <- create_profile_matrix(mtx = chrom_lung@assays$RNA@counts,            # cell x gene count matrix
                                    cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

saveRDS(custom_mtx, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_lung.rds")

heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)






## Breast only add lung immune stroma ------------------------------------------
chrom_breast_add_lung_immune_stroma <- merge(chrom_breast, 
                                      chrom_lung[, chrom_lung$level1_5 %in% c("B", "Fibro_muscle", "Granulocyte", "Myeloid", "T_NK", "Vessel")])

saveRDS(chrom_breast_add_lung_immune_stroma, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_breast_add_lung_healthy.rds")

annot_df <- data.frame(
  CellID = colnames(chrom_breast_add_lung_immune_stroma),
  LabeledCellType = chrom_breast_add_lung_immune_stroma$level1_5)

custom_mtx <- create_profile_matrix(mtx = chrom_breast_add_lung_immune_stroma@assays$RNA@counts,            # cell x gene count matrix
                                    cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

saveRDS(custom_mtx, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_breast_add_lung_healthy.rds")

heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)

## Lung only add breast immune stroma ------------------------------------------
chrom_lung_add_breast_immune_stroma <- merge(chrom_lung, 
                                             chrom_breast[, chrom_breast$level1_5 %in% c("B", "Fibro_muscle", "Granulocyte", "Myeloid", "T_NK", "Vessel")])

saveRDS(chrom_lung_add_breast_immune_stroma, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/chrom_lung_add_breast_healthy.rds")

annot_df <- data.frame(
  CellID = colnames(chrom_lung_add_breast_immune_stroma),
  LabeledCellType = chrom_lung_add_breast_immune_stroma$level1_5)

custom_mtx <- create_profile_matrix(mtx = chrom_lung_add_breast_immune_stroma@assays$RNA@counts,            # cell x gene count matrix
                                    cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

saveRDS(custom_mtx, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_lung_add_breast_healthy.rds")

heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)


# ## Breast only add lung immune ------------------------------------------
# chrom_breast_add_lung_immune <- merge(chrom_breast, 
#                                              chrom_lung[, chrom_lung$level1_5 %in% c("B", "Granulocyte", "Myeloid", "T_NK")])
# annot_df <- data.frame(
#   CellID = colnames(chrom_breast_add_lung_immune),
#   LabeledCellType = chrom_breast_add_lung_immune$level1_5)
# 
# custom_mtx <- create_profile_matrix(mtx = chrom_breast_add_lung_immune@assays$RNA@counts,            # cell x gene count matrix
#                                     cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
#                                     cellTypeCol = "LabeledCellType",  # column containing cell type
#                                     cellNameCol = "CellID",           # column containing cell ID/name
#                                     matrixName = "custom_mini_colon", # name of final profile matrix
#                                     outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
#                                     normalize = FALSE,                # Should data be normalized? 
#                                     minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
#                                     minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
#                                     scalingFactor = 5,                # what should all values be multiplied by for final matrix
#                                     discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
# 
# 
# heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
#         labRow = NA, margins = c(10, 5), cexCol = 0.7)

# ## Lung only add breast immune ------------------------------------------
# chrom_lung_add_breast_immune <- merge(chrom_lung, 
#                                              chrom_breast[, chrom_breast$level1_5 %in% c("B", "Granulocyte", "Myeloid", "T_NK")])
# annot_df <- data.frame(
#   CellID = colnames(chrom_lung_add_breast_immune),
#   LabeledCellType = chrom_lung_add_breast_immune$level1_5)
# 
# custom_mtx <- create_profile_matrix(mtx = chrom_lung_add_breast_immune@assays$RNA@counts, # cell x gene count matrix
#                                     cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
#                                     cellTypeCol = "LabeledCellType",  # column containing cell type
#                                     cellNameCol = "CellID",           # column containing cell ID/name
#                                     matrixName = "custom_mini_colon", # name of final profile matrix
#                                     outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
#                                     normalize = FALSE,                # Should data be normalized? 
#                                     minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
#                                     minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
#                                     scalingFactor = 5,                # what should all values be multiplied by for final matrix
#                                     discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.
# 
# 
# heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
#         labRow = NA, margins = c(10, 5), cexCol = 0.7)
rownames(custom_mtx)[grepl("^MT-", rownames(custom_mtx))]
table(grepl("^MT-", rownames(custom_mtx)))


rownames(custom_mtx)[grepl("^MT", rownames(custom_mtx))]
table(grepl("^MT", rownames(custom_mtx)))



# DLBCL -------------------------------------------------------------------
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

annot_df <- data.frame(
  CellID = colnames(chrom),
  LabeledCellType = chrom$level1_5)

custom_mtx <- create_profile_matrix(mtx = chrom@assays$RNA@counts,            # cell x gene count matrix
                                    cellAnnots = annot_df,  # cell annotations with cell type and cell name as columns 
                                    cellTypeCol = "LabeledCellType",  # column containing cell type
                                    cellNameCol = "CellID",           # column containing cell ID/name
                                    matrixName = "custom_mini_colon", # name of final profile matrix
                                    outDir = NULL,                    # path to desired output directory, set to NULL if matrix should not be written
                                    normalize = FALSE,                # Should data be normalized? 
                                    minCellNum = 5,                   # minimum number of cells of one type needed to create profile, exclusive
                                    minGenes = 10,                    # minimum number of genes expressed in a cell, exclusive
                                    scalingFactor = 5,                # what should all values be multiplied by for final matrix
                                    discardCellTypes = FALSE)# should cell types be filtered for types like mitotic, doublet, low quality, unknown, etc.

saveRDS(custom_mtx, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/sigmat_dlbcl.rds")

heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"),
        labRow = NA, margins = c(10, 5), cexCol = 0.7)

#################################################
#         Final Decon Matrix level1_5           #
#################################################

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


