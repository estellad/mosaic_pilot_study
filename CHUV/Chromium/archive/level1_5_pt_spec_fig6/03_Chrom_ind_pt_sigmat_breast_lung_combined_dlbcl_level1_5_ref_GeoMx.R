library(SpatialDecon)
library(dplyr)
library(Seurat)
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
sigmatpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/signature_matrix/level1_5/"
healthy_cells <- c("B", "Epithelia", "Fibro_Muscle", "Granulocyte", "Mast_cell", "Myeloid", "T_NK", "Vessel")

# helper --------------------------------------------------------------------
getSigmat <- function(chrom, pt_name = "L1"){
  if(pt_name == "B2"){ # Only B2, take tumor as all other pts tumor combined.
    chrom$Level2 <- ifelse(chrom$Level2 %in% c("Tu_B1", "Tu_B3", "Tu_B4"), "Tu_B2", chrom$Level2)
  }
  
  chrom_pt <- chrom[, chrom$Level2 %in% c(healthy_cells, paste0("Tu_", pt_name))]
  print(table(chrom_pt$level1_5))
  
  annot_df <- data.frame(
    CellID = colnames(chrom_pt),
    LabeledCellType = chrom_pt$level1_5
  )
  
  custom_mtx <- create_profile_matrix(mtx = chrom_pt@assays$RNA@counts,            # cell x gene count matrix
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
  
  saveRDS(custom_mtx, file.path(sigmatpath, paste0("sigmat_", pt_name, ".rds")))
  p <- heatmap(sweep(custom_mtx, 1, apply(custom_mtx, 1, max), "/"), labRow = NA, margins = c(10, 5), cexCol = 0.7)
  print(p)
  return(custom_mtx)
}


# ## Breast only add lung immune stroma -------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds")) # table(chrom$patient, chrom$Level2)
#       B Fibro_Muscle Granulocyte Mast_cell Myeloid T_NK Tu_B1 Tu_B3 Tu_B4 Vessel
# B1  714          818           2        98     247  974   944     0     0    140
# B2    5           15           0         0      24    0     0     0     0      6
# B3   27          791           1        30     284  215     0  2897     0    484
# B4   13          574           0         0     203   98     0     0   850    235
# L1   38           51           0        11      90  202     0     0     0     56
# L2  411           80           2         5     107  403     0     0     0     39
# L3 1567         2689         135       107    5569 6885     0     0     0    228
# L4 5582         1096          47       316    2222 3908     0     0     0    532
sigmatB1 <- getSigmat(chrom, pt_name = "B1")
sigmatB2 <- getSigmat(chrom, pt_name = "B2")
sigmatB3 <- getSigmat(chrom, pt_name = "B3")
sigmatB4 <- getSigmat(chrom, pt_name = "B4")


## Lung only add breast immune stroma -------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds")) # table(chrom$patient, chrom$Level2)
#       B Epithelia Fibro_Muscle Granulocyte Mast_cell Myeloid T_NK Tu_L1 Tu_L2 Tu_L3 Tu_L4 Vessel
# B1  714         0          818           2        98     247  974     0     0     0     0    140
# B2    5         0           15           0         0      24    0     0     0     0     0      6
# B3   27         0          791           1        30     284  215     0     0     0     0    484
# B4   13         0          574           0         0     203   98     0     0     0     0    235
# L1   38        48           51           0        11      90  202   306     0     0     0     56
# L2  411        72           80           2         5     107  403     0   637     0     0     39
# L3 1567        46         2689         135       107    5569 6885     0     0   578     0    228
# L4 5582       482         1096          47       316    2222 3908     0     0     0  1407    532
sigmatL1 <- getSigmat(chrom, pt_name = "L1")
sigmatL2 <- getSigmat(chrom, pt_name = "L2")
sigmatL3 <- getSigmat(chrom, pt_name = "L3")
sigmatL4 <- getSigmat(chrom, pt_name = "L4")


# DLBCL  ------------------------------------------------------------------
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")); chrom$patient <- chrom$sample_id # table(chrom$patient, chrom$Level2)
#            B Epithelia Fibro_Muscle Myeloid T_NK Tu_D1 Tu_D2 Tu_D3 Tu_D4 Tu_D5 Tu_D6 Vessel
# DLBCL_1   13         5         1112     862  535  6319     0     0     0     0     0     90
# DLBCL_2    0         5          191     748  494     0  5624     0     0     0     0     35
# DLBCL_3   27        12         1050    1876 2638     0     0  8641     0     0     0    131
# DLBCL_4   32      2541          295     304   79     0     0     0  2005     0     0     79
# DLBCL_5    2       127          287     638  314     0     0     0     0  1019     0     98
# DLBCL_6    3       255           19      45   10     0     0     0     0     0  1152      1
sigmatD1 <- getSigmat(chrom, pt_name = "D1")
sigmatD2 <- getSigmat(chrom, pt_name = "D2")
sigmatD3 <- getSigmat(chrom, pt_name = "D3")
sigmatD4 <- getSigmat(chrom, pt_name = "D4")
sigmatD5 <- getSigmat(chrom, pt_name = "D5")
sigmatD6 <- getSigmat(chrom, pt_name = "D6")


# # Sanity checks -----------------------------------------------------------
# list.files(sigmatpath)
# sigmat_dlbcl_combined <- readRDS(file.path(sigmatpath, "sigmat_dlbcl.rds"))
# heatmap(sweep(sigmat_dlbcl_combined, 1, apply(sigmat_dlbcl_combined, 1, max), "/"), labRow = NA, margins = c(10, 5), cexCol = 0.7) # Tumor not so clear as breast/lung, whether combine or pt-specific



