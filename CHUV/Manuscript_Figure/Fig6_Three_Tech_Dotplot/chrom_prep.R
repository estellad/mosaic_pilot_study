library(Seurat)
library(dplyr)
library(tidyr)
library(tidyverse)

CT_order_nontumor <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else",  "Macrophage", "T cells") # , "Tumor") # Fig 6
# Chromium ----------------------------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

## Separate tumor subtypes
# Fig 6d
chrom$level1_5_immune_tumor_subtypes <- ifelse(chrom$level1_5_immune == "Tumor", chrom$Harmonised_Level4, as.character(chrom$level1_5_immune))
table(chrom$level1_5_immune_tumor_subtypes) 
# B cells      Epithelia     Macrophage   Myeloid else             NK         Stroma        T cells     Tu_D1_LMO2    Tu_D1_RGS13   Tu_D1_SMIM14 
# 77           2945           3075           1398            349           3388           3721           3659            767           1893 
# Tu_D2_mito   Tu_D3_BCL2A1 Tu_D3_dividing    Tu_D3_FAM3C     Tu_D3_IGHD    Tu_D4_BCL7A      Tu_D4_PNN    Tu_D5_CCL22     Tu_D6_BCL2 
# 5624            325           1814           5795            707           1470            535           1019           1152 

chrom$level1_5_immune_tumor_subtypes <- factor(chrom$level1_5_immune_tumor_subtypes, 
                                      levels = c(CT_order_nontumor, "Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14",
                                                 "Tu_D2_mito", "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", 
                                                 "Tu_D3_IGHD", "Tu_D4_BCL7A", "Tu_D4_PNN", "Tu_D5_CCL22", "Tu_D6_BCL2"))

## Combined tumor subtypes
## Fig 6e
chrom$level1_5_immune_tumor <- ifelse(chrom$level1_5_immune == "Tumor", chrom$Level2, as.character(chrom$level1_5_immune))
table(chrom$level1_5_immune_tumor) 
# B cells    Epithelia   Macrophage Myeloid else           NK       Stroma      T cells        Tu_D1        Tu_D2        Tu_D3        Tu_D4 
# 77         2945         3075         1398          349         3388         3721         6319         5624         8641         2005 
# Tu_D5        Tu_D6 
# 1019         1152 

chrom$level1_5_immune_tumor <- factor(chrom$level1_5_immune_tumor, 
                                      levels = c(CT_order_nontumor, "Tu_D1", "Tu_D2", "Tu_D3", "Tu_D4", "Tu_D5", "Tu_D6"))


saveRDS(chrom, file.path(chrompath, "chrom_dlbcl.rds"))


