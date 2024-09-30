## Note 
# - Breast and lung for Figure 2&3 
# - DLBCL for Figure 6

chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon"
chrom_breast <- readRDS(file.path(chrompath, "chrom_breast.rds"))
table(chrom_breast$Harmonised_Level4)
# B_cell            B_plasma_IGHA1            B_plasma_IGHG1            B_plasma_IGHG3             B_plasma_IGHM             B_plasma_IGKC 
# 79                        86                        53                        55                        98                       365 
# B_plasma_IGLC1          DC_1                      DC_2              DC_activated                     DC_pc       Endothelia_vascular 
# 23                         7                        66                        12                         7                       582 
# Fibroblast             Fibroblast_B3               Granulocyte                Macrophage                 Mast_cell                  Monocyte 
# 732                       780                         3                       237                       128                       429 
# Muscle_smooth              NK                  Pericyte                     T_CD4           T_CD8_exhausted                     T_CTL 
# 686                       101                       283                       252                         9                       548 
# T_CXCL13                     T_reg              TNK_dividing               Tu_B1_MUCL1      Tu_B1_MUCL1_necrosis Tu_B1_MUCL1_transcription 
# 25                       339                        13                       254                       298                       392 
# Tu_B3_CYP4F8                Tu_B4_RHOB 
# 2897                       850 


## Level 1_5
"Epithelia" = "Epithelia"
"Stroma" = "Stroma"
"B cells" = "B"

# Level 4
"Myeloid else" = c("DC_1", "DC_2", "DC_activated", "DC_pc", "Granulocyte", "Mast_cell")
"NK cells" = "NK"
"Macrophage" = c("Macrophage", "Monocyte")
"T cells" = c("T_CD4", "T_CD8_exhausted", "T_CTL ", "T_CXCL13", "T_reg", "TNK_dividing")


# Lung -------------------------------------------------------------------------
chrom_lung <- readRDS(file.path(chrompath, "chrom_lung.rds"))
table(chrom_lung$Harmonised_Level4)
# B_cell                   B_plasma_IGHA1                   B_plasma_IGHG1                   B_plasma_IGHG3                    B_plasma_IGHM 
# 2015                             1048                             1510                              386                              254 
# B_plasma_IGKC                   B_plasma_IGLC1                             DC_1                             DC_2                     DC_activated 
# 2009                              376                              233                             1100                              284 
# DC_pc              Endothelia_vascular                         Epi_lung                       Fibroblast                      Granulocyte 
# 417                              654                              648                             2445                              184 
# Macrophage                        Mast_cell                         Monocyte                    Muscle_smooth                               NK 
# 4611                              439                             1343                             1471                              607 
# Pericyte                            T_CD4                  T_CD8_exhausted                            T_CTL                         T_CXCL13 
# 201                             2053                             1226                             2867                              700 
# T_reg                     TNK_dividing                      Tu_L1_SFTPB                      Tu_L2_FXYD2                       Tu_L3_G0S2 
# 3531                              414                              306                              637                              311 
# Tu_L3_G0S2_immune_signature     Tu_L4_KRT17_immune_signature               Tu_L4_KRT17_mucous             Tu_L4_KRT17_necrosis Tu_L4_KRT17_neutrophil_signature 
# 267                              404                               81                              600                              322 


# DLBCL -------------------------------------------------------------------------
chrom_dlbcl <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))

table(chrom_dlbcl$Harmonised_Level4)
# B_plasma                       DC_1                       DC_2                      DC_pc                 Endothelia               Epi_bronchus 
# 77                        291                        445                        662                        345                        121 
# Epi_mucous Epi_Mucous_surface_gastric       Epi_parietal_gastric               Fibro_Muscle                 Mono_Macro                         NK 
# 1071                       1375                        378                       2954                       3075                        349 
# Pericyte                      T_CD4                  T_CD4_reg                      T_CD8                 T_dividing                 Tu_D1_LMO2 
# 89                       1084                          9                       2384                        244                       3659 
# Tu_D1_RGS13               Tu_D1_SMIM14                 Tu_D2_mito               Tu_D3_BCL2A1             Tu_D3_dividing                Tu_D3_FAM3C 
# 767                       1893                       5624                        325                       1814                       5795 
# Tu_D3_IGHD                Tu_D4_BCL7A                  Tu_D4_PNN                Tu_D5_CCL22                 Tu_D6_BCL2 
# 707                       1470                        535                       1019                       1152 

# table(chrom_dlbcl$level1_5, chrom_dlbcl$Harmonised_Level4)

# -------------------------------------------------------------------------
## Level 1
"Epithelia" = "Epithelia"
"Stroma" = "Stroma"
"B cells" = "B"

# Level 4
"Myeloid else" = c("DC_1", "DC_2", "DC_pc")
"NK cells" = c("NK")
"Macrophage" = "Mono_Macro"
"T cells" = c("T_CD4", "T_CD4_reg", "T_CD8", "T_dividing")


Tumor = level4 # separated in level 4 decon -> Fig 6 separated dotplot
Tumor = level1_5 # renamed level 1.5 decon Tu_DLBCL to Tu_Dn -> Fig6 combined dotplot

# GeoMx ordering: Other, Macrophage, T cells 

# Signature matrix --------------------------------------------------------
if(disease == "breast"){
  # chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
  sigmat_name = "sigmat_breast_add_lung_healthy.rds"
  CF_order <- c("Malignant", "Other", "T cells", "Macrophage")                                            # Fig 2, 3 
}else if(disease == "lung"){
  # chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds")) 
  sigmat_name = "sigmat_lung_add_breast_healthy.rds"
  CF_order <- c("Malignant", "Other", "T cells", "Macrophage")                                            # Fig 2, 3 
}else{
  # chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds")); chrom$patient <- chrom$sample_id
  sigmat_name = "sigmat_dlbcl.rds"
  # CT_order <- c("Epithelia", "Stroma", "B cells", "NK", "Myeloid else",  "Macrophage", "T cells", "Tumor") # Fig 6
  CF_order <- c("B cells", "Other", "T cells", "Macrophage")                                               # Fig 2, 3 
}

CT_order <- c("Epithelia", "Stroma", "Tumor", "Macrophage", "T cells", "B cells", "NK", "Myeloid else")  # Fig 2, 3 