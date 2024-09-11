library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)

# After spotclean -----------------------------------------------------
# disease = "breast"
# disease = "lung"
disease = "dlbcl"
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")

seu <- readRDS(paste0(savepath.spotclean, save_rds_name))

###############################################################################
# Individual decon (level 4 - Macrophage)
###############################################################################
if(disease == "breast"){
  nsamples = 5
  vis_sample_list <- c("B1_2", "B1_4", "B2_2", "B3_2", "B4_2")
  foldername = str_to_title(disease)
}else if(disease == "lung"){
  nsamples = 5
  vis_sample_list <- c("L1_2", "L1_4", "L2_2", "L3_2", "L4_2")
  foldername = str_to_title(disease)
}else{ # dlbcl
  nsamples = 6
  vis_sample_list <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  foldername = toupper(disease)
}


level4_decon <- NULL
for(i in 1:nsamples){
  if(disease %in% c("breast", "lung")){
    level4_decon_i <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", 
                                      foldername, "/", vis_sample_list[i], "/", vis_sample_list[i], "_spot_Level4_decon.csv")) %>%
      select(X, Granulocyte, Macrophage, Mast_cell, Monocyte, DC_1, DC_2, DC_activated, DC_pc,   # Macrophage
             T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg,                                     # T cells
             TNK_dividing, NK) %>%                                                               # NK cells
      dplyr::rename(Barcode = X) %>%
      mutate(Section = vis_sample_list[i],
             `T cells` = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg)
  }else{ # dlbcl
    level4_decon_i <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/", 
                                      foldername, "/", vis_sample_list[i], "/", vis_sample_list[i], "_spot_Level4_decon.csv")) %>%
      select(X, Mono_Macro,                                                                      # Macrophage
             T_CD4, T_CD4_reg, T_CD8, T_dividing,                                                # T cells
             NK) %>%                                                                             # NK cells
      dplyr::rename(Barcode = X) %>%
      mutate(Section = vis_sample_list[i],
             Macrophage = Mono_Macro, 
             `T cells` = T_CD4 + T_CD4_reg + T_CD8 + T_dividing)
  }
  level4_decon <- rbind(level4_decon, level4_decon_i)
}

level4_decon_ <- level4_decon %>%
  mutate(idx = case_when(Section == vis_sample_list[1] ~ 1,
                         Section == vis_sample_list[2] ~ 2,
                         Section == vis_sample_list[3] ~ 3,
                         Section == vis_sample_list[4] ~ 4,
                         Section == vis_sample_list[5] ~ 5
                         ,Section == vis_sample_list[6] ~ 6
  )) %>%
  mutate(Barcode = paste0(Barcode, "_", idx))


# Merge in to seu -------------------------------
seuCD <- data.frame(Barcode = colnames(seu))

seuCD_ <- seuCD %>%
  left_join(level4_decon_ %>% select(Barcode, Macrophage, `T cells`))

seu$Macrophage_frac <- seuCD_$Macrophage
seu$Tcells_frac <- seuCD_$`T cells`

plot(density(level4_decon$Macrophage))
plot(density(level4_decon$`T cells`))
# summary(level4_decon$Macrophage)

saveRDS(seu, paste0(savepath.spotclean, save_rds_name))






