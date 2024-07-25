library(tidyr)

# Note: DLBCL Fibro_Muscle -> Fibro_muscle

# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
save_bs_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/"
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

decon_patho_gather <- NULL
for(i in 1:nsamples){
  sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  patho <- data.frame(Region = sce$Region,
                      Barcode = sce$Barcode)
  
  save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
  results <- read.csv(paste0(save_bs_path_sample, save_names[i], "_spot_level1_5_decon.csv"))
  results <- results %>% rename(Barcode = X)
  
  decon_patho_i <- results %>%
    left_join(patho, by = "Barcode")
  
  decon_patho_gather_i <- gather(decon_patho_i, key = "CellType", value = "Fraction", -c("Barcode", "Region"))
  decon_patho_gather_i$Section <- save_names[i]
  decon_patho_gather <- rbind(decon_patho_gather, decon_patho_gather_i)
}

decon_patho_gather_lung <- decon_patho_gather
decon_patho_gather_breast <- decon_patho_gather
decon_patho_gather_dlbcl <- decon_patho_gather
decon_patho_gather_dlbcl$CellType <- ifelse(decon_patho_gather_dlbcl$CellType == "Fibro_Muscle", 
                                            "Fibro_muscle", decon_patho_gather_dlbcl$CellType)
decon_patho_gather_dlbcl$Section <- case_when(decon_patho_gather_dlbcl$Section == "DLBCL_1" ~ "D1",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_2" ~ "D2",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_3" ~ "D3",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_4" ~ "D4",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_5" ~ "D5",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_6" ~ "D6"
                                              )

vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2c"
write.csv(decon_patho_gather_lung, file.path(vis_decon_path, "vis_lung_decon_long.csv"), row.names = FALSE)
write.csv(decon_patho_gather_breast, file.path(vis_decon_path, "vis_breast_decon_long.csv"), row.names = FALSE)
write.csv(decon_patho_gather_dlbcl, file.path(vis_decon_path, "vis_dlbcl_decon_long.csv"), row.names = FALSE)


