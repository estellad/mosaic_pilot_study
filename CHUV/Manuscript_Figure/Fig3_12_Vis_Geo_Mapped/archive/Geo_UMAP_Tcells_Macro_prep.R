library(tidyr)
library(dplyr)
# # Sanity (PanCK- samples included) -----------------------------------------
# breast_geo_level4_decon <- read.csv(file.path(deconresultpath, "breast_batched_decon_long.csv"))
# lung_geo_level4_decon <- read.csv(file.path(deconresultpath, "lung_batched_decon_long.csv"))
# dlbcl_geo_level4_decon <- read.csv(file.path(deconresultpath, "dlbcl_batched_decon_long.csv"))
# length(unique(breast_geo_level4_decon$sample))
# length(unique(lung_geo_level4_decon$sample))
# length(unique(dlbcl_geo_level4_decon$sample))

# -------------------------------------------------------------------------
# disease = "breast"
# disease = "lung"
disease = "dlbcl"

source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/GeoMx_init.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/mosaic_pilot_study/CHUV/GeoMx/00_GeoMx_Paths.R")

datapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/GeoMx/GeoMx_Normed_Batched/"
spe_ruv <- readRDS(file.path(datapath, paste0(disease, "_spe_ruv.rds")))

deconresultpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level4_decon_results_pt_specific"
geo_level4_decon <- read.csv(file.path(deconresultpath, 
                                       paste0(disease, "_batched_decon_long.csv")))

if(disease %in% c("breast", "lung")){
  test <- geo_level4_decon %>%
    filter(CellType %in% c("Macrophage", "T_CD4", "T_CD8_exhausted", "T_CTL", "T_CXCL13", "T_reg")) %>%
    select(sample, Fraction, CellType) %>%
    tidyr::spread(CellType, Fraction) %>%
    dplyr::rename(sample_id2 = sample) %>%
    mutate(`T cells` = T_CD4 + T_CD8_exhausted + T_CTL + T_CXCL13 + T_reg) %>%
    select(sample_id2, Macrophage, `T cells`)
}else{ # dlbcl
  test <- geo_level4_decon %>%
    filter(CellType %in% c("Mono_Macro", "T_CD4", "T_CD4_reg", "T_CD8", "T_dividing")) %>%
    select(sample, Fraction, CellType) %>%
    tidyr::spread(CellType, Fraction) %>%
    dplyr::rename(sample_id2 = sample) %>%
    mutate(`T cells` = T_CD4 + T_CD4_reg + T_CD8 + T_dividing) %>%
    dplyr::rename(Macrophage = "Mono_Macro") %>%
    select(sample_id2, Macrophage, `T cells`)
}


CD <- as.data.frame(colData(spe_ruv)) %>%
  select(sample_id2)

CD <- CD %>%
  left_join(test)

spe_ruv$Macrophage_frac <- CD$Macrophage
spe_ruv$Tcells_frac <- CD$`T cells`

saveRDS(spe_ruv, file.path(datapath, paste0(disease, "_spe_ruv.rds")))






