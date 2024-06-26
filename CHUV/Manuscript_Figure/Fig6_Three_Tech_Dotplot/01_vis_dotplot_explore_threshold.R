library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)
vis <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_Integrated_Owkin/DLBCL_Post_SpotClean/DLBCL-merge-SCTpostSpotClean.rds")
# table(chrom$Harmonised_Level4, chrom$Level1)

# # From level_1_5 healthy cells decon ------------------
# decon_combine <- case_when(
#   c("Fibro_Muscle", "Vessel") ~ "Stroma")
# 
# # From level_4 tumor cells decon ------------------
# decon_combine <- case_when(
#   c("Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14") ~ "Tu_D1",
#   "Tu_D2_mito" ~ "Tu_D2",
#   c("Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD") ~ "Tu_D3",
#   c("Tu_D4_BCL7A", "Tu_D4_PNN") ~ "Tu_D4", 
#   "Tu_D5_CCL22" ~ "Tu_D5",
#   "Tu_D6_BCL2" ~ "Tu_D6")

"Tu_D1" = c("Tu_D1_LMO2", "Tu_D1_RGS13", "Tu_D1_SMIM14");
"Tu_D2" = "Tu_D2_mito";
"Tu_D3" = c("Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD");
"Tu_D4" = c("Tu_D4_BCL7A", "Tu_D4_PNN"); 
"Tu_D5" = "Tu_D5_CCL22";
"Tu_D6" = "Tu_D6_BCL2"

# Decon levels ------------------------------------------------------------
# Healthy cells #################################
healthy_long_all <- NULL
for(i in 1:6){
  healthy <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
                      "/DLBCL_", i, "_spot_level1_5_decon.csv"))
  
  healthy_long <- healthy %>%
    mutate(Stroma = Fibro_Muscle + Vessel) %>%
    select(X, Tumor, Stroma, Myeloid, B, T_NK, Epithelia) %>%
    gather(key = "Cell_Type", value = "Fraction", -X) %>%
    mutate(X = paste0("DLBCL_", i, "_", X)) %>%
    filter(Cell_Type != "Tumor")
  
  healthy_long_all <- rbind(healthy_long_all, healthy_long)
}

# Tumor cells (separate) ###################################
tumor_long_all <- NULL
for(i in 1:6){
  tumor <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i, 
                             "/DLBCL_", i, "_spot_Level4_decon.csv"))
  
  tumor_long <- tumor %>%
    select(X, all_of(get(paste0("Tu_D", i)))) %>%
    gather(key = "Cell_Type", value = "Fraction", -X) %>%
    mutate(X = paste0("DLBCL_", i, "_", X)) 
  
  tumor_long_all <- rbind(tumor_long_all, tumor_long)
}


# Tumor cells (combined) ###################################
tumor_long_all_ <- NULL
for(i in 1:6){
  tumor <- read.csv(paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/DLBCL/DLBCL_", i,
                           "/DLBCL_", i, "_spot_Level4_decon.csv"))

  # Merge to one tumor label per patient
  tumor_long_ <- tumor %>%
    select(X, all_of(get(paste0("Tu_D", i)))) %>%
    column_to_rownames("X") %>%
    mutate(Tu_D = rowSums(.)) %>%
    select(Tu_D) %>%
    rownames_to_column("X") %>%
    gather(key = "Cell_Type", value = "Fraction", -X) %>%
    mutate(X = paste0("DLBCL_", i, "_", X)) %>%
    mutate(Cell_Type = paste0(Cell_Type, i))

  tumor_long_all_ <- rbind(tumor_long_all_, tumor_long_)
}

# ggplot(healthy_long_all %>% filter(Fraction > 0.7), aes(x = Fraction, y = Cell_Type)) +
# ggplot(tumor_long_all %>% filter(Fraction > 0.7), aes(x = Fraction, y = Cell_Type)) +
ggplot(tumor_long_all_ %>% filter(Fraction > 0.7), aes(x = Fraction, y = Cell_Type)) +
  geom_density_ridges(scale = 4) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_ridges()

healthy_long_all %>% filter(Fraction > 0.7) %>% group_by(Cell_Type) %>% summarise(n = n())
tumor_long_all %>% filter(Fraction > 0.7) %>% group_by(Cell_Type) %>% summarise(n = n())
tumor_long_all_ %>% filter(Fraction > 0.5) %>% group_by(Cell_Type) %>% summarise(n = n())
