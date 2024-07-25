library(ggplot2)
library(dplyr)
library(stringr)
library(Seurat)
# Without pt spec decon, take CellType, Fraction, Section
geo_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/GeoMx/Final_level1_5_decon_results"
geo_decon <- rbind(read.csv(file.path(geo_decon_path, "breast_batched_decon_long.csv")),
                   read.csv(file.path(geo_decon_path, "lung_batched_decon_long.csv")),
                   read.csv(file.path(geo_decon_path, "dlbcl_batched_decon_long.csv")) %>% rename(section_id = patient)
) %>% rename(Section = section_id)

vis_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Manuscript_Figures/Fig2c"
vis_decon <- rbind(read.csv(file.path(vis_decon_path, "vis_breast_decon_long.csv")) %>% 
                     mutate(CellType = ifelse(CellType == "Fibro_Muscle", "Fibro_muscle", CellType)),
                   read.csv(file.path(vis_decon_path, "vis_lung_decon_long.csv")),
                   read.csv(file.path(vis_decon_path, "vis_dlbcl_decon_long.csv"))
)


# Decon average across all spots in Visium versus decon average across all ROIs in GeoMx ----------
geo_test <- geo_decon %>% 
  group_by(cell_fraction, CellType, Section) %>%
  summarise(mean = mean(Fraction)) %>%
  filter(CellType == "Myeloid" & cell_fraction == "Macrophage")


vis_test <- vis_decon %>%
  group_by(Region, CellType, Section) %>%
  summarise(mean = mean(Fraction)) %>%
  filter(CellType == "Myeloid" & Region == "Immune_Cell_mix")


###############################################################################
# Overall Integrated UMAP (level 1.5 - Myeloid)
###############################################################################
# Breast show on UMAP ---------------------------------------------------------
vis_decon_breast <- read.csv(file.path(vis_decon_path, "vis_breast_decon_long.csv")) %>% 
  mutate(CellType = ifelse(CellType == "Fibro_Muscle", "Fibro_muscle", CellType)) 

vis_decon_breast_ <- vis_decon_breast %>%
  mutate(idx = case_when(Section == "B1_2" ~ 1,
                         Section == "B1_4" ~ 2,
                         Section == "B2_2" ~ 3,
                         Section == "B3_2" ~ 4,
                         Section == "B4_2" ~ 5
  )) %>%
  mutate(Barcode = paste0(Barcode, "_", idx))

# After spotclean -----------------------------------------------------
disease = "breast"
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")

seu <- readRDS(paste0(savepath.spotclean, save_rds_name))

seuCD <- data.frame(Barcode = colnames(seu))

seuCD_ <- seuCD %>%
  left_join(vis_decon_breast_ %>% filter(CellType == "Myeloid") %>% select(Barcode, Fraction))

seu$Myeloid_frac <- seuCD_$Fraction


FeaturePlot(seu, "Myeloid_frac") + scale_color_gradientn(colors = c("#F0F0F0","purple", "darkblue"))