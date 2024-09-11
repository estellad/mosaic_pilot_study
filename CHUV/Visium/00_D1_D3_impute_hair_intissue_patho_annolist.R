# impute D1 D3 hair region with patho annotation list ----------------------
D1_anno_list <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_1.csv")
D3_anno_list <- read.csv("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations/DLBCL_3.csv")
n = 1
n = 3

read_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot_2/visium/", 
                    sample_id[n], "/",
                    sample_id[n], "/",
                    "outs/spatial/tissue_positions.csv")
to_impute <- read.csv(read_path)

to_impute %>%
  mutate(in_tissue = ifelse(barcode %in% D1_anno_list$Barcode, 1, 0))


save_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/Visium_D1_D3_hair/",
                    paste0("D", n), "/",
                    "outs/spatial/tissue_positions.csv")
write.csv(to_impute, save_path)