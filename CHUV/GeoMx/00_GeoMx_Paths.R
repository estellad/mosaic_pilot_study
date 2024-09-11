# general paths
absolute_path_urb <- "/data/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/"
absolute_path_cur <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/"

data_path <- paste0(absolute_path_urb, "mosaic_pilot/")

if(disease == "breast"){
  # Breast ----------------------------------------------------
  sample_id <- c("1FHZ", "1GVR", "1256", "OPHI", "OPHI")
  section_id <- c("B4_1", "B3_1", "B2_1", "B1_1", "B1_3")
  geomx_name <- c("B4_1_1FHZ", "B3_1_1GVR", "B2_1_1256", "B1_1_OPHI", "B1_3_OPHI")
  nsamples = 5
  decon_file_pattern = "_B"
  # dims <- list(24, 23, 24, 24, 24) # +3 negative control templates
}else if(disease == "lung"){
  # Lung ----------------------------------------------------
  sample_id <- c("1GA2", "1G73", "1G73", "0WMU", "0PSV")
  section_id <- c("L4_3", "L3_1", "L3_3", "L2_1", "L1_1")
  geomx_name <- c("L4_3_1GA2", "L3_1_1G73", "L3_3_1G73", "L2_1_0WMU", "L1_1_0PSV")
  nsamples = 5
  decon_file_pattern = "_L"
  # dims <- list(26, 23, 24, 24, 26) # +3 negative control templates
}else if(disease == "dlbcl"){
  # DLBCL ---------------------------------------------------
  dcc_path <- paste0(absolute_path_urb, "mosaic_pilot_2/geomx/dcc")
  section_id <- c("D1", "D2", "D3", "D4", "D5", "D6")
  nsamples = 6
  decon_file_pattern = "_D"
  # 143 samples
}

save_path_raw <- paste0(absolute_path_urb, "env/geomx_rds/", disease, "/raw/")
save_path_qcd <- paste0(absolute_path_urb, "env/geomx_rds/", disease, "/qcd/")
save_path_all <- paste0(absolute_path_urb, "env/geomx_rds/")

