library(zellkonverter)
disease = "breast"
disease = "lung"
disease = "dlbcl"

rdspath <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
adatapath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd_anndata"

if(disease == "breast"){
  sample_names <- c("B1_2",      "B1_4",      "B2_2",      "B3_2",      "B4_2")  
  nsamples = 5
}else if(disease == "lung"){
  sample_names <- c("L1_2",      "L1_4",      "L2_2",      "L3_2",      "L4_2")   
  nsamples = 5
}else{
  sample_names <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
  nsamples = 6
}
 

for(i in 1:nsamples){
  sce <- readRDS(file.path(rdspath, paste0(sample_names[i], "_qcd.rds")))
  writeH5AD(sce, file = file.path(adatapath, paste0(sample_names[i], "_qcd.h5ad")))
}

