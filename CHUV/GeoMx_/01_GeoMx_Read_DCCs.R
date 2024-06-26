
source(paste0(absolute_path_urb, "env/ydong/GeoMx/GeoMx_init.R"))
######################################################################
#                              Merge all                             #
######################################################################

# DCCs --------------------------------------------------------------------
if(disease == "breast"){
  source(paste0(absolute_path_urb, "env/ydong/GeoMx/00_GeoMx_Paths.R"))
  
  DCCs_breast <- c()
  for(n in 1:5){
    # Read Data ---------------------------------------------------------------
    DCCFiles <- paste0(data_path, disease, "/", sample_id[n], "/geomx/",  section_id[n], "/dcc/")
    DCCFiles <- dir(DCCFiles, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
    DCCs_breast <- c(DCCs_breast, DCCFiles)
  }
  
  length(DCCs_breast) # 119
  DCCs <- DCCs_breast
}else if(disease == "lung"){
  source(paste0(absolute_path_urb, "env/ydong/GeoMx/00_GeoMx_Paths.R"))
  
  DCCs_lung <- c()
  for(n in 1:5){
    # Read Data ---------------------------------------------------------------
    DCCFiles <- paste0(data_path, disease, "/", sample_id[n], "/geomx/",  section_id[n], "/dcc/")
    DCCFiles <- dir(DCCFiles, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
    DCCs_lung <- c(DCCs_lung, DCCFiles)
  }
  
  length(DCCs_lung) # 123
  DCCs <- DCCs_lung
}else if(disease == "dlbcl"){
  source(paste0(absolute_path_urb, "env/ydong/GeoMx/00_GeoMx_Paths.R"))
  
  DCC_dlbcl <- dir(dcc_path, pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
  DCCs <- DCC_dlbcl
}


# LW ----------------------------------------------------------------------
if(disease %in% c("breast", "lung")){
  owkinLW <- readxl::read_xlsx(paste0(absolute_path_cur, "Owkin_Pilot_Data/GeoMx/owkinLW.xlsx"), sheet = "merged_LabWorkSheet")
  
  library(dplyr)
  owkinLW_new <- owkinLW %>% # 242 ROIs
    filter(!is.na(ROI_category)) 
  
  anno_path <- paste0(absolute_path_urb, "env/geomx_rds/", "owkinLW", ".xlsx")
  write_xlsx(owkinLW_new, path = anno_path)
  SampleAnnotationFile <- anno_path
  
  experimentDataColNames = "Panel"
}else if(disease == "dlbcl"){
  anno_path <-  paste0(absolute_path_urb, "env/geomx_rds/dlbcl/raw/", "DLBCL_LabWorkSheet_owkin_ED", ".xlsx")
  SampleAnnotationFile <- anno_path
  
  experimentDataColNames = "panel"
}

# PKC ---------------------------------------------------------------------
PKCFiles <- dir(data_path, pattern = ".pkc$", full.names = TRUE, recursive = TRUE)

# Read Data ---------------------------------------------------------------
Data = readNanoStringGeoMxSet(dccFiles = DCCs,
                              pkcFiles = PKCFiles,
                              phenoDataFile = SampleAnnotationFile,
                              phenoDataSheet = "Sheet1",
                              phenoDataDccColName = "Sample_ID",
                              protocolDataColNames = NULL, # each variable put here are removed from pData(object)
                              experimentDataColNames = experimentDataColNames)

saveRDS(Data, paste0(absolute_path_cur, "Owkin_Pilot_Data/GeoMx/", disease, "_raw/Data_", disease, ".rds"))


