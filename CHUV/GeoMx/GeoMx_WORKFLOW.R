disease = "lung"
# Read Data in with readNanoStringGeoMxSet(DCC, PKC) ----------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/01_GeoMx_Read_DCCs.R")
# GeoDiff count matrix ----------------------------------------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/02_Benchmarking_GeoDiff.R")
# NanostringQC template (GDR derivation, but do not filter on GDR) --------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/03_Nanostring_QC.R")
# standR normalization ----------------------------------------------------
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/04_standR_Normalization.R")




disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/02_Benchmarking_GeoDiff.R")

disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/02_Benchmarking_GeoDiff.R")





## Normalized SPE convert back to targetData for GDR ----------------------
# target_Data@assayData$logCPM <- assay(spe, "logCPM")
# Error in target_Data@assayData$logCPM <- assay(spe, "logCPM") : 
#   cannot add bindings to a locked environment
## Not possible to derive post norm GDR, because "NegGeoMean_" and "NegGeoSD_" are derived from before targetData,
## so now the only way is to add GDR to UMAP production for batch correction


# standR RLE -------------------------------------------------------------



# standR UMAP GDR --------------------------------------------------------
