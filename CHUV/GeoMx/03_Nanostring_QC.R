library(NanoStringNCTools)
library(GeomxTools)

Data <- readRDS(paste0(absolute_path_cur, "Owkin_Pilot_Data/GeoMx/", disease, "_raw/Data_", disease, ".rds"))

library(knitr)
pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

Data <- shiftCountsOne(Data, useDALogic = TRUE) # lung: 18815 123
                                                # breast: 18815 119
                                                # dlbcl: 18815 143

#########################################################################
#                               Segment QC                              #
#########################################################################
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
Data <-
  setSegmentQCFlags(Data, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(Data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(Data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(Data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(Data)[, negCols] <- sData(Data)[["NegGeoMean"]]


# detatch neg_geomean columns ahead of aggregateCounts call
pData(Data) <- pData(Data)[, !colnames(pData(Data)) %in% negCols]

kable(QC_Summary, caption = "QC Summary Table for each Segment") # lung: 1 low align
                                                                 # breast: 2 low saturation
                                                                 # dlbcl: 1 low reads, 2 low saturation, 4 low area

Data <- Data[, QCResults$QCStatus == "PASS"] # lung: 18815 122
                                             # breast: 18815 117
                                             # dlbcl: 18815 137


#########################################################################
#                                Probe QC                               #
#########################################################################
# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
Data <- setBioProbeQCFlags(Data, 
                           qcCutoffs = list(minProbeRatio = 0.1,
                                            percentFailGrubbs = 20), 
                           removeLocalOutliers = TRUE)

ProbeQCResults <- fData(Data)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

# > qc_df
## lung
# Passed Global Local
# 1  18798      0    17
## breast
# Passed Global Local
# 1  18795      0    20
## dlbcl
# Passed Global Local
# 1  18792      0    23

# Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(Data, 
         fData(Data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(Data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)
#> Features  Samples 
#>  18815      122 
#> Features  Samples 
#>  18815      117 
#> Features  Samples 
#>  18815      137 
Data <- ProbeQCPassed 


# aggregateCounts() -------------------------------------------------------
# Check how many unique targets the object has
length(unique(featureData(Data)[["TargetName"]]))
#> [1] 18677

# collapse to targets
target_Data <- aggregateCounts(Data)
dim(target_Data)
#> Features  Samples 
#>    18677      122
#>    Features  Samples 
#>    18677      117 


# LOQ ---------------------------------------------------------------------
# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_Data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_Data)[, vars[1]] * 
             pData(target_Data)[, vars[2]] ^ cutoff)
  }
}
pData(target_Data)$LOQ <- LOQ


LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_Data)$Module == module
  Mat_i <- t(esApply(target_Data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_Data)$TargetName, ]

# Segment GDR ------------------------------------------------------------
# Save detection rate information to pheno data
pData(target_Data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Data)$GeneDetectionRate <-
  pData(target_Data)$GenesDetected / nrow(target_Data)

# # GDR filtering out
# test <-
#   target_Data[, pData(target_Data)$GeneDetectionRate >= .1]
# 
# dim(test) 
# lung: post GDR 10% has 47 samples left
# breast: post GDR 10% has 115 samples left
# dlbcl: post GDR 10% has 132 samples left


# Convert to SPE for standR -----------------------------------------------
# SPE Coercion ------------------------------------------------------------
spe <- as.SpatialExperiment(target_Data, normData = "exprs", forceRaw = TRUE); dim(spe) # lung: 18677  122
                                                                                        # breast: 18677  117
                                                                                        # dlbcl: 18677 137

# Name count matrices -----------------------------------------------------
library(SpatialExperiment)
names(assays(spe))[which(names(assays(spe)) == "GeoMx")] <- "counts"
assay(spe, "logcounts", withDimnames = FALSE) <- edgeR::cpm(counts(spe), log = TRUE)
assay(spe, "log1p", withDimnames = FALSE) <- log1p(counts(spe))
assays(spe)

# Sanity checks
rownames(counts(spe))[which(grepl("Neg", rownames(counts(spe))))] # There is one neg control probe for spatial decon

colnames(spe) <- gsub(".dcc", "", colnames(spe))

if(disease %in% c("breast", "lung")){
  spe$sample_ids <- paste0(spe$Section_ID, "_", spe$Patient_pathID)
  table(spe$sample_ids)
  # # lung
  # L1_1_0PSV L2_1_0WMU L3_1_1G73 L3_3_1G73 L4_3_1GA2 
  # 25        24        23        24        26 
  # B1_1_OPHI B1_3_OPHI B2_1_1256 B3_1_1GVR B4_1_1FHZ 
  # 24        24        22        23        24 
}else if(disease == "dlbcl"){
  table(spe$Patient)
  # D1 D2 D3 D4 D5 D6 
  # 24 23 24 21 26 19 
}

save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx")
saveRDS(spe, file.path(save_path, paste0("GDR_spe_", disease,".rds"))) # lung GDR spe 18677 123 -> 122
                                                                       # breast GDR spe 18677 119 -> 115
                                                                       # dlbcl GDR spe 18677 147 -> 137






