# lung GDR spe 18677 123 -> 122
# breast GDR spe 18677 119 -> 117
# dlbcl GDR spe 18677 143 -> 137

# disease = "lung"
# disease = "breast"
disease = "dlbcl"
source(paste0(absolute_path_cur, "env/ydong/Owkin_Pilot/Code/GeoMx/Manuscript/GeoMx_init.R"))
library(ggspavis)

save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx")
spe <- readRDS(file.path(save_path, paste0("GDR_spe_", disease,".rds")))

CD <- data.frame(colData(spe))
CD <- janitor::clean_names(CD)

CD <- CD %>%
  rename(sample_id2 = sample_id, 
         sample_id = sample_id_2)

head(CD)
colData(spe) <- NULL
colData(spe) <- as(CD, "DFrame") 


# QC
## Gene-level
spe <- addPerROIQC(spe, rm_genes = TRUE)
# dim(spe) # lung 18677  122 # breast 18677   117

if(nrow(spe) == 18677) 
  print("No gene removed")


## AOI-level 
spe$log_lib_size <- log(spe$lib_size)
spe$log_area <- log(spe$area)

plotROIQC(spe, 
          x_axis = "gene_detection_rate", y_axis = "log_lib_size",
          x_lab = "Gene Detection Rate", y_lab = "Log Lib Size",
          color = slide_name) |
  plotROIQC(spe, 
            x_axis = "gene_detection_rate", y_axis = "log_lib_size",
            x_lab = "Gene Detection Rate", y_lab = "Log Lib Size",
            color = patient_id) |
  plotROIQC(spe, 
            x_axis = "gene_detection_rate", y_axis = "log_lib_size",
            x_lab = "Gene Detection Rate", y_lab = "Log Lib Size",
            color = cell_fraction)

# No lung samples need to be further removed

## Normalization
## scran Normalization
library(scran)
spe <- computeSumFactors(spe, assay.type="counts")
data_norm_pre <- sweep(assays(spe)$counts, 2, spe$sizeFactor,'/')
assays(spe, withDimnames=FALSE)$scrannormcounts <- log(data_norm_pre + 1)

## logNormCounts - mean count normalization
assays(spe, withDimnames=FALSE)$meannormcounts <- assay(logNormCounts(spe, assay.type = "counts"), "logcounts")

## TMM
spe_tmm <- geomxNorm(spe, method = "TMM", log = FALSE)
assays(spe, withDimnames=FALSE)$tmmcounts <- assay(spe_tmm, "logcounts")

## Upperquartile (Q3)
spe_upperquatile <- geomxNorm(spe, method = "upperquartile", log = FALSE)
assays(spe, withDimnames=FALSE)$upperquartile <- assay(spe_upperquatile, "logcounts")

## Quantile
assays(spe, withDimnames=FALSE)$quantile <- normalize.quantiles(assay(spe, "log1p")) # preprocessCore::normalize.quantile on top of log1p

## scTransform
assay(spe, "scTransform", withDimnames = FALSE) <- sctransform::vst(counts(spe), min_cells = 0)$y

# > assays(spe)
# List of length 9
# names(10): counts logcounts log1p scrannormcounts meannormcounts tmmcounts upperquartile quantile scTransform 

save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoMx_test_norm/")
saveRDS(spe, file.path(save_path, paste0(disease, ".rds")))

# spe <- readRDS(file.path(save_path, paste0(disease, ".rds")))

# TODO: add in GeoDiff rds -----------------
test <- readRDS(paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx/GeoDiff_", disease, "_countmat.rds")) # lung 18676 101 # breast 18676 110 # dlbcl 115
colnames(test) <- str_replace(colnames(test), ".dcc", "")

common_samples <- intersect(colnames(spe), colnames(test))

spe <- spe[, colnames(spe) %in% common_samples] 
test <- test[, common_samples]

# check colnames are the same order
all(colnames(spe) == colnames(test))

setdiff(rownames(spe), rownames(test)) # "NegProbe-WTX", remove from normalization benchmarking (might add back for decon)

spe <- spe[rownames(spe) %in% rownames(test), ]

# check rownames are the same order
all(rownames(spe) == rownames(test))

assay(spe, "GeoDiff", withDimnames = FALSE) <- test

# > assays(spe)
# List of length 10
# names(10): counts logcounts log1p scrannormcounts meannormcounts tmmcounts upperquartile quantile scTransform GeoDiff


saveRDS(spe, file.path(save_path, paste0(disease, "_geodiff_included.rds")))














