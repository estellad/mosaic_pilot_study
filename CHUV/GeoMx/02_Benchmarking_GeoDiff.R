# BiocManager::install("GeoDiff")

library(GeoDiff)
library(dplyr)
library(ggplot2)
library(NanoStringNCTools)
library(GeomxTools)
library(Biobase)
library(reshape2)


Data <- readRDS(paste0(absolute_path_cur, "Owkin_Pilot_Data/GeoMx/", disease, "_raw/Data_", disease, ".rds"))

# Background --------------------------------------------------------------
featureType(Data)
#> [1] "Probe"

paste("## of Negative Probes:", sum(fData(Data)$Negative))


Data <- fitPoisBG(Data)
# summary(fData(Data)$featfact[fData(Data)$Negative])


set.seed(123)
Data_diag <- diagPoisBG(Data)
notes(Data_diag)$disper                 # lung: 2.660215 -> there is batch effect
                                        # breast: 1.23529 -> there is no batch effect already, but also try remove it
                                        # dlbcl: 1.203253 -> there is no batch effect already, but also try remove it
which(assayDataElement(Data_diag, "low_outlier") == 1, arr.ind = TRUE)
which(assayDataElement(Data_diag, "up_outlier") == 1, arr.ind = TRUE)


# Batch effect ------------------------------------------------------------
Data <- fitPoisBG(Data, groupvar = "slide name")
set.seed(123)
Data_diag <- diagPoisBG(Data, split = TRUE)
notes(Data_diag)$disper_sp       # lung: 1.873139, disper is now < 2, no more batch effect
                                 # breast: 1.218413
                                 # dlbcl: 1.058326


# Aggregate ---------------------------------------------------------------
all0probeidx <- which(rowSums(exprs(Data))==0)
if (length(all0probeidx) > 0) {
  Data <- Data[-all0probeidx, ]
}
Data <- aggreprobe(Data, use = "cor") # lung: 18815 123 (18676 genes + 139 RTS neg probes)
                                      # breast: 18815 119 (18676 genes + 139 RTS neg probes)
                                      # dlbcl: 18815 143 (18676 genes + 139 RTS neg probes)


# Target QC ---------------------------------------------------------------
## Score test
# filter the data to only targets above background, using a suggested pvalue threshold of 1e-3
# Not excluding anything here
Data <- BGScoreTest(Data)
sum(fData(Data)[["pvalues"]] < 1e-3, na.rm = TRUE) # lung 13995 
                                                   # breast 15791 
                                                   # dlbcl 17091 

## Estimate the size factor
set.seed(123)
Data <- fitNBth(Data, split = TRUE)
features_high <- rownames(fData(Data))[fData(Data)$feature_high_fitNBth == 1] # lung: 1500
                                                                              # breast: 1500
                                                                              # dlbcl: 1500
length(features_high)

# To estimate the signal size factor, we use the fit negative binomial threshold function. 
# This size factor represents technical variation between ROIs like sequencing depth
bgMean <- mean(fData(Data)$featfact, na.rm = TRUE) 
# notes(Data)[["threshold"]] # lung: 1253.047
                             # breast: 1928.293
                             # dlbcl: 2236.028
# bgMean                     # lung: 1782.245
                             # breast: 1739.079
                             # dlbcl: 1130.784


# SampleQC ---------------------------------------------------------------
ROIs_high <- sampleNames(Data)[which((quantile(fData(Data)[["para"]][, 1],
                                                 probs = 0.90, na.rm = TRUE) -
                                          notes(Data)[["threshold"]])*Data$sizefact_fitNBth>2)] # lung: 123 -> 101 samples are roi high
                                                                                                # breast: 119 -> 110 samples are roi high
                                                                                                # dlbcl: 143 -> 115 samples are roi high

# get only biological probes, so 139 RTS are also removed from GeoDiff norm benchmarking.
posdat <- Data[-which(fData(Data)$CodeClass == "Negative"), ]
posdat <- exprs(posdat)
features_all <- rownames(posdat) # 18676


# Norm, not split by slides -----------------------------------------------
names(assayData(Data))
#> [1] "exprs"
set.seed(123)
Data <- fitPoisthNorm(object = Data,
                      features_high = features_high, # 1500
                      ROIs_high = ROIs_high,
                      threshold_mean = bgMean,
                      sizescalebythreshold = TRUE)

names(assayData(Data))

# View(assayData(Data)$normmat)


# Norm, split by slides ---------------------------------------------------
set.seed(123)
Data <- fitPoisthNorm(object = Data,
                      split = TRUE,
                      features_high = features_high,
                      ROIs_high = ROIs_high,
                      threshold_mean = bgMean,
                      sizescalebythreshold = TRUE)
#> The results are based on stored `groupvar`, slide name
#> probe finished
#> Iteration = 1, squared error = 0.0312931060333072
#> probe finished
#> Iteration = 2, squared error = 5.04355329507063e-06
#> Model converged.
#> probe finished
#> Iteration = 1, squared error = 0.0481496582408976
#> probe finished
#> Iteration = 2, squared error = 1.55472736038387e-06
#> Model converged.

names(assayData(Data))
#> [1] "exprs"       "normmat0_sp" "normmat"     "normmat_sp"  "normmat0"
# normmat-sp - normalization after iteration 2, with consideration of slide name batch effect

# # Compare norm methods ----------------------------------------------------
# # get only biological probes
# posdat <- Data[-which(fData(Data)$CodeClass == "Negative"), ]
# posdat <- exprs(posdat)
# 
# quan <- sapply(c(0.75, 0.8, 0.9, 0.95), function(y)
#   apply(posdat, 2, function(x) quantile(x, probs = y)))
# 
# corrs <- apply(quan, 2, function(x) cor(x, Data$sizefact_fitNBth))
# names(corrs) <- c(0.75, 0.8, 0.9, 0.95)
# 
# corrs
# #>      0.75       0.8       0.9      0.95 
# #> 0.9807486 0.9843474 0.9840221 0.9762114
# 
# quan75 <- apply(posdat, 2, function(x) quantile(x, probs = 0.75))
# 
# # features_all <- rownames(posdat)
# # TODO
# 
# 
# norm_dat_backqu75 <- sweep(posdat[, ROIs_high], 2,
#                            (Data[, ROIs_high]$sizefact * bgMean),
#                            FUN = "-") %>%
#   sweep(., 2, quan75[ROIs_high], FUN = "/") %>%
#   pmax(., 0) %>%
#   `+`(., 0.01) %>%
#   log2()
# 
# 
# 
# # Q3 plot -----------------------------------------------------------------
# dat_plot <- cbind(pData(Data)[ROIs_high, c("slide name", "Indication")],
#                   t(norm_dat_backqu75[features_high, ]))
# 
# dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)
# 
# dat_plot <- melt(dat_plot, id.vars = c("ROI_ID", "slide name", "Indication"))
# colnames(dat_plot) <- c("ROI_ID", "slide_name", "Indication", "variable", "value")
# 
# 
# ggplot(dat_plot, aes(x = value)) +
#   geom_density(aes(fill = slide_name, group = ROI_ID, alpha = 0.01)) +
#   facet_wrap(~slide_name) +
#   ggtitle("Q3 Normalization")+
#   labs(x = "Q3 Normalized Value (log2)")
# 
# 
# 
# # GeoDiff plot ------------------------------------------------------------
# annot <- pData(Data)
# 
# dat_plot <- cbind(annot[ROIs_high, c("slide name", "Indication")],
#                   t(assayDataElement(Data[features_high, ROIs_high], "normmat_sp")))
# dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)
# dat_plot <- melt(dat_plot, id.vars = c("ROI_ID", "slide name", "Indication"))
# colnames(dat_plot) <- c("ROI_ID", "slide_name", "Indication", "variable", "value")
# 
# ggplot(dat_plot, aes(x = value)) +
#   geom_density(aes(fill = slide_name, group = ROI_ID, alpha = 0.01)) +
#   facet_wrap(~`slide_name`) +
#   ggtitle("Poisson threshold normalization")+
#   labs(x = "Poisson Threshold Normalized Value (log2)")
# 
# 
# annot <- pData(Data)
# 
# dat_plot <- cbind(annot[ROIs_high, c("slide name", "Indication")],
#                   t(assayDataElement(Data[features_high, ROIs_high], "normmat")))
# dat_plot <- cbind(dat_plot, ROI_ID = ROIs_high)
# dat_plot <- melt(dat_plot, id.vars = c("ROI_ID", "slide name", "Indication"))
# colnames(dat_plot) <- c("ROI_ID", "slide_name", "Indication", "variable", "value")
# 
# ggplot(dat_plot, aes(x = value)) +
#   geom_density(aes(fill = slide_name, group = ROI_ID, alpha = 0.01)) +
#   facet_wrap(~`slide_name`) +
#   ggtitle("Poisson threshold normalization")+
#   labs(x = "Poisson Threshold Normalized Value (log2)")


# GeoDiff -----------------------------------------------------------------
# View(assayDataElement(Data, "normmat_sp")) # lung: 18815 123 samples -> save for benchmarking? no too much NAs everywhere
# 
# mat <- assayDataElement(Data, "normmat_sp") # all RTS genes have NA post GeoDiff norm
# length(which(is.na(colSums(mat)))) # lung: 123 NAs
# 
# mat <- assayDataElement(Data[features_all, ], "normmat_sp")
# length(which(is.na(colSums(mat)))) # lung: 22 NAs, only ROI high is kept for normalization

# So the final candidate for saving is, just biological probes and ROIs past high threshold
mat <- assayDataElement(Data[features_all, ROIs_high], "normmat_sp")
dim(mat)            # lung: 18676   101
                    # breast: 18676   110
                    # dlbcl: 18676   115
save_path <- paste0(absolute_path_cur, "Owkin_Pilot_Intermediate/GeoMx")
saveRDS(mat, file.path(save_path, paste0("GeoDiff_", disease, "_countmat.rds")))









