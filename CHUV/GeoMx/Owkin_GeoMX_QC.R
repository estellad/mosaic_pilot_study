# env : geomx.yaml

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(NanoStringNCTools)
library(GeomxTools)
library(dplyr)
library(ggforce)
library(ggplot2)
library(scales)
library(GeoDiff)
library(Biobase)
library(reshape2) 
library(data.table)
library(knitr)
library(ggalluvial)

##### Load GeoMx object
load(file=snakemake@input[["obj"]]) # Data object
Data = shiftCountsOne(Data, useDALogic = TRUE)

##### Load QC parameters
QC_params =
    list(minSegmentReads = as.numeric(snakemake@params[["minSegmentReads"]]),
         percentTrimmed = as.numeric(snakemake@params[["percentTrimmed"]]),
         percentStitched = as.numeric(snakemake@params[["percentStitched"]]),
         percentAligned = as.numeric(snakemake@params[["percentAligned"]]),
         percentSaturation = as.numeric(snakemake@params[["percentSaturation"]]),
         minNegativeCount = as.numeric(snakemake@params[["minNegativeCount"]]),
         maxNTCCount = as.numeric(snakemake@params[["maxNTCCount"]]),
         minNuclei = as.numeric(snakemake@params[["minNuclei"]]),
         minArea = as.numeric(snakemake@params[["minArea"]]))

#################################
### Segment QC
#################################

##### Every ROI/AOI segment will be tested for:
# * Raw sequencing reads: segments with <1000 raw reads are removed.
# * % Aligned,% Trimmed, or % Stitched sequencing reads: segments below ~80% for one or more of these QC parameters are removed.
# * % Sequencing saturation ([1-deduplicated reads/aligned reads]%): segments below ~50% require additional sequencing to capture full sample diversity and are not typically analyzed until improved.
# * Negative Count: this is the geometric mean of the several unique negative probes in the GeoMx panel that do not target mRNA and establish the background count level per segment; segments with low negative counts (1-10) are not necessarily removed but may be studied closer for low endogenous gene signal and/or insufficient tissue sampling.
# * No Template Control (NTC) count: values >1,000 could indicate contamination for the segments associated with this NTC; however, in cases where the NTC count is between 1,000- 10,000, the segments may be used if the NTC data is uniformly low (e.g. 0-2 counts for all probes).
# * Nuclei: >100 nuclei per segment is generally recommended; however, this cutoff is highly study/tissue dependent and may need to be reduced; what is most important is consistency in the nuclei distribution for segments within the study.
# * Area: generally correlates with nuclei; a strict cutoff is not generally applied based on area.

Data <- setSegmentQCFlags(Data, qcCutoffs = QC_params)    

##### Collate QC Results
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

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         facet_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
    if(grepl(" ", fill_by)){
        colnames(assay_data)[which(colnames(assay_data) == fill_by)] = gsub(" ", "_", fill_by)
        fill_by = gsub(" ", "_", fill_by)
    }
    if(grepl(" ", facet_by)){
        colnames(assay_data)[which(colnames(assay_data) == facet_by)] = gsub(" ", "_", facet_by)
        facet_by = gsub(" ", "_", facet_by)
    }
    plt <- ggplot(assay_data,
                  aes_string(x = paste0("unlist(`", annotation, "`)"),
                             fill = fill_by)) +
                #   aes(x = !!!ensyms(annotation),
                #              fill = fill_by)) +
        geom_histogram(bins = 50) +
        geom_vline(xintercept = thr, lty = "dashed", color = "black") +
        theme_bw() + #guides(fill = "none") +
        facet_wrap(as.formula(paste("~", facet_by)), ncol = 1) +
        labs(x = annotation, y = "Segments, #", title = annotation)
    if(!is.null(scale_trans)) {
        plt <- plt +
            scale_x_continuous(trans = scale_trans)
    }
    plt
}

########### QC plots
pdf(snakemake@output[["segQCplt"]], width=10, height=5)

fill_by <- "slide name"
facet_by <- "segment"

QC_histogram(sData(Data), "Trimmed (%)", fill_by, facet_by, QC_params[["percentTrimmed"]])
QC_histogram(sData(Data), "Stitched (%)", fill_by, facet_by, QC_params[["percentStitched"]])
QC_histogram(sData(Data), "Aligned (%)", fill_by, facet_by, QC_params[["percentAligned"]])
QC_histogram(sData(Data), "Saturated (%)", fill_by, facet_by, QC_params[["percentSaturation"]]) +
    labs(title = "Sequencing Saturation (%)",
         x = "Sequencing Saturation (%)")
print(head(sData(Data)))
QC_histogram(sData(Data), "NTC", fill_by, facet_by, QC_params[["maxNTCCount"]], scale_trans = "log10")

if(snakemake@params[["plot_QC_area"]]){
    QC_histogram(sData(Data), "area", fill_by, facet_by, QC_params[["minArea"]], scale_trans = "log10")
}

if(snakemake@params[["plot_QC_nuclei"]]){
    QC_histogram(sData(Data), "nuclei", fill_by, facet_by, QC_params[["minNuclei"]])
}

##### Focus on negative count probe set (used to identify wells with low DNA material)

# calculate the negative geometric means for each module
negativeGeoMeans <- 
    esBy(negativeControlSubset(Data), 
         GROUP = "Module", 
         FUN = function(x) { 
             assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
         }) 
protocolData(Data)[["NegGeoMean"]] <- negativeGeoMeans

pkcs <- annotation(Data)
modules <- gsub(".pkc", "", pkcs)

# explicitly copy the Negative geoMeans from sData to pData

negCols <- paste0("NegGeoMean_", modules)
pData(Data)[, negCols] <- sData(Data)[["NegGeoMean"]]
for(ann in negCols) {
    plt <- QC_histogram(pData(Data), ann, fill_by, facet_by, QC_params[["minNegativeCount"]], scale_trans = "log10")
    print(plt)
}
dev.off()

# save(Data, file=snakemake@output[["postNegGeoMean"]])


# detatch neg_geomean columns ahead of aggregateCounts call
# pData(Data) <- pData(Data)[, !colnames(pData(Data)) %in% negCols]

print("####################################################")
# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(Data)$NTC),
      col.names = c("NTC Count", "# of Segments"),
      caption = "Non-Template Control count per segment")
kable(QC_Summary, caption = "QC Summary Table for each Segment")

## remove flagged segments
First_segm_filter = Data[, QCResults$QCStatus != "PASS"]
Data <- Data[, QCResults$QCStatus == "PASS"]


#################################
### Probe QC
#################################

# > Before we summarize our data into gene-level count data, we will remove low-performing probes.
# In short, this QC is an outlier removal process, whereby probes are either removed entirely from the study (global) or from specific segments (local). 
# >The QC applies to gene targets for which there are multiple distinct probes representing the count for a gene per segment.
# In WTA data, one specific probe exists per target gene; thus, Probe QC does not apply to the endogenous genes in the panel.
# Rather, it is performed on the negative control probes; there are multiple probes representing our negative controls, which do not target any sequence in the genome.
# These probes enable calculation of the background per segment and will be important for determining gene detection downstream.
# >After Probe QC, there will always remain at least one probe representing every gene target. In other words, Probe QC never removes genes from your data.
# >A probe is removed globally from the dataset if either of the following is true:
# > * the geometric mean of that probe’s counts from all segments divided by the geometric mean of all probe counts representing the target from all segments is less than 0.1
# > * the probe is an outlier according to the Grubb’s test in at least 20% of the segments
# > A probe is removed locally (from a given segment) if the probe is an outlier according to the Grubb’s test in that segment.
# > We do not typically adjust these QC parameters.

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
Data <- setBioProbeQCFlags(Data, 
                               qcCutoffs = list(minProbeRatio = snakemake@params[["minProbeRatio"]],
                                                percentFailGrubbs = snakemake@params[["GrubbOutPercent"]]), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(Data)[["QCFlags"]]


# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

## Exclude outliers

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
    subset(Data, 
           fData(Data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
               fData(Data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
Data <- ProbeQCPassed 

print("####################################################")
kable(qc_df, caption="Probe QC filter table:")

print("Feature number:")
dim(ProbeQCPassed)

print("Target number:")
length(unique(featureData(Data)[["TargetName"]]))
print("####################################################")

save(Data, file=snakemake@output[["unAggregated_object"]])
# collapse to targets

print("Aggregate counts from probe to genes")
print(dim(Data))
target_Data <- aggregateCounts(Data)
print(dim(target_Data))

if(length(unique(featureData(Data)[["TargetName"]])) != dim(target_Data)[1]){
    print("WARNING: Gene was drop by probe filtering")
}

#################################
### Limit Of Quantification (LOQ) QC
#################################

# > In addition to Segment and Probe QC, we also determine the limit of quantification (LOQ) per segment.
# The LOQ is calculated based on the distribution of negative control probes and is intended to approximate
# the quantifiable limit of gene expression per segment. Please note that this process is more stable in
# larger segments. Likewise, the LOQ may not be as accurately reflective of true signal detection rates
# in segments with low negative probe counts (ex: <2). The formula for calculating the LOQ in the ith segment is:

# LOQ_i=geomean(NegProbe_i)∗geoSD(NegProbe_i)^n 

# > We typically use 2 geometric standard deviations (n=2) above the geometric mean as the LOQ,
# which is reasonable for most studies. We also recommend that a minimum LOQ of 2 be used if the
# LOQ calculatedin a segment is below this threshold.

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_Data))
for(module in modules) {
    vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                   module)
    if(all(vars[1:2] %in% colnames(pData(target_Data)))) {
        LOQ[, module] <-
            pmax(snakemake@params[["minLOQ"]],
                 pData(target_Data)[, vars[1]] * 
                     pData(target_Data)[, vars[2]] ^ snakemake@params[["n_geoMean"]])
        pData(target_Data)[,paste0("LOQ_", module)] <- LOQ[, module]
    }
}


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

# Save detection rate information to pheno data
pData(target_Data)$GenesDetected <- 
    colSums(LOQ_Mat, na.rm = TRUE)
pData(target_Data)$GeneDetectionRate <-
    pData(target_Data)$GenesDetected / nrow(target_Data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_Data)$DetectionThreshold <- 
    cut(pData(target_Data)$GeneDetectionRate,
        breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
        labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

pdf(snakemake@output[["LOQQCplt"]], width=10, height=5)

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
options(repr.plot.width=6, repr.plot.height=4)
ggplot(pData(target_Data),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = segment)) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")

ggplot(pData(target_Data)) +
    geom_density(aes(x = GeneDetectionRate)) +
    geom_histogram(aes(x = GeneDetectionRate, fill=segment), binwidth=0.01) +
    theme_bw() +
    labs(x = "Gene Detection Rate")

ggplot(pData(target_Data),
       aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = `slide name`)) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate",
         y = "Segments, #",
         fill = "Segment Type")

ggplot(pData(target_Data)) +
    geom_density(aes(x = GeneDetectionRate)) +
    geom_histogram(aes(x = GeneDetectionRate, fill=`slide name`), binwidth=0.01) +
    theme_bw() +
    labs(x = "Gene Detection Rate")

dev.off()


##### Per well sensitivity plot
well_plt = as.data.frame(cbind(sData(target_Data)$Plate_ID,
                               sData(target_Data)$Well,
                               pData(target_Data)$GeneDetectionRate,
                               pData(target_Data)$`slide name`,
                               rownames(pData(target_Data)),
                               pData(target_Data)$`ROI Coordinate X`,
                               pData(target_Data)$`ROI Coordinate Y`,
                               pData(target_Data)$Patient
                               ))
well_first_filter = as.data.frame(cbind(sData(First_segm_filter)$Plate_ID,
                                        sData(First_segm_filter)$Well,
                                        rep(NA, length(sData(First_segm_filter)$Well)),
                                        pData(First_segm_filter)$`slide name`,
                                        rownames(pData(First_segm_filter)),
                                        pData(First_segm_filter)$`ROI Coordinate X`,
                                        pData(First_segm_filter)$`ROI Coordinate Y`,
                                        pData(First_segm_filter)$Patient
                                        ))

well_plt = rbind(well_plt, well_first_filter)
colnames(well_plt) = c("Plate", "Well", "GeneDetectionRate",
                       "slide_name", 'sample', "ROI_X_coord",
                       "ROI_Y_coord", "Patient")

well_plt$Filter = "PASS"
well_plt$Filter[which(well_plt$GeneDetectionRate < 0.1)] = "GeneDetectionRate"
well_plt$Filter[which(well_plt$sample %in% rownames(pData(First_segm_filter)))] = "FirstSegmentFilter"


print(head(well_plt))
print(table(well_plt$Filter))

well_plt$row = gsub("[0-9]", "", well_plt$Well)
well_plt$row = factor(well_plt$row, levels = rev(LETTERS[1:8]))
well_plt$col = as.numeric(gsub("[A-Z]", "", well_plt$Well))

well_plt$GeneDetectionRate = as.numeric(well_plt$GeneDetectionRate)
well_plt$ROI_X_coord = as.numeric(well_plt$ROI_X_coord)
well_plt$ROI_Y_coord = as.numeric(well_plt$ROI_Y_coord)


n_pat = length(unique(well_plt$Patient))
bigPalet = RColorBrewer::brewer.pal(10, "Paired")
multiple = floor(n_pat/length(bigPalet))
rest = n_pat %% length(bigPalet)
pat_palet = c(rep(bigPalet, multiple), bigPalet[1:rest])


n_slides = length(unique(well_plt$slide_name))
bigPalet = RColorBrewer::brewer.pal(10, "Paired")
multiple = floor(n_slides/length(bigPalet))
rest = n_slides %% length(bigPalet)
slide_palet = c(rep(bigPalet, multiple), bigPalet[1:rest])

pdf(file=snakemake@output[["well_plt"]], width=20, height=10)

print(ggplot(well_plt, aes(x=col, y=row, col=slide_name, fill = GeneDetectionRate)) +
    geom_point(size=7, shape=21, stroke=1.5) + theme_bw() + 
    geom_point(data = well_plt[which(well_plt$Filter == "FirstSegmentFilter"),], size=5, shape=4, col="darkred") +
    geom_point(data = well_plt[which(well_plt$Filter == "GeneDetectionRate"),], size=2, col="darkred") +
    scale_x_continuous(breaks=seq(1, 12, 1)) +
    scale_color_manual(values = slide_palet) +
    scale_fill_gradient2(limits = c(0, 1)) +
    facet_wrap(".~Plate")
)

print(ggplot(well_plt, aes(x=col, y=row, col=Patient, fill = GeneDetectionRate)) +
    geom_point(size=7, shape=21, stroke=1.5) + theme_bw() + 
    geom_point(data = well_plt[which(well_plt$Filter == "FirstSegmentFilter"),], size=5, shape=4, col="darkred") +
    geom_point(data = well_plt[which(well_plt$Filter == "GeneDetectionRate"),], size=2, col="darkred") +
    scale_x_continuous(breaks=seq(1, 12, 1)) +
    scale_color_manual(values = pat_palet) +
    scale_fill_gradient2(limits = c(0, 1)) +
    facet_wrap(".~Plate")
)

print(ggplot() +
    geom_point(data = well_plt[which(well_plt$Filter != "FirstSegmentFilter"),],
               aes(x=ROI_X_coord, y=ROI_Y_coord, fill = GeneDetectionRate, col=Filter),
               shape=21, size=3) +
    scale_color_manual(values=c("PASS"="royalblue", "GeneDetectionRate"="darkorange")) +
    theme_bw() +
    geom_point(data = well_plt[which(well_plt$Filter == "FirstSegmentFilter"),],
               aes(x=ROI_X_coord, y=ROI_Y_coord, fill = GeneDetectionRate),
               size=3, shape=4, col="red") +
    scale_fill_gradient2(limits = c(0, 1)) + theme_bw() +
    facet_wrap(".~Plate")
)

print(ggplot() +
    geom_point(data = well_plt[which(well_plt$Filter != "FirstSegmentFilter"),],
               aes(x=ROI_X_coord, y=ROI_Y_coord, fill = GeneDetectionRate, col=Filter),
               shape=21, size=3) +
    scale_color_manual(values=c("PASS"="royalblue", "GeneDetectionRate"="darkorange")) +
    theme_bw() +
    geom_point(data = well_plt[which(well_plt$Filter == "FirstSegmentFilter"),],
               aes(x=ROI_X_coord, y=ROI_Y_coord, fill = GeneDetectionRate),
               size=3, shape=4, col="red") +
    scale_fill_gradient2(limits = c(0, 1)) + theme_bw() +
    facet_wrap(".~Patient")
)


dev.off()
write.table(well_plt, file=snakemake@output[["qc_table"]], sep="\t", quote=F, row.names=FALSE)

### We first filter out segments with exceptionally low signal.
print("Filtering segments based on GeneDetectionRate threshold")
print("Shape before:")
print(dim(pData(target_Data)))
target_Data <-
    target_Data[, pData(target_Data)$GeneDetectionRate >= snakemake@params[["min_gene_DetectionRate"]]]
print("Shape after:")
print(dim(pData(target_Data)))


# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_Data)]
fData(target_Data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_Data)$DetectionRate <-
    fData(target_Data)$DetectedSegments / nrow(pData(target_Data))

# Gene of interest detection table
goi <- snakemake@params[["goi"]]
goi_df <- data.frame(
    Gene = goi,
    Number = fData(target_Data)[goi, "DetectedSegments"],
    DetectionRate = percent(fData(target_Data)[goi, "DetectionRate"]))


print("####################################################")
kable(goi_df, capton="Detection rate of genes of interest")
print("####################################################")

# We will graph the total number of genes detected in different percentages of segments.
# Based on the visualization below, we can better understand global gene detection in
# our study and select how many low detected genes to filter out of the dataset. 
# Gene filtering increases performance of downstream statistical tests and improves
# interpretation of true biological signal.

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
    unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                  function(x) {sum(fData(target_Data)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_Data))
rownames(plot_detect) <- plot_detect$Freq

pdf(snakemake@output[["detecRatePlt"]], width=10, height=5)

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
              vjust = 1.6, color = "black", size = 4) +
    scale_fill_gradient2(low = "orange2", mid = "lightblue",
                         high = "dodgerblue3", midpoint = 0.65,
                         limits = c(0,1),
                         labels = scales::percent) +
    theme_bw() +
    scale_y_continuous(labels = scales::percent, limits = c(0,1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = "Threshold of % of segment with gene signal above DetectionRate",
         y = "Genes Detected, % of Panel > LOQ")

dev.off()

# Subset to target genes detected in at least X% (default = 10%) of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_Data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)


print("Filtering genes based on GeneDetectionRate threshold")
print("Shape before:")
print(dim(fData(target_Data)))
target_Data <- 
    target_Data[fData(target_Data)$DetectionRate >= snakemake@params[["min_segment_DetectionRate"]] |
                        fData(target_Data)$TargetName %in% neg_probes, ]
print("Shape after:")

print(dim(fData(target_Data)))

save(target_Data, file=snakemake@output[["QC_passed_object"]])


#Let’s re-plot the Sankey diagram showing our current working dataset.
#This is now a dataset that no longer contains segments flagged by Segment QC or that have low gene detection rates.
png(snakemake@output[["sankey_postQC"]], 700, 700)

alu = as.data.frame(pData(target_Data))
alu = alu[which(!is.na(alu$segment)),]
if(snakemake@params[["Alluvial_Axis3"]] == 0 & snakemake@params[["Alluvial_Axis4"]] == 0){
  ggplot(alu,
        aes(axis1 = !!sym(snakemake@params[["Alluvial_Axis1"]]),
            axis2 = !!sym(snakemake@params[["Alluvial_Axis2"]]))) +
    geom_alluvium(aes(fill = !!sym(snakemake@params[["Alluvial_Fill"]])), width = 1/12) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") + 
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(snakemake@params[["Alluvial_Axis1"]],
                                snakemake@params[["Alluvial_Axis2"]]), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme_bw()  
}else if(snakemake@params[["Alluvial_Axis4"]] == 0){
  ggplot(alu,
        aes(axis1 = !!sym(snakemake@params[["Alluvial_Axis1"]]),
            axis2 = !!sym(snakemake@params[["Alluvial_Axis2"]]),
            axis3 = !!sym(snakemake@params[["Alluvial_Axis3"]]))) +
    geom_alluvium(aes(fill = !!sym(snakemake@params[["Alluvial_Fill"]])), width = 1/12) +
    geom_stratum(width = 1/12, fill = "white", color = "grey") + 
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c(snakemake@params[["Alluvial_Axis1"]],
                                snakemake@params[["Alluvial_Axis2"]],
                                snakemake@params[["Alluvial_Axis3"]]), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    theme_bw()  
}else{
  ggplot(alu,
      aes(axis1 = !!sym(snakemake@params[["Alluvial_Axis1"]]),
          axis2 = !!sym(snakemake@params[["Alluvial_Axis2"]]),
          axis3 = !!sym(snakemake@params[["Alluvial_Axis3"]]),
          axis4 = !!sym(snakemake@params[["Alluvial_Axis4"]]))) +
  geom_alluvium(aes(fill = !!sym(snakemake@params[["Alluvial_Fill"]])), width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c(snakemake@params[["Alluvial_Axis1"]],
                              snakemake@params[["Alluvial_Axis2"]],
                              snakemake@params[["Alluvial_Axis3"]],
                              snakemake@params[["Alluvial_Axis4"]]), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_bw()   
}
dev.off()
