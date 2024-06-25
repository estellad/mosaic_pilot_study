args <- commandArgs(trailingOnly = TRUE)

disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/05_Baye_Spot_raw/spot_params.R")

i = as.numeric(args[1])

print(i)
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
save_bs_path <- paste0(baye_savepath, foldername, "/", save_names[i], "/")
sce <- readRDS(file.path(save_bs_path, paste0(save_names[i], "_PCAed.rds")))

assay(sce, "log1p") <- log1p(counts(sce))

# spatialCluster ----------------------------------------------------------
print("spot clustering")
set.seed(123)
sce = spatialCluster(sce, use.dimred = "PCA", q = nclus[i], nrep = 50000, d = 50, burn.in = 1000, init.method = spcl.initmethod[i])
saveRDS(sce, file.path(save_bs_path, paste0(save_names[i], "_baye_clustered.rds")))


###########################################################################
#                    (spatialEnhance on full genes)                       #
###########################################################################

# spatialEnhance Full PCs with all genes ----------------------------------
# sce <- readRDS(file.path(save_bs_path, paste0(save_names[i], "_baye_clustered.rds")))
print("spot -> subspot getting PCs all genes")
set.seed(123)
sce_enhanced <- spatialEnhance(
  sce,
  q = nclus[i], platform = "Visium",
  use.dimred = "PCA",
  model = "t", gamma = 3,
  jitter_prior = 0.3, jitter_scale = 0, adapt.before = 0,
  nrep = 100000, burn.in = 10000,
  save.chain=TRUE,
  chain.fname=file.path(save_bs_path, paste0(save_names[i], "_enhanced_all_genes.h5")),
  verbose = TRUE, cores = 2
)

saveRDS(sce_enhanced, file.path(save_bs_path, paste0(save_names[i], "_baye_clustered_all_enhanced.rds"))) 


# ####################################### Enhance Feature #################################
# Helper enhance a subset --------------------------------------------------
print("loading helper on enhance a subset")
enhance_a_subset <- function(sce_enhanced, sce){
  print(system.time(sce_enhanced <- enhanceFeatures(sce_enhanced, sce,
                                                    model="lm",
                                                    assay.type = "log1p",
                                                    feature_names=rownames(sce),
                                                    nrounds=0)))
  
  
  if(min(exp(assay(sce_enhanced, "log1p") - 1), na.rm = TRUE) >= 0){# sanity: check if this will have negative
    assay(sce_enhanced, "raw_subspot") <-  exp(assay(sce_enhanced, "log1p") - 1)
  }else{
    assay(sce_enhanced, "raw_subspot") <-  exp(assay(sce_enhanced, "log1p"))
  }
  
  return(sce_enhanced)
}

# Enhance expr the full --------------------------------------------------
print("clustered_enhanced_expr on all genes")
sce_enhanced <- enhance_a_subset(sce_enhanced, sce)
saveRDS(sce_enhanced, file.path(save_bs_path, paste0(save_names[i], "_baye_clustered_all_enhanced_expr.rds")))





