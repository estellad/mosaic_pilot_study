# env : visu.yaml

log = file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(R.utils)
library(ggforce)
library(ggplot2)
library(data.table)
library(reshape2)
library(cowplot)
library(patchwork)
library(ggpubr)
library(tibble)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(readxl)

plotHeatmap = function(res.ccp, annot, clustering_col, annotation_columns,out_hmap_file = "~/Heatmap.png", 
                       show_legend=NULL, legd=TRUE, plt_width=1200, plt_height=800, width=unit(8, "cm"), 
                       height=unit(8,"cm"), titl="", rownames_hmap=F, heatmp=NULL){
  
  annot = annot[,which(colnames(annot) %in% c(clustering_col, names(annotation_columns)) )]
  
  k = max(as.numeric(annot[,clustering_col]))
  co.mat = res.ccp[[k]][[1]]
  colnames(co.mat) = rownames(co.mat) = names(res.ccp[[k]][[3]])
  clust = res.ccp[[k]]$consensusTree
  
  cols_interest = c(
    clustering_col,
    names(annotation_columns)[names(annotation_columns) %in% colnames(annot)]
  )
  
  annotation_columns[[clustering_col]] = setNames(res.ccp[[k]]$clrs[[3]], 1:k)
  
  df = annot[, cols_interest]
  
  if(is.null(show_legend)){
    show_legend = rep(T, ncol(df))
  }else{
    show_legend = colnames(df) %in% show_legend
    names(show_legend) = colnames(df)
  }
  top_annot = HeatmapAnnotation(
    df=df,
    col=annotation_columns,
    annotation_name_side = "left",
    na_col="ghostwhite",
    show_legend = show_legend
  )
  
  if(!is.null(heatmp)){
    heatmp = t(scale(t(heatmp)))
    H = ComplexHeatmap::Heatmap(
      heatmp,
      top_annotation=top_annot, 
      cluster_columns=clust,
      clustering_method_rows = "ward.D2",
      clustering_distance_rows = "pearson",
      show_column_names=T,
      show_heatmap_legend = legd,
      col = colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
      height = height,
      width=width,
      show_row_names = rownames_hmap,
      column_title = titl)
  }else{
    H = ComplexHeatmap::Heatmap(
      co.mat,
      top_annotation=top_annot, 
      cluster_columns=clust,
      cluster_rows=clust,
      show_column_names=T,
      show_heatmap_legend = legd,
      show_row_names=rownames_hmap,
      col = colorRamp2(c(0, 1), c("white", "blue")), 
      height = height,
      width=width,
      column_title = titl)
  }
  
  # options(repr.plot.width = plt_width, repr.plot.height = plt_height)
  png(file=out_hmap_file, width = plt_width, height = plt_height)
  ComplexHeatmap::draw(H, annotation_legend_side = "right")
  dev.off()
}

clusterSamples = function(d, pI, maxK, plot){
  # median absolute deviation over samples per transcripts 
  
  cat("Number of genes took in account for clustering : ", nrow(d), "\n")
  print("Before as.dist function, print the class of d")
  print(class(d))
  print("Get number of NA in d")
  print(sum(is.na(d)))  
  d = as.dist(1 - cor(d, method = "pearson"))
  
  pdf(file=plot, width=10, height=10)
  res.ccp = ConsensusClusterPlus(
    d=d,
    maxK=maxK,
    pItem=pI,
    pFeature=1,
    reps=1000,
    seed=42,
    plot=NULL,
    verbose=FALSE,
    innerLinkage="ward.D2",
    finalLinkage="ward.D2",
    clusterAlg="hc", # pam / km
    distance="pearson" # spearman euclidian, binary, maximum, caberra, minkowski
  )
  dev.off()
  return(res.ccp)
} 

generate_annot_colors = function(annotation, cols_to_skip){
  annotation = annotation[, which(!colnames(annotation) %in% cols_to_skip)]
  cols = list()
  i=1
  bigColSet = c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(10, "Set3"))
  for(ann in colnames(annotation)){
    print(ann)
    annotation[which(annotation[,ann] == ""), ann] = NA
    if(is.numeric(annotation[,ann])){
      col = colorRamp2(c(quantile(annotation[,ann], probs = seq(.1, .9, by = .1), na.rm=T)[1], 
                         quantile(annotation[,ann], probs = seq(.1, .9, by = .1), na.rm=T)[9]), 
                       c("khaki1", "red3"))
      cols[[i]] = col
    }else if(is.character(annotation[,ann]) | is.factor(annotation[,ann])){
      col = bigColSet[seq_len(length(unique(na.omit(annotation[,ann]))))]
      names(col) = unique(na.omit(annotation[,ann]))
      cols[[i]] = col
    }else{
      print("Error, not numeric nor character annotation")
      print(ann)
    }
    i=i+1
  }
  names(cols) = colnames(annotation)
  return(cols)
}


run_clustering = function(d, an, minK, maxK, output_dir, pI=1){
  ## Clustering ~K
  print(names(an))
  for(clmn in names(an)){
    if(length(unique(an[,clmn]))==1){
      an = an[,!names(an) == clmn]
    }
  }
  print("Names of the metadata columns with more than one unique valeu suppled to run_clustering")
  print(names(an))
  print("Now execute clusterSamples function")
  resCCP = clusterSamples(d=d, pI=pI, maxK=min(maxK, ncol(d)-1),
                          plot=paste0(output_dir, "/K_choice.pdf"))
  
  for(k in minK:min(c(maxK, ncol(d)-1))){
    cs = resCCP[[k]]$consensusClass
    clusters = as.factor(cs)
    an$Cluster = clusters[match(rownames(an), names(clusters))]
    cols = generate_annot_colors(an, c("Sample", "Cluster"))
    
    if(nrow(d) > 20){
      rownames_hmap=F
    }else{
      rownames_hmap=T
    }
    plotHeatmap(resCCP, an, "Cluster", cols, out_hmap_file = paste0(output_dir, "/Clustering_k", k, "_", "CCF.png"),
                legd=F,plt_width=1500, plt_height=1000, width=unit(30, "cm"), height=unit(20, "cm"),
                titl=paste0("#Genes=", nrow(d)), rownames_hmap=rownames_hmap)
    
    plotHeatmap(resCCP, an, "Cluster", cols, out_hmap_file = paste0(output_dir, "/Clustering_k", k, "_", "exp.png"),
                legd=F, heatmp = d,plt_width=1500, plt_height=1000, width=unit(30, "cm"), height=unit(20, "cm"),
                titl=paste0("#Genes=", nrow(d)), rownames_hmap=rownames_hmap)
    an[,paste0("Cluster_k_", k)] = an$Cluster
    an$Cluster = NULL
  }
  return(an)
}



#expr = as.data.frame(fread("results/CHUV2_2_comparison/exp_matrices/CHUV2_2_comparison_quartil3_norm.tsv.gz", sep="\t"))
# read the input expression matrix - this input can be one of many expression matrices
# thanks to the norm wilcard - the expresssion matrices are defined in the config
expr = as.data.frame(fread(snakemake@input[["exp"]], sep="\t"))
print("Input expression matrix read into R.")

# set the genes as rownames
rownames(expr) = expr$V1
expr$V1 = NULL

# ensure the sample names are in the correct format (no dcc extension and no periods
if(all(grepl("\\.dcc", colnames(expr)))){
  print("All samples in expr are ending in .dcc. Remove this extension to match the metadata sample IDs")
  colnames(expr) = gsub("\\.dcc", "", colnames(expr))
}

if(!all(grepl("\\.", colnames(expr)))){
  print("WARNING: You have some column names in your expression matrix that contain periods")
}

# now log transform as long as not one of the norms provided
if(!snakemake@wildcards[["norm"]] %in% c("GD_poiss_norm", "GD_poiss_split_norm")){
  print(paste0("Log applied for ", snakemake@wildcards[["norm"]]))
  expr = log2(expr)
}else{
  print(paste0("Log NOT applied for ", snakemake@wildcards[["norm"]]))
}

### Read genes selected by previous rule (on variance)
gene_selection = as.data.frame(fread(snakemake@input[['gene_selection']], header=T))[,1]

nans = apply(expr, 2, function(x) sum(is.na(x)))
nans = names(nans[which(nans>0)])
print("Length of nans vector:")
print(length(nans))             
expr = expr[which(rownames(expr) %in% gene_selection), which(! colnames(expr) %in% nans)]

### Take in account situation where there is only 0 or 1 gene in the filtered dataframe
### (ie: if QC removed almose all genes in user list of genes for clustering) 

if(nrow(expr) < snakemake@params[["maxK"]]){
  file.create(snakemake@output[["done"]])
  failed_df = as.data.frame(matrix(c("Failed, not enough genes to cluster samples")))
  write.table(failed_df, snakemake@output[["updated_meta"]], sep="\t", row.names=T, quote=F)
}else{
  # read in the meta data
  meta = as.data.frame(read_excel(snakemake@input[["meta"]]))
  # make sample IDs rownames of meta and print the top left corner  
  rownames(meta) = meta$Sample_ID
  print(meta[1:5,1:5])
  # find out which samples are no template controls and remove these from meta and expr
  find_NTC = apply(apply(meta, 1, function(x) grepl("emplate", x)), 2, sum)
  find_NTC = names(find_NTC[which(find_NTC==1)])
  print("Printing contents of find_NTC:")
  print(find_NTC)
  meta = meta[which(!rownames(meta) %in% find_NTC),
              intersect(colnames(meta), snakemake@params[["meta_heatmap"]])]
  expr = expr[,which(colnames(expr) %in% rownames(meta))]
  meta = meta[which(rownames(meta) %in% colnames(expr)),]
  
  # reorder expr based on meta rownames
  expr = expr[, rownames(meta)]
  
  output_dir = unlist(strsplit(snakemake@output[["done"]], "\\/"))
  output_dir = paste0(output_dir[1:length(output_dir)-1], collapse="/")
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # execute the run_clustering function
  print("Now execute the run_clustering function")                       
  meta_w_clust = run_clustering(d=expr, an=meta, snakemake@params[["minK"]], snakemake@params[["maxK"]], output_dir)
  
  write.table(meta_w_clust, snakemake@output[["updated_meta"]], sep="\t", row.names=T, quote=F)
  
  file.create(snakemake@output[["done"]])
}