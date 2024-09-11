# env : visu.yaml

log = file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library(dplyr)
library(ggforce)
library(ggplot2)
library(data.table)
library(reshape2)
library(cowplot)
library(patchwork)
library(ggpubr)
library(umap)
library(readxl)
library(tibble)

bigColSet = c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(10, "Set3"))

col_variables = snakemake@params[["col_variables"]]
meta = as.data.frame(read_excel(snakemake@input[["meta"]]), stringsAsFactors=FALSE)

col_variables = col_variables[which(col_variables %in% colnames(meta))]

bigColSet = c(RColorBrewer::brewer.pal(9, "Set1")[1:5],
              RColorBrewer::brewer.pal(8, "Dark2"),
              RColorBrewer::brewer.pal(8, "Set2"))

if("Sample_ID" %in% colnames(meta)){
  rownames(meta) = meta$Sample_ID
}else{
  samp_id = colnames(meta)[grepl("^DSP", as.vector(meta[1,]))]
  rownames(meta) = meta[,samp_id]
}
if(!any(grepl("\\.dcc", rownames(meta)))){
  rownames(meta) = paste0(rownames(meta), ".dcc")
}

find_NTC = apply(apply(meta, 1, function(x) grepl("emplate", x)), 2, sum)
find_NTC = names(find_NTC[which(find_NTC==1)])
meta = meta[which(!rownames(meta) == find_NTC),]

if(snakemake@params[["create_ROI_type_variable"]]){
  meta = meta %>%
    rownames_to_column('ID') %>%
    group_by(`slide name`, roi) %>%
    mutate(ROI_content = paste(sort(gsub("-aoi-001", "", aoi)), collapse=" | ")) %>%
    as.data.frame() %>%
    column_to_rownames('ID')
  
  col_variables = c(col_variables, "ROI_content")
}
print(table(meta$ROI_content))
print("#################################################@")

sum_0 = function(x){
  return(sum(x==0 | is.na(x)) / length(x))
}

# select_exp_features = function(exp, variance_select, prevalence_filter){
#     ### remove genes expressed in less than <prevalence filter> % of samples
#     if(!is.null(prevalence_filter)){
#         prev_vect = as.numeric(apply(exp, 1, sum_0))
#         exp = exp[which(prev_vect < prevalence_filter),]
#     }
#     ### Select genes based on variance
#     if(!is.null(variance_select)){
#         perGene_var = apply(exp, 1, var, na.rm=T)
#         perGene_var = sort(perGene_var, index.return=TRUE, decreasing=TRUE)
#         perGene_var = perGene_var$ix[1:variance_select]
#         exp = exp[perGene_var,]
#     }
#     return(exp)
# }

UMAP_exp = function(expr, meta, output_file_plot, output_file_df, logged=F,
                    annots = c("Section", "segment", "slide name"),
                    variance_select = 2000, prevalence_filter = 0.1,
                    title=""){
  
  custom_umap <- umap::umap.defaults
  custom_umap$random_state <- 42
  annot = meta[,c(annots, "roi")]
  expr = expr[,colSums(is.na(expr))<nrow(expr)]
  #expr = select_exp_features(expr, variance_select, prevalence_filter)
  dat <- t(expr)
  
  if(logged){
    if(min(dat, na.rm=T)<=0){
      return()
    }
    dat = log2(dat) 
    title = paste0(title, "(log expr.)")
  }
  dat_umap = umap(dat, config=custom_umap)
  dat <- as.data.frame(dat)
  dat[, c("UMAP1", "UMAP2")] <- dat_umap$layout[, c(1,2)]
  
  dat = merge(dat, annot, by=0, all.x=TRUE)
  write.table(dat[,c(c(colnames(dat)[1], "UMAP1", "UMAP2"), annots)], file=output_file_df, sep="\t", row.names=F, quote=F)
  
  
  list_plots = list()
  for(ann in annots){
    print(ann)
    
    p <- ggplot(data = dat, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(fill = !!sym(ann)), shape=21, size=2.5, colour="white") +
      theme_bw()+
      theme(legend.position="bottom") + guides(fill=guide_legend(ncol=2,byrow=TRUE))+
      labs(title = title)
    
    if(! is.numeric(dat[,ann])){
      p = p + scale_color_manual(values = bigColSet[1:length(unique(dat[,ann]))] )
    }
    list_plots <- c(list_plots,list(p))
  }
  if(length(annots)>=2){
    width_ = 750
    height = ceiling(length(annots)/2) * 400
  }else{
    width_ = 375
    height = 400
  }
  
  png(file=output_file_plot, width=width_*4, height=height*4, res = 200)
  print(patchwork::wrap_plots(list_plots, ncol=2))
  dev.off()
}

plot_gene = function(expr, plt_dat, genes, col_x, col_fill){
  genes_df = as.data.frame(t(expr))
  if(length(setdiff(genes, colnames(genes_df)))>0){
    print("Warning, missing genes:")
    print(setdiff(genes, colnames(genes_df)))
    genes = intersect(colnames(genes_df), genes)
  }
  
  if(length(genes) == 1){genes=genes[1]}
  
  if(length(genes) == 0){return(ggplot())}
  genes_df = genes_df[,genes, drop=F]
  
  plt_dat = merge(plt_dat,genes_df,by="row.names",all.x=TRUE)
  plt_dat = reshape2::melt(plt_dat, measure.vars = genes)
  
  colnames(plt_dat)[which(colnames(plt_dat) == "variable")] = "Gene"
  colnames(plt_dat)[which(colnames(plt_dat) == "value")] = "Norm_exp"
  
  if(any(is.na(plt_dat[,col_x]))){
    print("Warning, some line of the metadata have NA values in segment column:")
    print(plt_dat[which(is.na(plt_dat[,col_x])),])
    plt_dat = plt_dat[which(!is.na(plt_dat[,col_x])),]
  }
  
  number_of_rows = ceiling(length(genes) / 10) ## max 10 panels per line
  
  plt = ggplot(data = plt_dat, aes(x=!!sym(col_x), y=Norm_exp)) + 
    facet_wrap(.~Gene, nrow=number_of_rows, scales="free_y") + theme_bw() + 
    geom_boxplot() +
    # geom_violin() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) + 
    geom_jitter(width=.07, shape=21, size=1.5, colour="white", aes(fill=!!sym(col_fill)))
  
  if(! is.numeric(plt_dat[,col_fill])){
    plt = plt + scale_fill_manual(values = bigColSet[1:length(unique(plt_dat[,col_fill]) )] )
  }
  
  return(plt)
}


# now read in the expression matrix and make genes rownames
expr = as.data.frame(fread(snakemake@input[["exp"]]), header=T)
rownames(expr) = expr$V1
expr$V1 = NULL

# run a check on the expression matrix to ensure the sample names
# do not have periods. Should have dashes instead
if(!all(grepl("\\.", colnames(expr)))){
  print("WARNING: You have some column names in your expression matrix that contain periods")
}

##### Filter using gene_selection file
gene_selection = as.data.frame(fread(snakemake@input[['gene_selection']], header=T))[,1]


##### Plot cell population specific markers
pdf(file=snakemake@output[['marker']], width=15, height=15)
for(col_fill in cols){
  a = plot_gene(expr, meta, snakemake@params[["immune_mark"]], col_x, col_fill) + ggtitle("Immune markers")
  b = plot_gene(expr, meta, snakemake@params[["epith_mark"]], col_x, col_fill) + ggtitle("Epithelial and cancer markers")
  c = plot_gene(expr, meta, snakemake@params[["stroma_mark"]], col_x, col_fill) + ggtitle("Stromal markers")
  print(a/b/c)
}
dev.off()