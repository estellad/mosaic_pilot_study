# install.packages("enrichR")
library(enrichR)

contrast <- "1&5&9_14"
save_path_DE <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures/Fig_DE_Volcano/B3_2"
DE_result <- read.csv(file.path(save_path_DE, paste0(contrast, ".csv")))

results <- DE_result4
results <- as.data.frame(results)
# results <- results %>% filter(cluster == "159_TME")
results <- results %>% filter(cluster == "159")
# results <- results %>% filter(cluster == "14_TME")
# results <- results %>% filter(cluster == "14")
gene_list <- results$gene # [1:100]

# Perform enrichment analysis using Enrichr with the specified gene list and database
# databases <- c("KEGG_2021_Human", "Reactome_2021", "GO_Biological_Process_2021", 
#                "GO_Molecular_Function_2021", "GO_Cellular_Component_2021", 
#                "MSigDB_Hallmark_2020", "WikiPathways_2021_Human", "Panther_2021")
databases <- c("KEGG_2021_Human", "MSigDB_Hallmark_2020", "MSigDB_Computational", "MSigDB_Oncogenic_Signatures")

enrichment_results <- enrichr(gene_list, databases)

# enrichment_results4 <- enrichment_results$KEGG_2021_Human
# enrichment_results <- enrichment_results$GO_Biological_Process_2021
# enrichment_results <- enrichment_results$WikiPathways_2021_Human
enrichment_results <- enrichment_results$MSigDB_Hallmark_2020 
# enrichment_results2 <- enrichment_results$MSigDB_Computational 
# enrichment_results3 <- enrichment_results$MSigDB_Oncogenic_Signatures 

head(enrichment_results)

enrichment_results_plt_df <- data.frame(enrichment_results) %>%
  arrange(Odds.Ratio) %>%
  tail(8)

head(enrichment_results_plt_df)

enrichment_results_plt_df$Term <- factor(enrichment_results_plt_df$Term,
                                         levels = enrichment_results_plt_df$Term)

# Plot barplot of adjusted p-values
p <- ggplot(data = enrichment_results_plt_df, aes(x=Term, y=Odds.Ratio, fill=-log10(P.value))) +
  geom_bar(stat="identity") + 
  ggtitle("Enrichment Analysis Results") + 
  labs(fill = expression("-log"[10]*"p-value")) + 
  xlab("Enriched Terms") + 
  ylab("Odds Ratio") + coord_flip() + 
  scale_fill_gradient(low =  "#fc98fc", # "#9393ff",  
                      high = "#ff00ff" # "#0000ff" 
                      ) + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), 
        plot.title = element_blank(), 
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.position = "bottom"
        ) # + ylim(0,30)

plot_title = "Vis_B3_2_DE_Pathway_159.pdf"
# plot_title = "Vis_B3_2_DE_Pathway_14.pdf"
pdf(file = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Poster_Figures/", plot_title),
    width = 8, # 7
    height = 4)
print(p)
dev.off()

# TODO: dotplot for tumor subtype cluster
# TODO: from volcano plot, and check drug targets










