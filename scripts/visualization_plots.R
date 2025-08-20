# Set cluster identities for downstream analysis
Idents(seurat.integrated) = "seurat_clusters"
# Define marker gene sets for hepatocyte, HB, and stemness signatures
hep_genes  = c("ALB","APOA1","TTR","CYP2B6")   
hb_genes   = c("AFP","DLK1","GPC3")            
stem_genes = c("POU5F1","NANOG","SOX2")        
# Add module scores (per-cell enrichment) for the defined gene sets
seurat.integrated = AddModuleScore(
  seurat.integrated,
  features = list(hep_genes, hb_genes, stem_genes),
  name = "RefScore"
)
view(seurat.integrated@meta.data)
# List cluster IDs
clusters = levels(Idents(seurat.integrated))
# Define marker genes of interest
genes_use = c("CYP2B6","PZP","POU5F1","AFP","ALB")
# FeaturePlot per cluster 
# Loop through clusters and generate FeaturePlots for each selected gene
for (cl in clusters) {
  cells_use = WhichCells(seurat.integrated, idents = cl)
  # store plots for current cluster
  plots = list()
  for (gene in genes_use) {
    p = FeaturePlot(
      seurat.integrated,
      cells = cells_use,
      features = gene
    ) + ggtitle(gene)  
    
    plots[[gene]] = p
  }
  
  # Combine plots into a grid for the current cluster
  final_plot = wrap_plots(plots, ncol = 3) + 
    plot_annotation(title = paste("Cluster", cl))
  
  print(final_plot) 
}
#  Violin plots per cluster 
# Generate violin plots of expression for selected genes in each cluster

for (cl in clusters) {
  
  cluster_obj = subset(seurat.integrated, idents = cl)
  plots = list()  
  genes_use = c("CYP2B6","PZP","POU5F1","AFP","ALB")
  
  for (gene in genes_use) {
    p = VlnPlot(
      cluster_obj,
      features = gene,
      pt.size = 0
    ) + ggtitle(gene)           
    plots[[gene]] = p
  }
  # Combine violin plots for each cluster
  final_plot <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = paste("Cluster", cl),  
      theme = theme(plot.title = element_text(hjust = 0, face = "bold", size = 14)) 
    )
  print(final_plot)
}
# DotPlot for reference gene sets
DotPlot(seurat.integrated, features = c(hep_genes, hb_genes, stem_genes)) +
  RotatedAxis()+
  ggtitle("Reference markers across clusters")
head(seurat.integrated@meta.data[, c("Type","patient")])

# Set identities to "Type" (tumor, background, PDX)
Idents(seurat.integrated) = "Type"
# Differential Expression Analysis 
# Define pairwise comparisons of interest
comparisons = list(
  c("tumor", "background"),
  c("tumor", "PDX"),
  c("PDX", "background")
)

# Loop through comparisons and run DE analysis
for(comp in comparisons){
  message("Running DE: ", comp[1], " vs ", comp[2])
  
  de_res = FindMarkers(
    object = seurat.integrated,
    ident.1 = comp[1],
    ident.2 = comp[2],
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    latent.vars = "patient" 
  )
  
  # Save DE results
  write.csv(de_res, paste0("DE_", comp[1], "_vs_", comp[2], ".csv"))
}
# Volcano Plots for DE results 

for(comp in comparisons){
  # Load DE results from file
  de_file = paste0("DE_", comp[1], "_vs_", comp[2], ".csv")
  de_res = read.csv(de_file, row.names = 1)
  # Define significance based on adjusted p-value and fold change
  de_res$significant <- with(de_res, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25, "yes", "no"))
  
  p = ggplot(de_res, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
    ggtitle(paste(comp[1], "vs", comp[2])) +
    xlab("log2 Fold Change") +
    ylab("-log10 adjusted p-value") +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(p)
}
# Heatmaps for top DE genes per comparison 

heatmap_list = list()

for(comp in comparisons){
  # Load DE results for the comparison
  de_file = paste0("DE_", comp, ".csv")
  de_res = read.csv(de_file, row.names = 1)
  # Filter significant genes (adjusted p-value + fold change cutoff)
  sig_res = de_res %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  # Select top 20 most differentially expressed genes by |log2FC|
  top_genes = sig_res %>% arrange(desc(abs(avg_log2FC))) %>% head(20) %>% rownames()
  
  if(length(top_genes) > 1){
    expr_mat = GetAssayData(seurat.integrated, assay = "RNA", slot = "data")[top_genes, ]
    # Add annotation   
    annotation_col = data.frame(Type = seurat.integrated$Type)
    rownames(annotation_col) = colnames(expr_mat)
    # Draw heatmap of top DE genes    
    hm = pheatmap(as.matrix(expr_mat),
                  scale = "row",
                  annotation_col = annotation_col,
                  show_rownames = TRUE,
                  show_colnames = FALSE,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",
                  main = comp)
    
    heatmap_list[[comp]] = hm
  }
}
#######################################
# GO Enrichment Analysis for Up/Downregulated genes

comparisons = c("tumor_vs_background", "tumor_vs_PDX", "PDX_vs_background")

for(comp in comparisons){
  # Load DE results
  
  de_file = paste0("DE_", comp, ".csv")
  de_res = read.csv(de_file, row.names = 1)
  
  # Upregulated genes
  
  up_res = subset(de_res, p_val_adj < 0.05 & avg_log2FC > 0.25)
  up_genes = rownames(up_res)
  up_ids = bitr(up_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  if(nrow(up_ids) > 0){
    ego_up = enrichGO(gene = up_ids$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)
    
    # Save dotplot of enriched terms
    pdf(paste0("GO_", comp, "_Up_dotplot.pdf"))
    print(dotplot(ego_up, showCategory=15) + ggtitle(paste("GO Enrichment (Up):", comp)))
    dev.off()
    # Save enrichment results as CSV
    write.csv(as.data.frame(ego_up), paste0("GO_", comp, "_Up_enrichment.csv"), row.names=FALSE)
  }
  
  # Downregulated genes
  down_res = subset(de_res, p_val_adj < 0.05 & avg_log2FC < -0.25)
  down_genes = rownames(down_res)
  down_ids = bitr(down_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  
  if(nrow(down_ids) > 0){
    ego_down = enrichGO(gene = down_ids$ENTREZID,
                        OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID",
                        ont = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05)
    
    # Save dotplot of enriched terms
    pdf(paste0("GO_", comp, "_Down_dotplot.pdf"))
    print(dotplot(ego_down, showCategory=15) + ggtitle(paste("GO Enrichment (Down):", comp)))
    dev.off()
    # Save enrichment results as CSV
    write.csv(as.data.frame(ego_down), paste0("GO_", comp, "_Down_enrichment.csv"), row.names=FALSE)
  }
}



