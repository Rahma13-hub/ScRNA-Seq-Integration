# Integrated scRNA-seq Analysis of Hepatoblastoma: Tumor, Background, and PDX
This project contains a full analysis workflow for scRNA-seq data using Seurat, including integration, clustering, differential expression (DE) analysis, visualization, and GO enrichment.

## Table of Contents
1. [Data Integration & Clustering](#-1-data-integration--clustering)  
2. [Marker Gene Visualization](#-2-marker-gene-visualization)  
3. [Differential Expression (DE) Analysis](#-3-differential-expression-de-analysis)  
4. [Heatmaps of Top DE Genes](#-4-heatmaps-of-top-de-genes)  
5. [GO Enrichment Analysis](#-5-go-enrichment-analysis-biological-processes)  
6. [Outputs Summary](#-outputs-summary)

## 1. Data Integration & Clustering
- Samples (Patients + Types) were integrated using Seurat’s integration pipeline.  
- Batch effect (patient effect) was regressed out.  
- UMAP visualization of clusters and sample grouping.  

 Example:

![DimPlot by cluster](results/DimPlot_clusters.png)  
![DimPlot by patient](results/DimPlot_patient.png)

## 2. Marker Gene Visualization
- FeaturePlot and ViolinPlot were used for selected marker genes (CYP2B6, PZP, POU5F1, AFP, ALB).  
- DotPlot was generated to compare reference markers across clusters.

Examples:

![FeaturePlot markers](results/FeaturePlot_markers.png)  
![ViolinPlot markers](results/ViolinPlot_markers.png)  
![DotPlot markers](results/DotPlot_reference.png)

## 3. Differential Expression (DE) Analysis
Main pairwise comparisons:
- Tumor vs Background  
- Tumor vs PDX  
- PDX vs Background  

For each comparison:
- DE was performed using FindMarkers (Wilcoxon test, patient effect controlled).  
- Results saved as CSV files (DE_<comparison>.csv).  

 Example Volcano plots:

![Volcano Tumor vs Background](results/Volcano_tumor_vs_background.png)  
![Volcano Tumor vs PDX](results/Volcano_tumor_vs_PDX.png)  
![Volcano PDX vs Background](results/Volcano_PDX_vs_background.png)

---

## 4. Heatmaps of Top DE Genes
- Selected Top 20 DE genes (based on |log2FC|) for each comparison.  
- Heatmaps generated with sample annotation by Type.  

 Example:

![Heatmap Tumor vs Background](results/Heatmap_tumor_vs_background.png)

---

## 5. GO Enrichment Analysis (Biological Processes)
- Performed separately for Upregulated and Downregulated genes.  
- Enrichment done with clusterProfiler::enrichGO.  

Outputs include:
- DotPlot PDFs (GO_<comp>_Up_dotplot.pdf, GO_<comp>_Down_dotplot.pdf)  
- CSV files with enrichment results (GO_<comp>_Up_enrichment.csv)  

 Example DotPlots:

![GO Up Tumor vs Background](results/GO_tumor_vs_background_Up.png)  
![GO Down Tumor vs Background](results/GO_tumor_vs_background_Down.png)

---

## Outputs Summary

Figures:
- DimPlots, FeaturePlots, ViolinPlots, DotPlots  
- Volcano plots  
- Heatmaps  
- GO enrichment dotplots  

Tables:
- DE results per comparison (DE_*.csv)  
- GO enrichment results (GO_*_enrichment.csv)  
