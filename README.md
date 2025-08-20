#  Integrated scRNA-seq Analysis of Hepatoblastoma: Tumor, Background, and PDX

This project contains a full analysis workflow for scRNA-seq data using Seurat, including integration, clustering, differential expression (DE) analysis, visualization, and GO enrichment.

---

## Table of Contents
1. [Data Integration & Clustering](#-1-data-integration--clustering)  
2. [Marker Gene Visualization](#-2-marker-gene-visualization)  
3. [Differential Expression (DE) Analysis](#-3-differential-expression-de-analysis)  
4. [Heatmaps of Top DE Genes](#-4-heatmaps-of-top-de-genes)  
5. [GO Enrichment Analysis](#-5-go-enrichment-analysis-biological-processes)  
6. [Outputs Summary](#-outputs-summary)  
7. [Notes](#-notes)  

---
## 1. Scripts
- [Integration & Clustring](scripts/Integrate_ScRNA-Seq)
- [Visualization Plots](scripts/visualization_plots)

## 2. Data Integration & Clustering
- Samples (Patients + Types) were integrated using Seurat’s integration pipeline.  
- Batch effect (patient effect) was regressed out.  
- UMAP visualization of clusters and sample grouping.  

 Example:

![DimPlot_patient_type.png](results/clustering/DimPlot_clusters.png)  

---

## 2. Marker Gene Visualization
- FeaturePlot and ViolinPlot were used for selected marker genes (CYP2B6, PZP, POU5F1, AFP, ALB).  
- DotPlot was generated to compare reference markers across clusters.  

 Examples:

![FeaturePlot markers](results/markers/FeaturePlot_markers.png)  
![ViolinPlot markers](results/markers/ViolinPlot_markers.png)  
![DotPlot markers](results/markers/DotPlot_reference.png)

---

## 3. Differential Expression (DE) Analysis
Main pairwise comparisons:
- Tumor vs Background  
- Tumor vs PDX  
- PDX vs Background  

For each comparison:
- DE was performed using FindMarkers (Wilcoxon test, patient effect controlled).  
- Results saved as CSV files (DE_<comparison>.csv).  

 Example Volcano plots:

![Volcano Tumor vs Background](results/volcano/Volcano_tumor_vs_background.png)  
![Volcano Tumor vs PDX](results/volcano/Volcano_tumor_vs_PDX.png)  
![Volcano PDX vs Background](results/volcano/Volcano_PDX_vs_background.png)

---

## 4. Heatmaps of Top DE Genes
- Selected Top 20 DE genes (based on |log2FC|) for each comparison.  
- Heatmaps generated with sample annotation by Type.  

 Example:

![Heatmap Tumor vs Background](results/heatmaps/Heatmap_tumor_vs_background.png)

---
## 5. GO Enrichment Analysis (Biological Processes)
- Performed separately for Upregulated and Downregulated genes.  
- Enrichment done with clusterProfiler::enrichGO.  

Outputs include:
- DotPlot PDFs (GO_<comp>_Up_dotplot.pdf, GO_<comp>_Down_dotplot.pdf)  
- CSV files with enrichment results (GO_<comp>_Up_enrichment.csv)  

 Example DotPlots (PDF links):  

Tumor vs Background
- [Upregulated](results/enrichment/GO_tumor_vs_background_Up_dotplot.pdf)  
- [Downregulated](results/enrichment/GO_tumor_vs_background_Down_dotplot.pdf)  

Tumor vs PDX
- [Upregulated](results/enrichment/GO_tumor_vs_PDX_Up_dotplot.pdf)  
- [Downregulated](results/enrichment/GO_tumor_vs_PDX_Down_dotplot.pdf)  

PDX vs Background
- [Upregulated](results/enrichment/GO_PDX_vs_background_Up_dotplot.pdf)  
- [Downregulated](results/enrichment/GO_PDX_vs_background_Down_dotplot.pdf)
