**README for GO_analysis_mouse_DRG.R**

**Step 0: Initialization and Library Loading**
-Libraries: The script starts by loading necessary R packages (topGO, Seurat, tidyverse and org.Mm.eg.db) which are essential for analysis.


**Step 1: Data Loading**
-Data Loading: Seurat object and findallmarker result from that object


**Step 2: Define background gene and input gene list**
-background gene: all genes in seurat object
-input gene list: 100 marker genes in each cell type sort by log2FC

**Step 3: GO analysis**
- Genes were annotated for their biological process and associated gene ontology terms. 

