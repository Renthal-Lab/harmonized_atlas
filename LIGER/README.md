**Data Loading and Conversion**
Data Loading: Loads a Seurat object from a specified path, which contains scRNA-seq data ready for integration.
Conversion: Converts the Seurat object into a LIGER object for integrative analysis, specifying the use of RNA assay and including metadata variables for the dataset.

**Data Processing**
- Normalization: Normalizes the LIGER object data.
- Variable Features Selection: Selects genes for integration based on variance, specifying to use 2000 genes
- Scaling: Scales the data without centering.
- Optimization: Performs optimization using Alternating Least Squares (ALS) with the specified number of dimensions.
- Quantile Normalization: Applies quantile normalization across datasets to align the distribution of gene expression values.
- Clustering: Clusters cells using the Louvain algorithm with a specified resolution.
- Dimensionality Reduction and Visualization: Generates a UMAP representation of the integrated data and saves the plots to a PDF file named liger_umap.pdf.
- Gene Expression Visualization
- Plots expression levels of selected marker genes on UMAP coordinates and saves the plots to PDF files. This section demonstrates how to visualize gene expression across the integrated dataset for specific genes of interest.

**Differential Expression Analysis**
Performs differential expression analysis comparing either clusters or datasets and saves the results to CSV files.

**Saving Results**
- LIGER Object: Saves the processed LIGER object for future use.
- Conversion to Seurat: Converts the LIGER object back to a Seurat object and saves it, allowing for further analysis or visualization using Seurat's functions.
