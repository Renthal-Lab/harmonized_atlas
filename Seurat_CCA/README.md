**Step 0: Initialization and Library Loading**

Libraries: The script starts by loading necessary R packages (Seurat, tibble, dplyr, gridExtra, grid, ggplot2, hdf5r, TeachingDemos, and others) which are essential for data manipulation, visualization, and scRNA-seq data analysis.
Parameters Setting: It sets several parameters like res (resolution for clustering), minFeature (minimum number of features), maxpercent.mt (maximum percentage of mitochondrial genes allowed), and others which are used throughout the script.

**Step 1: Data Loading and Merging**
Data Loading: Loads scRNA-seq datasets from different studies stored as .Rds files.
Data Merging: Merges these datasets into a single Seurat object, seurat_mat, for integrated analysis.

**Step 2: Quality Control (QC) and Normalization**
Mitochondrial Gene Percentage Calculation: Calculates the percentage of mitochondrial gene expression as a quality control metric.
Data Subsetting: Filters cells based on the number of detected features and mitochondrial gene content.
Normalization: Performs data normalization and identifies variable features for each dataset independently.

**Step 3: Integration and Dimensionality Reduction**
Feature Selection: Identifies features for integration across different datasets.
Data Integration: Integrates data using the identified features to correct for batch effects.
PCA and Clustering Preparation: Performs PCA and estimates the number of principal components to use for further analysis.

**Step 4: Clustering and Visualization**
UMAP and t-SNE: Generates UMAP and t-SNE plots for visualizing the data in a reduced dimensionality space.
Clustering: Performs clustering on the integrated and scaled data.
Marker Gene Visualization: Generates plots for marker genes using both UMAP and t-SNE projections.

**Step 5: Saving Results and Finding Marker Genes**
Saving Files: Saves various files, including merged Seurat object, QC plots, UMAP, and t-SNE plots.
Marker Gene Analysis: Identifies and exports marker genes for each cluster. Enrichment scores for each clusters are relative to all other clusters in UMAP space

**Additional Notes**
Directory Management: The script creates directories for saving output files and sets the working directory to manage outputs systematically.
Logging: Initiates logging of console output to a file, allowing tracking of the script's progress and debugging if needed.
