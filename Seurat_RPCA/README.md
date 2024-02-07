**Step 1: Data Loading**
Seurat Object Creation: Reads in scRNA-seq data stored in an RDS file into a Seurat object, setting the default assay to RNA.

**Step 2: Quality Control and Normalization**
Mitochondrial Content Calculation: Calculates the percentage of mitochondrial genes for each cell to assess cell quality.
Quality Control Visualization: Generates violin plots for the number of features, total counts, and mitochondrial content per cell, and saves the plots as "QC.pdf".
Data Subsetting: Filters cells based on feature RNA counts and mitochondrial content.
Data Processing: Normalizes the data, identifies variable features, and scales the data. PCA is performed to reduce dimensionality.

**Step 3: PCA, Clustering, and Dimension Reduction**
Determining Number of PCs: Calculates the cumulative percentage of variance explained to select the number of principal components.
Harmony Integration: Applies the Harmony algorithm for batch correction using PCA features, specifying datasets as batches.
Dimensionality Reduction: Performs UMAP and t-SNE for visualization, constructs a neighbor graph based on Harmony dimensions, and identifies clusters.

**Step 4: Generating Results**
Visualization: Generates UMAP and t-SNE plots showing the integrated data colored by dataset, and saves these visualizations.
Marker Gene Visualization: Visualizes expression patterns of predefined marker genes for neuron and non-neuron cell types using UMAP plots.

**Step 5: Saving Files and Finding Marker Genes**
Saving Processed Data: Saves the integrated Seurat object as an RDS file for future analysis.
Marker Gene Identification: Identifies marker genes for each cluster and saves the results.
