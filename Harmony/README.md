**Data Preparation**
- Data Loading: Reads a Seurat object containing scRNA-seq data from an RDS file.
- Quality Control: Calculates the percentage of mitochondrial gene expression (percent.mt) for each cell as an initial quality control measure. It also generates violin plots for the number of detected features, total counts, and mitochondrial content per cell, saving these plots as "QC.pdf".
- Data Filtering: Subsets the data to retain cells with more than 200 detected features and less than 10% mitochondrial content.
  
**Data Normalization and Integration**
- Normalization and Feature Selection: Normalizes gene expression data, identifies variable features, scales the data, and performs principal component analysis (PCA) as preparatory steps for data integration.
- Harmony Integration: Integrates the scRNA-seq datasets using Harmony to correct for batch effects, specifying datasets as batches. PCA dimensions are used for integration, and the number of dimensions to use is determined based on the explained variance.
- Dimensionality Reduction and Clustering: After integration, the script runs UMAP for dimensionality reduction, constructs a neighbor graph based on the Harmony dimensions, and identifies clusters using an unsupervised clustering algorithm.

**Visualization and Saving Results**
- UMAP Visualization: Generates and saves a UMAP plot that visualizes the integrated data, colored by dataset, to assess the effectiveness of the batch correction and the overall integration.
- Saving the Integrated Object: Saves the processed and integrated Seurat object as an RDS file for future analysis or sharing.

**Marker Gene Analysis**
- Marker Gene Lists: Defines lists of marker genes for neurons and non-neuron cell types to facilitate detailed cellular composition analysis within the integrated dataset.
- Feature Plot Generation: Iterates over the marker gene lists to generate and save UMAP plots for each set of genes, highlighting their expression patterns across the dataset.

