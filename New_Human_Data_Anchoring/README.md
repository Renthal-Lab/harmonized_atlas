# README

## Inputs

'atlas': Object created by subsetting neuronal atlas Seurat object for Primate species only. Command used was:
```R
atlas <- subset(neuron_seurat, subset = Species %in% c("C. Macaque","R. Macaque"))
```
obj: Seurat object containing counts from raw data given by other Pain Consortium Centres. The genes for humans have been converted to their respective mice orthologs as per code given elsewhere in this repository. This is the object to be anchored to the atlas object.

## Steps

1. Setup environment (Lines 3-4)
2. Load input data (Lines 6-7)
3. Find the top 2000 variable features and the top 50 principal components from the raw counts of the object to be anchored to the atlas (Lines 9-12)
4. Process and generate PCs from the raw RNA counts (Lines 14-18)
5. Anchor, transfer labels and project the cells from the human data to the primate neuronal atlas using RNA-to-RNA anchoring but projecting the cells in the UMAP space generated from the integrated assay of the atlas object (Lines 20-25)
6. Save object (Line 27)  
