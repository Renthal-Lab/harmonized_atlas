import os
import tempfile
import warnings
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
#import scrublet as scr
import scvi
import torch

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)
sc.set_figure_params(figsize=(4, 4))
torch.set_float32_matmul_precision("high")
save_dir = "/mnt/data/projects/Shams_working/"


tg_data = sc.read("Seurat_annotated.h5ad")

print(tg_data)
adata = tg_data
adata.layers["counts"] = adata.X.copy()

print("# cells, # genes ", adata.shape)
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=3)
print("# cells, # genes after filtering:", adata.shape)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.raw = adata  # keep full dimension safe


query_mask = np.array(
        [s in ["Nguyen 2019_10x_Mouse", "Nguyen_2019_Drop-seq_Mouse","Pei 2023_10x_Mouse","Yang 2022_10x_Human","Yang 2022_inDrops_Mouse","Yang_2022_10x_Mouse"] for s in adata.obs["dataset"]]
)

adata_ref = adata[~query_mask].copy()
adata_query = adata[query_mask].copy()


sc.pp.highly_variable_genes(
    adata_ref,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="dataset",
    subset=True,
)

adata_query = adata_query[:, adata_ref.var_names].copy()
scvi.model.SCVI.setup_anndata(adata_ref, layer="counts", batch_key="dataset")
scvi_ref = scvi.model.SCVI(
    adata_ref,
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)
scvi_ref.train()

SCVI_LATENT_KEY = "X_scVI"

adata_ref.obsm[SCVI_LATENT_KEY] = scvi_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_ref)
sc.pl.umap(adata_ref, color=['dataset'],save='ref_level1_TG.pdf')
scvi_ref.save("/mnt/data/projects/Shams_working/TG_level1_ref", overwrite=True)
scvi.model.SCVI.prepare_query_anndata(adata_query, scvi_ref)


scvi_query = scvi.model.SCVI.load_query_data(
    adata_query,
    scvi_ref,
)
scvi_query.train(max_epochs=200, plan_kwargs={"weight_decay": 0.0})
adata_query.obsm[SCVI_LATENT_KEY] = scvi_query.get_latent_representation()
sc.pp.neighbors(adata_query, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)

sc.pl.umap(
    adata_query,
    color=["Species", "dataset"],
    frameon=False,
    ncols=1,
    save="TG_level1_query.pdf"
)
adata_full = anndata.concat([adata_query, adata_ref])
adata_full

adata_full.obsm[SCVI_LATENT_KEY] = scvi_query.get_latent_representation(
    adata_full
)

sc.pp.neighbors(adata_full, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata_full)
sc.tl.umap(adata_full)

sc.pl.umap(
    adata_full,
    color=["Species", "dataset"],
    frameon=False,
    ncols=1,
    save="TG_integrated_level1.pdf"
)

sc.write('adata_tg_level1.h5ad', adata_full)
marker_genes = {
            "General": ["Snap25","Rbfox3","Sparc","Mpz"],
            "Th": ["Tafa4","Th","Gfra2","Pou4f2"],
            "Calca": ["Tac1","Calca","Nefh","Sstr2","Adra2a","Oprk1","Bmpr1b","Smr2","Chrna7"],
            "Nefh": ["Nefh","Pvalb","Calb1","Cadps2","Ntrk3","Ntrk2","S100a16"],
            "Mrgpr": ["Scn11a","Lpar3","Mrpgrd","Mrgpra3","Mrgprb4","Trpv1"],
            "Trpm8": ["Trpm8","Hpca","Trpv1","Rxfp1"],
            "Sst": ["Sst","Nppb","Il31ra"],
            "Injured": ["Atf3","Sox11","Jun"],
            "Satglia": ["Apoe","Fabp7","Ednrb"],
            "Schwann": ["Mpz","Mbp","Scn7a"],
            "Immune": ["Lyz2","Mrc1","Csf1r","Ccr2","Cd3e","Cd8a","Ccr7","Cd4","Cd79a","Cd79b","Ighd","Ighm","Ms4a1","S100a8","S100a9","Retnlg"],
            "Endothelial": ["Pecam1","Cldn5","Egfl7"],
            "Pericyte": ["Pdgfrb","Notch3","Kcnj8"],
            "Fibroblast": ["Dcn","Pdgfra","Mgp","Ptgds","Fxyd5"]
}

celltypes=["General","Th","Calca","Nefh","Trpm8","Sst","Injured","Satglia","Schwann","Immune","Endothelial","Pericyte","Fibroblast"]
#marker genes in our data
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
                marker_genes_in_data[ct] = markers_found

for ct in B_plasma_cts:
        print(f"{ct.upper()}:")  # print cell subtype name
        sc.pl.umap(
                adata,
                color=marker_genes_in_data[ct],
                min=0,
                max="p99",  # set vmax to the 99th percentile of the gene count instead of the maximum, to prevent outliers from making expression in other cells invisible. Note that this can cause problems for extremely lowly expressed genes.
                sort_order=False,  # do not plot highest expression on top, to not get a biased view of the mean expression among cells
                frameon=False,
                cmap="Reds",  # or choose another color map e.g. from here: https://matplotlib.org/stable/tutorials/colors/colormaps.html
                save="markers.pdf"
                )
        print("\n\n\n")  # print white space for legibility
exit()



SCVI_LATENT_KEY = "X_scVI"
adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
sc.tl.leiden(adata)
SCVI_MDE_KEY = "X_scVI_MDE"
adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY])


sc.pl.embedding(
    adata,
    basis=SCVI_MDE_KEY,
    color=["batch", "leiden"],
    frameon=False,
    ncols=1,
)

sc.write('adata_tg_level1.h5ad', adata)
                                                 
