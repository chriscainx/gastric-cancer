# Processed Data in R Compatible Formats

This folder contains the processed single cell data compatible with R analysis, to facilitate the demonstrative usage in this repository.

- `gcdata` is a 10X formatted sparse matrix bundle containing cell expression matrix in UMI counts, cell IDs and gene IDs. This folder can be read by the Seurat package.
- `cell_metadata.csv` is the cell metadata file containing sample IDs and cell type labels.
- `bulk_rna_tpm.txt` is the bulk RNA-seq data in TPM unit.

An R Seurat object can be generated from the gene expression counts matrix and metadata using following commands:

```r
gcmatrix <- Read10X("gcmatrix")
gcmeta <- read_csv("cell_metadata.csv")
gcdata <- CreateSeuratObject(gcmatrix, meta.data = gcmeta)
gcdata <- NormalizeData(gcdata, verbose = FALSE)
gcdata <- ScaleData(gcdata, verbose = FALSE)
```

Reduced dimensions such as PCA and UMAP can be generated from the normalized expression counts using following commands:

```r
gcdata <- FindVariableFeatures(gcdata, verbose = FALSE)
gcdata <- RunPCA(gcdata, npcs = 30, verbose = FALSE)
gcdata <- RunUMAP(gcdata, reduction = "pca", dims = 1:30, verbose = FALSE)
```
