# This file contains steps for data normalization from unfiltered cell-gene 
# expression matrices using Seurat and SCTransform

# Load libraries
library(dplyr)
library(readr)
library(BiocParallel)
library(Seurat)
library(sctransform)

# Load sample metadata table
meta <- read_csv("sample_metadata.csv")

# Initialize seurat objects (filter cells, keep genes)
raw.data <- bplapply(1:nrow(meta), function(i){
  data <- Read10X_h5(meta$File[i], use.names = FALSE)
  colnames(data) <- paste0(meta$Sample[i], "-", colnames(data))
  data
}, BPPARAM = MulticoreParam(32)) %>% do.call(cbind, .)

gcall <- CreateSeuratObject(raw.data, min.cells = 3, min.features = 200)
gcall$Sample <- substr(colnames(gcall), 1, 7)
rownames(meta) <- meta$Sample

gcall$Patient <- meta[gcall$Sample,]$Patient
gcall$Tissue <- meta[gcall$Sample,]$Tissue
gcall$Platform <- meta[gcall$Sample,]$Platform

# Filter empty and low quality droplets
molecule.info <- read_tsv("data/features.tsv.gz", col_names = FALSE)
mito.features <- molecule.info$X1[grep(pattern = "^MT-", x = molecule.info$X2)]
percent.mito <- Matrix::colSums(x = GetAssayData(object = gcall, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = gcall, slot = 'counts'))
gcall$percent.mito <- percent.mito
gcall <- subset(x = gcall, subset = nCount_RNA < 75000 & nFeature_RNA < 7500 & percent.mito < 0.1)

# Normalize data using platform as a blocking var, get corrected umi counts and pearson residual from SCTransform for
# downstream ICA.
# Note Seurat and SCTransform behaviors may have changed due to new versions released, please
# proceed according to the latest guidelines from Seurat and SCTransform.
options(mc.cores = 12)
set.seed(64)
gcall$log_umi <- log10(gcall$nCount_RNA)
sct_result <- sct(GetAssayData(object = gcall, slot = 'counts'), return_corrected_umi = TRUE, batch_var = 'Platform', 
  n_genes = 3000, latent_var_nonreg="percent.mito", cell_attr=gcall@metadata)

# The data normalization using SCTransform ends here.
# Save data to h5ad format, for next step processing with python

library(reticulate)
anndata = import("anndata",convert=FALSE)
bbknn = import("bbknn", convert=FALSE)
sc = import("scanpy.api",convert=FALSE)
np = import("numpy")

np$save("residuals.npy", sct_result$y)
adata = anndata$AnnData(X=sct_result$umi_corrected, obs=gcall@metadata)
adata$write("gcall.correctedumi.h5ad"))