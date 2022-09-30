# Code for Data Processing and Clustering

This folder contains the code for data processing and clustering on the original dataset.

- `01.DataNormalization.R` is the data normalization step conducted using the Seurat package, after alignment with the 10X CellRanger software.
- `02.ICA.py` is the data adjustment step to detect and remove manifestation of the HSP gene module activation in the single cell transcriptomes caused by enzyme stimulation during single cell isolation.
- `03.Clustering.py` contains the steps and parameters used for identification of stable cell clusters across sequencing chemistries.

Due to rapid updating of the softwares and packages used in the single cell RNA-seq analysis field, the behaviors of APIs or function calls contained in the code may change in newly assembled analysis environments. However, our conclusions will be stable when analyzed in alternative environments given same mathematical and statistical framework is employed.