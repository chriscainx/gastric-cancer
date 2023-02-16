# Supporting materials for gastric cancer project

Here, we deposit the supporting code and data of our gastric cancer research (Kang et al., 2022) for availability to the broad scientific community. The code and data deposited here are free for academic use. Any other types of use are prohibited.

Please note that GitLFS charges a fee for data downloading bandwidth. If you encounter the `This repository is over its data quota` notification, it indicates the free bandwidth has been used up by other downloaders. Please consider downloading the gene expression data directly from GEO. Code deposited in this repository will not be affected by the quota.

## Introduction

### Background
The tumor microenvironment (TME) has been shown to strongly influence treatment outcome for cancer patients in various indications and to influence overall survival. For gastric cancer (GC) the cells forming the TME have not been extensively characterized.

### Results
Here, we combined bulk and single-cell RNA-sequencing from tumors and matched normal tissue of 24 treatment-naïve GC patients to better understand which cell types and transcriptional programs are associated with malignant transformation of the stomach. Clustering 96,623 cell of non-epithelial origin revealed 81 well-defined TME cell types. Activated fibroblasts and endothelial cells were most prominently overrepresented in tumors. Intercellular network reconstruction as well as survival analysis of an independent cohort implied the importance of these cell types together with immunosuppressive myeloid cell subsets and regulatory T cells in establishing an immunosuppressive microenvironment correlating with worsened prognosis and lack of response in anti-PD1 treated patients. In contrast, a subset of IFNγ activated T cells and HLA-II expressing macrophages were found to be linked to response and increased overall survival.

### Conclusions
Our gastric cancer single cell TME compendium together with the matched bulk transcriptome data provides a unique resource for the identification of new potential biomarkers for patien stratification and helps further elucidate gastric cancer biology and insights for therapy.

## Contents of this repository

We organized this repository to host important code used in the research. Processed data and demo for explorative use (such as integrating with other tumor microenvironment data) are also deposited to facilitate further scientific investigations. In this repository, the following folders can be found.

```
project
├───00.Processed Data
├───01.Data Processing and Clustering
├───02.Cell Communication
└───03.Integration with Other Atlas
```

- `00.Processed Data` stores cell-gene normalized expression data and cell metadata.
- `01.Data Processing and Clustering` stores the code for data processing and cell clustering.
- `02.Cell Communication` stores the code for cell communication analysis, and supplementary information.
- `03.Integration with Other Atlas` stores an example of integrating this dataset with other tumor microenvironment datasets.

Other supporting code, scripts, or data are available for qualified groups from the authors upon request. Please contact Boxi Kang (kbxchrs#gmail.com) or Jordi Camps (jordi.camps#bayer.com) for such inquiries. 
