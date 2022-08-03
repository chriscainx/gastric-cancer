# Use ICA to identify the batch component in global Pearson Residuals

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import bbknn
from matplotlib import rcParams
from ica import ica1

sc.set_figure_params(dpi=160, color_map='viridis', fontsize=6)
sc.settings.verbosity=2
sc.settings.n_jobs=64

adata = sc.read_h5ad("gcall.correctedumi.h5ad")
pears_res = np.load("residuals.npy")

# Use Pearson residuals for ICA

n_ics = 128
A,S,W = ica1(pois_res.T, n_ics)

# All Independent Components are tested for correlation against the a HSP gene, for quick identification of the 
# dissociation-driven component

hsp = np.array(adata[:, 'HSPA1A'].X)

def modcor(b):
    c, p = pearsonr(hsp, b)
    return abs(c)

with Pool(64) as p:
    corr = p.map(modcor, [A[:, i] for i in range(n_ics)])
max(corr)

# Max correlation was identified for Independent Component 15, in our case.

adata = sc.read_h5ad("check/02.correctedumi.h5ad")
adata.obs['IC_HSP'] = A[:, 15]

# Identify the top-weighted genes in IC_15

adata.var['HSP_loading'] = S[15,:]
adata.var.nlargest(100, 'HSP_loading')

# --------------------------------------------------------------------------------
#           | gene_ids          |  n_cells   |  highly_variable   |  HSP_loading
# index                        
# --------------------------------------------------------------------------------
# HSPA6     | ENSG00000173110   |  9942.0    |  True              |  100.647754
# HSPA1A    | ENSG00000204389   |  50455.0   |  True              |  64.528467
# HSPA1B    | ENSG00000204388   |  48008.0   |  True              |  37.909546
# DNAJB1    | ENSG00000132002   |  66949.0   |  True              |  29.365036
# HSP90AA1  | ENSG00000080824   |  91604.0   |  True              |  27.476427
# HSPH1     | ENSG00000120694   |  34183.0   |  True              |  24.872009
# HSPB1     | ENSG00000106211   |  48888.0   |  True              |  20.291679
# IGLL5     | ENSG00000254709   |  16041.0   |  True              |  18.908379
# HSPD1     | ENSG00000144381   |  44129.0   |  True              |  18.229984
# BAG3      | ENSG00000151929   |  16355.0   |  True              |  15.629638
# HSPE1     | ENSG00000115541   |  56849.0   |  True              |  15.141479
# ZFAND2A   | ENSG00000178381   |  15781.0   |  True              |  13.962453
# DNAJA1    | ENSG00000086061   |  56059.0   |  True              |  11.626205
# HSPA8     | ENSG00000109971   |  78454.0   |  True              |  11.282519
# DNAJA4    | ENSG00000140403   |  7575.0    |  True              |  10.991242
# HSP90AB1  | ENSG00000096384   |  76211.0   |  True              |  10.394976
# CACYBP    | ENSG00000116161   |  31686.0   |  True              |  10.151749
# --------------------------------------------------------------------------------

# Remove loading of the IC_15 from the pearson residuals (SCTransform result)

IC_batch = np.matmul(A[:,[15]], S[[15],:])
pears_res = pears_res - IC_batch.T
np.save("residual.corrected.npy", pears_res)

# Scale the single cell data to remove impact of IC_15 on all gene expressions
# Select highly variable genes using SCTransform result

adata.var['variance'] = np.var(pears_res, 1)
sc.pp.log1p(adata)
adata.raw = adata
adata = adata[:, adata.var.nlargest(1000, 'variance').index]
sc.pp.regress_out(adata, ['percent_mito', 'IC_HSP'])
adata.X = adata.X - np.mean(adata.X, axis=0)

# Data correction complete.
adata.write("02.ICA_corrected.h5ad")

