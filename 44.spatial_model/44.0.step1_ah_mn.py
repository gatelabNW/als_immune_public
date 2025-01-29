import os
import scanpy as sc
import numpy as np
import pandas as pd
import scipy as sp
import torch
import cell2location as c2l
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

print("INFO: Started!")
all_sample_sc = sc.read_h5ad("/projects/b1042/Gate_Lab/projects/als-project/spatial/03.seurat_process_final/data/seurat_42.4__04.h5ad")
gw_ref = sc.read("/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/neuron_mn_ref.h5ad")
results_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/04.c2l_model/c2l1211_ah_mn"

MAP2_s = all_sample_sc[(all_sample_sc.obs['manual_3'] == 'Anterior_horns')].copy()
# excluded_samples = ['NEUDW867LF9_SCC_10.sl3___V52L19.016___C10', 'NEUHD481VCL_SCC_19.sl14___V52L06.340___HCBOT',
#                     'NEUHD481VCL_SCC_19.sl14___V52L06.340___HCBOT']
# MAP2_s_subset = MAP2_s[~MAP2_s.obs['sample'].isin(excluded_samples)].copy()


# remove mitochondria-encoded (MT) genes
MAP2_s.var['MT_gene'] = [gene.startswith('MT-') for gene in MAP2_s.var_names]

# remove MT genes for spatial mapping (keeping their counts in the object)
MAP2_s.obsm['MT'] = MAP2_s[:, MAP2_s.var['MT_gene'].values].X.toarray()
MAP2_s = MAP2_s[:, ~MAP2_s.var['MT_gene'].values]

# remove hemoglobin-encoded (HB) genes
MAP2_s.var['HB_gene'] = [gene.startswith('HB') for gene in MAP2_s.var_names]

# remove HB genes for spatial mapping (keeping their counts in the object)
MAP2_s.obsm['HB'] = MAP2_s[:, MAP2_s.var['HB_gene'].values].X.toarray()
MAP2_s = MAP2_s[:, ~MAP2_s.var['HB_gene'].values]


# Ensure same gene sets
shared_features = [
    feature for feature in MAP2_s.var_names if feature in gw_ref.var_names
]

# Make sure plot device is off
plt.close()

# Subset to shared genes
MAP2_s = MAP2_s[:, shared_features].copy()
gw_ref = gw_ref[:, shared_features].copy()

# permissive gene filtering per C2L's recommendation
selected = filter_genes(gw_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# save plot
plt.savefig(f"{results_folder}/plots/filter_genes.pdf", bbox_inches = "tight")
plt.close()

# filter the object
gw_ref = gw_ref[:, selected].copy()

# set up NB regression
c2l.models.RegressionModel.setup_anndata(adata=gw_ref,
                                         labels_key='AnnotationForDeconvolution',
                                         batch_key='sample'
                                         )

# NB reg
mod = RegressionModel(gw_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
print("INFO: Training!")
# train model
mod.train(max_epochs=500, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# Plot training loss
mod.plot_history(0)
plt.savefig(f"{results_folder}/plots/training_loss.pdf")
plt.close()
plt.cla()
plt.clf()

# Export posterior
gw_ref = mod.export_posterior(
    gw_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in gw_ref.varm.keys():
    inf_aver = gw_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                       for i in gw_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = gw_ref.var[[f'means_per_cluster_mu_fg_{i}'
                           for i in gw_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = gw_ref.uns['mod']['factor_names']
# gw_ref.__dict__['_raw'].__dict__['_var'] = gw_ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(gw_ref.var_names, inf_aver.index)
MAP2_s = MAP2_s[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# save predictions to CSV
inf_aver.to_csv(f"{results_folder}/inf_aver.csv")

# save plot - the plot device on HPC prints both plots in the same space on top of each other!!
# github issue is open but unresolved
# https://github.com/BayraktarLab/cell2location/issues/341
mod.plot_QC()
plt.savefig(f"{results_folder}/plots/training_qc.pdf")
plt.close()

# print to log file
print("INFO: Saving!")

# Save model
mod.save(f"{results_folder}/models/step1_model", overwrite=True)

# Save anndata objects with results
adata_file = f"{results_folder}/adata_objects/ref.h5ad"
gw_ref.write(adata_file)
sc_file = f"{results_folder}/adata_objects/MAP2_s_after_training.h5ad"
MAP2_s.write(sc_file)