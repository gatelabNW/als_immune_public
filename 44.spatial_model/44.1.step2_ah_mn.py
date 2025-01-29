# Written by: Thomas Watson
# Summary: Cell2Location Step 2
#
#-------------------------------------------------------------------------------

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
results_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/04.c2l_model/c2l1211_ah_mn"
run_name = f'{results_folder}/Cell2Location_res'

all_sample_sc = sc.read_h5ad(f"{results_folder}/adata_objects/MAP2_s_after_training.h5ad")

# load training posterior probabilities
inf_aver = pd.read_csv(f"{results_folder}/inf_aver.csv")

# index vals saved in column 1 but not index so make sure inf_aver.index = gw_ref_downsamled.var_names before running the model
inf_aver.index = inf_aver.iloc[:]['Unnamed: 0'].astype("str")
inf_aver.drop(columns=inf_aver.columns[0], axis=1,  inplace=True)
all_sample_sc.var_names = all_sample_sc.var_names.astype("str")

# prepare anndata for cell2location model
c2l.models.Cell2location.setup_anndata(adata=all_sample_sc, batch_key="sample", continuous_covariate_keys =["cdr_centered", "nCount_Spatial"])

mod = c2l.models.Cell2location(
    all_sample_sc, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=7,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)

mod.train(max_epochs=5000,
          # train using full data (batch_size=None)
          batch_size=5000,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
          )

# plot ELBO loss history during training
plt.close()
mod.plot_history(0)
plt.legend(labels=['full data training'])
plt.savefig(f"{results_folder}/plots/step2_training_loss.pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()
###########################################################################################################

all_sample_sc = mod.export_posterior(
    #all_sample_sc, sample_kwargs={'num_samples': 1000, 'batch_size': 5000, 'use_gpu': True}
    all_sample_sc, sample_kwargs={'batch_size': 2500, 'use_gpu': True},
    add_to_obsm=["q05","q50", "q95", "q0001"], use_quantiles=True
)

# Save model
mod.save(f"{results_folder}/models/step2_model", overwrite=True)
plt.close()
mod.plot_QC(summary_name = "q05")
plt.savefig(f"{results_folder}/plots/step2_training_QC.pdf", bbox_inches = "tight")
plt.close()
#############################################################################################################

all_sample_sc.obs[all_sample_sc.uns['mod']['factor_names']] = all_sample_sc.obsm['q05_cell_abundance_w_sf']

# Save anndata object with results
adata_file = f"{results_folder}/adata_objects/q05final.h5ad"

all_sample_sc.write(adata_file)