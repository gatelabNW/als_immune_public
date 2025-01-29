import os
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
import scanpy as sc
results_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref"
output_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/04.c2l_model"
all_h5ad_files = [file for file in os.listdir(results_folder) if file.endswith("h5ad") if file != "ref.h5ad" and file != 'ref_processed.h5ad']

for cur_h5ad in all_h5ad_files:
    cur_h5ad_id, _ = os.path.splitext(cur_h5ad)
    sc_file = f"{output_folder}/anndata/{cur_h5ad_id}__after_training.h5ad"
    if os.path.exists(sc_file):
        print(f"File exists: {sc_file}, skipping...")
        continue
    dat = sc.read_h5ad(f"{results_folder}/{cur_h5ad}")
    inf_aver = pd.read_csv(f"{results_folder}/inf_aver.csv")

    inf_aver.index = inf_aver.iloc[:]['Unnamed: 0'].astype("str")
    inf_aver.drop(columns=inf_aver.columns[0], axis=1, inplace=True)
    dat.var_names = dat.var_names.astype("str")

    c2l.models.Cell2location.setup_anndata(adata=dat)
    mod = c2l.models.Cell2location(
        dat, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=7,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )

    mod.train(
        max_epochs=15000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
        use_gpu=True,
    )
    plt.close()
    mod.plot_history(1000)
    plt.legend(labels=['full data training'])
    plt.savefig(f"{output_folder}/plots/{cur_h5ad_id}__ss_curve.png")
    plt.close()

    dat = mod.export_posterior(
        dat, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )

    mod.save(f"{output_folder}/model/{cur_h5ad_id}__model.pt", overwrite = True)

    dat.write(sc_file)