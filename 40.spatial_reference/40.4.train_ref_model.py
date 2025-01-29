import scanpy as sc
import numpy as np
import cell2location as c2l
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import pandas as pd
import os

results_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/scanpy_obj"
barcodes_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/03.seurat_process/barcodes"
print("INFO: Started!")

ref = sc.read("/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/ref.h5ad")
root_folder = "/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v4_0404"
search_pattern = os.path.join(root_folder, '**', 'outs')
subfolders = [path for path in glob.glob(search_pattern, recursive=True)]

print("INFO: Started!")

ref.X = ref.X.astype('int')
selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
ref = ref[:, selected].copy()
# set up NB regression
c2l.models.RegressionModel.setup_anndata(adata=ref,
                                         # cell type, covariate used for constructing signatures
                                         labels_key='AnnotationForDeconvolution',
                                         batch_key='sample'
                                         )

mod = RegressionModel(ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()

print("INFO: Training!")
mod.train(max_epochs=500, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)
plt.close()
mod.plot_history(20)
plt.savefig(f"{results_folder}/plots/training_loss.png")
plt.close()

print("INFO: Saving!")
# Save model
mod.save(f"{results_folder}/model.pt", overwrite=True)

ref = mod.export_posterior(
    ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

mod.plot_QC()
plt.savefig(f"{results_folder}/plots/training_qc.png")
plt.close()

if 'means_per_cluster_mu_fg' in ref.varm.keys():
    inf_aver = ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                                    for i in ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = ref.var[[f'means_per_cluster_mu_fg_{i}'
                        for i in ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = ref.uns['mod']['factor_names']
ref.__dict__['_raw'].__dict__['_var'] = ref.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
ref.write(f"{results_folder}/ref_processed.h5ad")
inf_aver.to_csv(f"{results_folder}/inf_aver.csv")

for cur_folder in subfolders:
    parts = cur_folder.split('/')
    ID = parts[-2]
    print(ID)
    ID_no_dash = ID.replace("-", ".")
    # get post QC barcodes
    cur_sample_bc_file = f"{barcodes_folder}/{ID_no_dash}.csv"
    if not os.path.exists(cur_sample_bc_file):
        print(f"Skip {cur_sample_bc_file}.")
        continue
    cur_post_qc_barcodes_df = pd.read_csv(cur_sample_bc_file)
    cur_post_qc_barcodes = cur_post_qc_barcodes_df["barcode"].tolist()
    all_sample_sc = sc.read_visium(cur_folder)
    all_sample_sc = all_sample_sc[cur_post_qc_barcodes, :]

    # remove mitochondria-encoded (MT) genes
    all_sample_sc.var['MT_gene'] = [gene.startswith('MT-') for gene in all_sample_sc.var_names]

    # remove MT genes for spatial mapping (keeping their counts in the object)
    all_sample_sc.obsm['MT'] = all_sample_sc[:, all_sample_sc.var['MT_gene'].values].X.toarray()
    all_sample_sc = all_sample_sc[:, ~all_sample_sc.var['MT_gene'].values]

    # remove hemoglobin-encoded (HB) genes
    all_sample_sc.var['HB_gene'] = [gene.startswith('HB') for gene in all_sample_sc.var_names]

    # remove HB genes for spatial mapping (keeping their counts in the object)
    all_sample_sc.obsm['HB'] = all_sample_sc[:, all_sample_sc.var['HB_gene'].values].X.toarray()
    all_sample_sc = all_sample_sc[:, ~all_sample_sc.var['HB_gene'].values]

    all_sample_sc.X = all_sample_sc.X.astype('int')
    shared_features = [
        feature for feature in all_sample_sc.var_names if feature in ref.var_names
    ]
    print("{} genes shared!".format(str(len(shared_features))))

    plt.savefig(f"{results_folder}/plots/filter_genes.png")
    plt.close()

    intersect = np.intersect1d(all_sample_sc.var_names, inf_aver.index)
    all_sample_sc.var_names_make_unique()
    all_sample_sc = all_sample_sc[:, intersect].copy()
    inf_aver_cur_sample = inf_aver.loc[intersect, :].copy()
    inf_aver_cur_sample.to_csv(f"{results_folder}/{ID}__inf_aver.csv")

    sc_file = f"{results_folder}/{ID}.h5ad"
    all_sample_sc.write(sc_file)
