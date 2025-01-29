import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
import os
# Needed for some plotting
import matplotlib.pyplot as plt

# msigdb.to_csv("/projects/b1169/zzhang/enrichment_ref/msigdb.csv", index=False)
msigdb = pd.read_csv("/projects/b1169/zzhang/enrichment_ref/msigdb.csv")
progeny = dc.get_progeny(organism='human', top=500)
all_DEG_root = [
    "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG/degs/diagnosis_general/SCT/age_sex",
    "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG/degs/female_c9_hc/SCT/age"
    "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/05.DEG/degs/diagnosis/SCT/age_sex"
]
output_root = "/projects/b1042/Gate_Lab/projects/als-project/crispr_clean_final/08.enrichment/decoupleR"
lfc_thres = 0.585
min_DEG_num = 3
pathways_used = [
    'cell_type_signatures',
    'immunesigdb',
    'go_molecular_function',
    'go_biological_process',
    'go_cellular_component',
    'cell_type_signatures',
    'reactome_pathways',
    'hallmark',
    'kegg_pathways'
]


def list_files(directory):
    output = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            output.append(os.path.join(root, file))
    return output


def run_ora(pathway_names, pathway_db, DEG_df, output_prefix, min_feature_set=2):
    for cur_pathway_name in pathway_names:
        cur_output_file = output_prefix + "_" + f"{cur_pathway_name}.csv"
        pathway_db_sub = pathway_db[pathway_db['collection'] == cur_pathway_name]

        pathway_db_sub = pathway_db_sub[~pathway_db_sub.duplicated(['geneset', 'genesymbol'])]
        cur_db_enr_pvals = dc.get_ora_df(
            df=DEG_df,
            net=pathway_db_sub,
            source='geneset',
            target='genesymbol'
        )
        cur_db_enr_pvals['feature_count'] = cur_db_enr_pvals['Features'].str.split(';').apply(len)
        cur_output_plot_file = output_prefix + "_" + f"{cur_pathway_name}.pdf"
        enr_pvals_sig = cur_db_enr_pvals[(cur_db_enr_pvals["FDR p-value"] < 0.05)]
        enr_pvals_sig = enr_pvals_sig[enr_pvals_sig["feature_count"] > min_feature_set]
        if (enr_pvals_sig.shape[0] == 0):
            print(f"No significantly enriched pathway for {cur_pathway_name}. Next.")
        else:
            # dc.plot_dotplot(enr_pvals_sig, x='Combined score', y = 'Term', s='Odds ratio', c = 'FDR p-value',
            # scale = 0.25, figsize=(7,10)) plt.savefig(cur_output_plot_file) plt.close()
            enr_pvals_sig.to_csv(cur_output_file, index=False)
    return enr_pvals_sig


for DEG_root in all_DEG_root:
    all_DEG_files = list_files(DEG_root)
    parts = DEG_root.split('/')
    print(parts)
    cur_comp_level = parts[-3]
    if cur_comp_level == "female_c9_hc":
        BH_thres = 0.001
    else:
        BH_thres = 0.01
    for cur_DEG_file in all_DEG_files:
        print(f"INFO: Processing {cur_DEG_file}")
        cur_filename = os.path.basename(cur_DEG_file)
        # cur_comp_level = cur_filename.split("___")[0]
        cur_comp_str = cur_filename.split("___")[1]
        cur_comp, _ = os.path.splitext((cur_comp_str.split("__")[1]))
        cur_cell_type = cur_comp_str.split("__")[0]
        cur_out_dir = f"{output_root}/{cur_comp_level}/{cur_cell_type}"
        if not os.path.exists(cur_out_dir):
            os.makedirs(cur_out_dir)

        # running ORA
        cur_comp_ct_de_all = pd.read_csv(cur_DEG_file, index_col=0)
        cur_comp_ct_de_up = cur_comp_ct_de_all[(cur_comp_ct_de_all["BH"] < BH_thres) & \
                                               (cur_comp_ct_de_all["avg_log2FC"] > lfc_thres)]
        cur_comp_ct_de_dn = cur_comp_ct_de_all[(cur_comp_ct_de_all["BH"] < BH_thres) & \
                                               (cur_comp_ct_de_all["avg_log2FC"] < -lfc_thres)]
        cur_comp_ct_de_all['sig'] = cur_comp_ct_de_all['avg_log2FC'] * (
            -np.log10(cur_comp_ct_de_all['p_val_adj'] + 1e-3))
        if cur_comp_ct_de_up.shape[0] < min_DEG_num:
            print("{} has only {} DEG.\nSkipped overrepresentation analysis!".format(cur_filename,
                                                                                     str(cur_comp_ct_de_up.shape[0])))
            continue
        else:
            cur_comp_ct_de_prefix = cur_out_dir + "/" + os.path.splitext(cur_comp_str)[0] + "_up"
            _ = run_ora(pathway_names=pathways_used, pathway_db=msigdb, DEG_df=cur_comp_ct_de_up,
                        output_prefix=cur_comp_ct_de_prefix)
            cur_comp_ct_de_prefix = cur_out_dir + "/" + os.path.splitext(cur_comp_str)[0] + "_dn"
            _ = run_ora(pathway_names=pathways_used, pathway_db=msigdb, DEG_df=cur_comp_ct_de_dn,
                        output_prefix=cur_comp_ct_de_prefix)

        # run signaling perturbation database
        try:
            mat = cur_comp_ct_de_all[['sig']].T.rename(index={'sig': cur_cell_type})
            pathway_acts, pathway_pvals = dc.run_mlm(mat=mat, net=progeny)
        except ValueError:
            print(f"No sources with more than min_n=5 targets for {cur_cell_type}")
        print("#==========================================================#\n\n")
