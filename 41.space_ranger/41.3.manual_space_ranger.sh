#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomicsguestA
#SBATCH --job-name spaceranger
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 128G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1042/Gate_Lab/zzhang/cluster_logs/%x_oe%j.log
#SBATCH --verbose





#id="GWF_19-47_10B_sl13"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/GWF_19-47_10B_sl13___V52L19-048___HCBOT/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/CAVG10001_2023-10-13_09-54-40_2023-10-13_09-36-08_V52L19-048_A_GWF_19-47_10B_sl13.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/GWF19-47_multipage_small.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/V52L19-048-A1.IF.json"
#slide="V52L19-048"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/GWF_19-47_10B_sl13___V52L19-048___HCBOT"



#id="GWF_19-47_10B_sl13"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/GWF_19-47_10B_sl13___V52L19-048___HCTOP/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCTOP/CAVG10001_2023-10-13_09-54-40_2023-10-13_09-36-08_V52L19-048_A_GWF_19-47_10B_sl13.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCTOP/GWF19-47_multipage_small.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCTOP/V52L19-048-A1.IF.json"
#slide="V52L19-048"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/GWF_19-47_10B_sl13___V52L19-048___HCTOP"



#id="NEUHD481VCL_SCC_19-sl14"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/NEUHD481VCL_SCC_19-sl14___V52L06-340___HCBOT/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCBOT/CAVG10001_2023-09-29_10-03-59_2023-09-29_09-33-50_V52L06-340_A_NEUHD481VCL_SCC_19-sl14.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCBOT/NEUDH481VCL-large-multipage_small.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCBOT/V52L06-340-A1.IF.json"
#slide="V52L06-340"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/NEUHD481VCL_SCC_19-sl14___V52L06-340___HCBOT"



#id="NEUHD481VCL_SCC_19-sl14"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/NEUHD481VCL_SCC_19-sl14___V52L06-340___HCTOP/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCTOP/CAVG10001_2023-09-29_10-03-59_2023-09-29_09-33-50_V52L06-340_A_NEUHD481VCL_SCC_19-sl14.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCTOP/NEUDH481VCL-small_multipage_small.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_NEUHD481VCL_SCC_19-sl14___V52L06-340___1___A___HCTOP/V52L06-340-A1.IF.json"
#slide="V52L06-340"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/NEUHD481VCL_SCC_19-sl14___V52L06-340___HCTOP"


#id="GWF_19-47_10B_sl13"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_original_final/GWF_19-47_10B_sl13___V52L19-048___HCBOT/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/CAVG10001_2023-10-13_09-54-40_2023-10-13_09-36-08_V52L19-048_A_GWF_19-47_10B_sl13.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/GWF19-47_original_8bit_map2_only.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/V52L19-048-A1.8bit_map2_only.original.json"
#slide="V52L19-048"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_original_final/GWF_19-47_10B_sl13___V52L19-048___HCBOT"

#id="GWF_19-47_10B_sl13"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_original_final/GWF_19-47_10B_sl13___V52L19-048___HCBOT/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/CAVG10001_2023-10-13_09-54-40_2023-10-13_09-36-08_V52L19-048_A_GWF_19-47_10B_sl13.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/GWF19-47_original_8bit_map2_only.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF_19-47_10B_sl13___V52L19-048___1___A___HCBOT/V52L19-048-A1.8bit_map2_only.v4.original.json"
#slide="V52L19-048"
#area="A1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_original_final/GWF_19-47_10B_sl13___V52L19-048___HCBOT/effective_spot_selection"



# TODO: These two samples are reran manually due to error
#id="NEUDW867LF9_SCC_10-sl3"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v6_original_0801/NEUDW867LF9_SCC_10-sl3___V52L19-016___C10/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEUDW867LF9_SCC_10-sl3___V52L19-016___2___B___HC/CAVG10001_2023-09-22_10-41-55_2023-09-22_10-30-01_V52L19-016_B_NEUDW867LF9_SCC_10-sl3.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEUDW867LF9_SCC_10-sl3___V52L19-016___2___B___HC/H2_NLF9_restained-multipage-no-motor-neuron_binned.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEUDW867LF9_SCC_10-sl3___V52L19-016___2___B___HC/V52L19-016-B1.original.json"
#slide="V52L19-016"
#area="B1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v6_original_0801/NEUDW867LF9_SCC_10-sl3___V52L19-016___C10"




#id="NEULM733WR7-10-sl3"
#transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
#probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
#feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
#lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v6_original_0801/NEULM733WR7SCC-10-sl3___V52L18-387___C6/lib.csv"
#cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEULM733WR7SCC-10-sl3___V52L18-387___2___B___HC/CAVG10001_2023-09-15_11-11-46_2023-09-15_11-04-57_V52L18-387_B_NEULM733WR7SCC-10-sl3.tif"
#dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEULM733WR7SCC-10-sl3___V52L18-387___2___B___HC/NEULM733WR7_DAPI20j_GFP20j_RFP20j_AF64720j_Seq0000-binned.tif"
#json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/B_NEULM733WR7SCC-10-sl3___V52L18-387___2___B___HC/V52L18-387-B1.original.json"
#slide="V52L18-387"
#area="B1"
#output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v6_original_0801/NEULM733WR7SCC-10-sl3___V52L18-387___C6"


id="GWF-18-33-10B-sl11"
transcriptome_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/refdata-gex-GRCh38-2020-A"
probe_set="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv"
feature_ref="/projects/b1042/Gate_Lab/projects/als-project/spatial/00.ref/space_ranger/Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv"
lib="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_v9_IF_1004/GWF-18-33-10B-sl11___V52L19-028___C7/lib.csv"
cyt_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF-18-33-10B-sl11___V52L19-028___1___A___HC/CAVG10001_2023-09-22_09-52-30_2023-09-22_09-33-25_V52L19-028_A_GWF-18-33-10B-sl11.tif"
dk_img="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF-18-33-10B-sl11___V52L19-028___1___A___HC/S3_GW18.binned.tif"
json="/projects/b1042/Gate_Lab/projects/als-project/spatial/01.spaceranger_input_processed_v3_dual_alignment/A_GWF-18-33-10B-sl11___V52L19-028___1___A___HC/V52L19-028-A1.original.json"
slide="V52L19-028"
area="A1"
output_dir="/projects/b1042/Gate_Lab/projects/als-project/spatial/02.spaceranger_output_original_final/GWF-18-33-10B-sl11___V52L19-028___C7"



# Navigate to output directory
cd $output_dir

# Run spaceranger
/projects/b1169/zzhang/software/spaceranger-2.1.1/bin/spaceranger count \
--id="${id}" \
--transcriptome="${transcriptome_ref}" \
--probe-set="${probe_set}" \
--feature-ref="${feature_ref}" \
--libraries="${lib}" \
--cytaimage="${cyt_img}" \
--darkimage="${dk_img}" \
--loupe-alignment="${json}" \
--slide="${slide}" \
--area="${area}" \
--localcores=16 \
--localmem=128


#--colorizedimage="${dk_img}" \

# use colorized image for that one special sample me and francesco couldn't figure out why : GWF_19-47_10B_sl13___V52L19-048___HCBOT
