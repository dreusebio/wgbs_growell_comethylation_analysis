#!/bin/bash
#SBATCH --job-name=Cometh_growell
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --time=14-00:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err                   	
#SBATCH --cpus-per-task=59  #use 30 or 50 cores for certain things
#SBATCH --account=genome-center-grp #genome-center-grp
#SBATCH --mem=950G #Need up to 500G for certain tasks
#SBATCH --mail-type=ALL  
#SBATCH --chdir=/quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George

###################
# Run Information #
###################

start=`date +%s`

hostname
####################
# Load Environment #
####################

# Initialize conda (works in non-interactive SLURM shells)
source ~/.bashrc

# Activate the conda environment
conda activate Comethyl_v1.3.0

# Unset R user library paths
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# Quick sanity check (writes to your .out log)
Rscript --vanilla -e 'cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n"); print(.libPaths()); library(SparseArray); library(comethyl); cat("OK\n")'

##############
# Comethyl Scripts#
##############

# # Script 00
# pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset_label Baseline \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_baseline.xlsx \
#   --cyto_dir  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/processed/Baseline \
#   --workers 8

#  pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset_label 36-38wk \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_36_38wk.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/processed/36-38wk \
#   --workers 8

# pixi run Rscript  scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset_label Postpartum \
#   --meta_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_postpartum.xlsx \
#   --cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/processed/Postpartum \
#   --workers 8

# # Script 01
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/01_filter_cpgs_call_regions_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset_label Baseline \
#   --input_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/00_cpg_extraction/Baseline/Unfiltered_BSseq.rds \
#   --cov 3 \
#   --per_sample 0.75

# Script 02
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/02_filter_reference_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --dataset_label Baseline \
  --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Regions.txt \
  --cov_min 20 \
  --meth_sd 0.06

# Script 03
pixi run Rscript /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/03_build_reference_region_methylation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --dataset_label Baseline \
  --input_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
  --input_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin20_methSD0p06/Filtered_Regions.rds

# Script 04
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --reference_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
  --dataset1_label 36-38wk \
  --dataset1_meta /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_36_38wk.xlsx \
  --dataset1_cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/processed/36-38wk \
  --dataset2_label Postpartum \
  --dataset2_meta  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_postpartum.xlsx \
  --dataset2_cyto_dir /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/processed/Postpartum \
  --cov 1 \
  --per_sample 0.01 \
  --workers 8  

# Script 05
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --reference_regions /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin20_methSD0p06/Filtered_Regions.rds \
  --dataset1_label 36-38wk \
  --dataset1_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/04_align_to_reference_cpgs/36-38wk/Filtered_BSseq.rds \
  --dataset2_label Postpartum \
  --dataset2_bs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/04_align_to_reference_cpgs/Postpartum/Filtered_BSseq.rds

# Script 05b
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/05b_filter_shared_complete_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --dataset1_label Baseline \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/03_region_methylation/Baseline/cov3_75pct/covMin20_methSD0p06/Region_Methylation.rds \
  --dataset2_label 36-38wk \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/03_region_methylation/36-38wk/36-38wk/covMin20_methSD0p06/Region_Methylation.rds \
  --dataset3_label Postpartum \
  --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/03_region_methylation/Postpartum/Postpartum/covMin20_methSD0p06/Region_Methylation.rds \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin20_methSD0p06/Filtered_Regions.txt

# Script 06

pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/Baseline_Methylation_complete.rds \
  --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_baseline.xlsx \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --selected_trait_mode   protected_and_technical

pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/36-38wk_Methylation_complete.rds \
  --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_36_38wk.xlsx \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --selected_trait_mode   protected_and_technical
  
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/Postpartum_Methylation_complete.rds \
  --sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_postpartum.xlsx \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --selected_trait_mode   protected_and_technical

#Script 07
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/Baseline_Methylation_complete.rds \
  --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Baseline/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Baseline/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 TRUE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05

pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/36-38wk_Methylation_complete.rds \
  --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/36-38wk/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/36-38wk/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 TRUE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05

pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/07_methylation_adjustment_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --input_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/Postpartum_Methylation_complete.rds \
  --input_pcs /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Postpartum/PCs.rds \
  --pc_trait_stats /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/06_pc_diagnostics/Postpartum/PC_Trait_Correlation_Stats_Bicor.tsv \
  --protected_traits_file config/protected_traits.txt \
  --technical_traits_file config/technical_traits.txt \
  --run_v1 TRUE \
  --run_v2 TRUE \
  --run_v3 TRUE \
  --v2_cor_method bicor \
  --v2_p_thresh 0.05

# # script 08
# #v1
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset1_label Baseline \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v1_all_pcs/Baseline_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label 36-38wk \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v1_all_pcs/36-38wk_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset3_label Postpartum \
#   --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v1_all_pcs/Postpartum_Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59 \
#   --scale_free_cutoff 0.8

#v2
pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --dataset1_label Baseline \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v2_exclude_protected_pcs/Baseline_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label 36-38wk \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v2_exclude_protected_pcs/36-38wk_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset3_label Postpartum \
  --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v2_exclude_protected_pcs/Postpartum_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --softpower_cor pearson \
  --power_min 1 \
  --power_max 20 \
  --threads 59 \
  --scale_free_cutoff 0.8

# #v3
#pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/08_softpower_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --dataset1_label Baseline \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v3_technical_pcs_only/Baseline_Adjusted_Region_Methylation_technicalPCs_only.rds \
#   --dataset2_label 36-38wk \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v3_technical_pcs_only/36-38wk_Adjusted_Region_Methylation_technicalPCs_only.rds \
#   --dataset3_label Postpartum \
#   --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v3_technical_pcs_only/Postpartum_Adjusted_Region_Methylation_technicalPCs_only.rds \
#   --adjustment_version v3_technical_pcs_only \
#   --softpower_cor pearson \
#   --power_min 1 \
#   --power_max 20 \
#   --threads 59 \
#   --scale_free_cutoff 0.8

# #script 09
# #v1
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
#   --dataset1_label Baseline \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v1_all_pcs/Baseline_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label 36-38wk \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v1_all_pcs/36-38wk_Adjusted_Region_Methylation_allPCs.rds \
#   --dataset3_label Postpartum \
#   --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v1_all_pcs/Postpartum_Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/08_soft_power/v1_all_pcs/chosen_power_pearson.txt \
#   --consensus_cor pearson \
#   --network_type signed \
#   --tom_type signed \
#   --deep_split 1 \
#   --min_module_size 10 \
#   --merge_cut_height 0.2 \
#   --max_block_size 1000 \
#   --network_calibration "single quantile" \
#   --save_consensus_toms TRUE \
#   --threads 59 \
#   --max_p_outliers 0.1

#v2

pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
  --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
  --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
  --dataset1_label Baseline \
  --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v2_exclude_protected_pcs/Baseline_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset2_label 36-38wk \
  --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v2_exclude_protected_pcs/36-38wk_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --dataset3_label Postpartum \
  --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v2_exclude_protected_pcs/Postpartum_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
  --adjustment_version v2_exclude_protected_pcs \
  --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/08_soft_power/v2_exclude_protected_pcs/chosen_power_pearson.txt \
  --consensus_cor pearson \
  --network_type signed \
  --tom_type signed \
  --deep_split 1 \
  --min_module_size 10 \
  --merge_cut_height 0.2 \
  --max_block_size 50000 \
  --network_calibration "single quantile" \
  --save_consensus_toms TRUE \
  --threads 16 \
  --max_p_outliers 0.1

# #v3

# # pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/09_modules_consensus.R \
# #   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
# #   --regions_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/shared_complete_regions_4col.tsv \
# #   --dataset1_label Baseline \
# #   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v3_technical_pcs_only/Baseline_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --dataset2_label 36-38wk \
# #   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/05b_shared_complete_regions/36-38wk_Methylation_complete.rds \
# #   --dataset3_label Postpartum \
# #   --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v3_technical_pcs_only/Postpartum_Adjusted_Region_Methylation_technicalPCs_only.rds \
# #   --adjustment_version v3_technical_pcs_only \
# #   --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/08_soft_power/v3_technical_pcs_only/chosen_power_pearson.txt \
# #   --consensus_cor pearson \
# #   --network_type signed \
# #   --tom_type signed \
# #   --deep_split 1 \
# #   --min_module_size 10 \
# #   --merge_cut_height 0.1 \
# #   --max_block_size 1000 \
# #   --network_calibration "single quantile" \
# #   --save_consensus_toms TRUE \
# #   --threads 59 \
# #   --max_p_outliers 0.1

# # --chosen_power_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/08_soft_power/v2_exclude_protected_pcs/chosen_power_pearson.txt \
# # --power 15

# # #Script 10
# #v2
# pixi run Rscript  scripts/comethyl_scripts/consensus/10_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_outcome_exposure_pcs \
#   --consensus_cor pearson \
#   --max_p_outliers 0.1

# # #Script 11
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/11_membership_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --dataset1_label Baseline \
#   --dataset1_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Baseline/v2_exclude_protected_pcs/Baseline_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset2_label 36-38wk \
#   --dataset2_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/36-38wk/v2_exclude_protected_pcs/36-38wk_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --dataset3_label Postpartum \
#   --dataset3_meth /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/07_methylation_adjustment/Postpartum/v2_exclude_protected_pcs/Postpartum_Adjusted_Region_Methylation_excluding_protected_PCs_bicor.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --membership_cor bicor \
#   --max_p_outliers 0.1

# #script 12
# #v2
# pixi run Rscript  /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/12a_me_trait_analysis_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --dataset1_label Baseline \
#   --dataset1_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_baseline.xlsx \
#   --dataset2_label 36-38wk \
#   --dataset2_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_36_38wk.xlsx \
#   --dataset3_label Postpartum \
#   --dataset3_sample_info /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/data/metadata/sample_info_postpartum.xlsx \
#   --run_bicor TRUE \
#   --run_pearson FALSE \
#   --p_thresh 0.05 \
#   --top_n 50 \
#   --module_dendro_distance bicor \
#   --max_p_outliers 0.1

# # Script 12b
# pixi run Rscript   /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/scripts/comethyl_scripts/consensus/12b_me_trait_presentation_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George \
#   --consensus_modules_rds /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/09_consensus_modules/v2_exclude_protected_pcs/shared/Consensus_Modules.rds \
#   --adjustment_version v2_exclude_protected_pcs \
#   --dataset1_label Baseline \
#   --dataset1_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/v2_exclude_protected_pcs/Baseline/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --dataset2_label 36-38wk \
#   --dataset2_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/v2_exclude_protected_pcs/36-38wk/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --dataset3_label Postpartum \
#   --dataset3_stats_file /quobyte/lasallegrp/projects/GROWELL/WGBS/2025_DBS_comethyl_George/comethyl_output/consensus/12a_me_trait_analysis/v2_exclude_protected_pcs/Postpartum/ME_Trait_Correlation_Stats_Bicor.tsv \
#   --set_dir config/presentation_sets \
#   --module_dendro_distance bicor \
#   --p_thresh 0.05 \
#   --top_n 50 \
#   --full_width 12 \
#   --full_height 12 \
#   --top_width 9 \
#   --top_height 6