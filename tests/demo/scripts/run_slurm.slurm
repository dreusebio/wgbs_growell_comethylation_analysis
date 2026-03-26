#!/bin/bash
#SBATCH --job-name=Cometh_CHDS_Breast_f1
#SBATCH --mail-user=gekuodza@ucdavis.edu
#SBATCH --time=23:00:00
#SBATCH --partition=high
#SBATCH --ntasks=1
#SBATCH --output=00_std_err_logs/comethyl_pearson_report.out
#SBATCH --error=00_std_err_logs/comethyl_pearson_report.err                   	
#SBATCH --cpus-per-task=59  #use 30 or 50 cores for certain things
#SBATCH --account=genome-center-grp #genome-center-grp
#SBATCH --mem=700G #Need up to 500G for certain tasks
#SBATCH --mail-type=ALL  
#SBATCH --chdir=/quobyte/lasallegrp/George/CHDS_Breast_Cancer_3_Gen/Comethyl_MK/f1

###################
# Run Information #
###################

start=`date +%s`

hostname

####################
# Load Environment #
####################


# Source your config file so conda can be used (likely .profile, could be .bashrc)

source ~/.profile
conda activate /quobyte/lasallegrp/George/envs/Comethyl_v1.3.0_BACKUP/
# conda activate /quobyte/lasallegrp/programs/.conda/envs/Comethyl_v1.3.0_UPDATED/
unset R_LIBS_USER
export R_LIBS_USER=""
export R_PROFILE_USER=/dev/null
export R_ENVIRON_USER=/dev/null

# quick sanity check (writes to your .out log)
Rscript --vanilla -e 'cat("R_LIBS_USER=", Sys.getenv("R_LIBS_USER"), "\n"); print(.libPaths()); library(SparseArray); library(comethyl); cat("OK\n")'

##############
# Comethyl Scripts#
##############

# Script 00
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label Baseline \
  --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Baseline/cytosine_reports \
  --workers 8

Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label 36-38wk \
  --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/36-38wk/cytosine_reports \
  --workers 8

Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/00_import_cpg_reports_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label Postpartum \
  --meta_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
  --cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Postpartum/cytosine_reports \
  --workers 8

# Script 01
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/01_filter_cpgs_call_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label Baseline \
  --input_rds /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/00_cpg_extraction/Baseline/Unfiltered_BSseq.rds \
  --cov 3 \
  --per_sample 0.75

#Script 02
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/02_filter_reference_regions_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label Baseline \
  --input_regions /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Regions.txt \
  --cov_min 4 \
  --meth_sd 0.05

# Script 03
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/03_build_reference_region_methylation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --dataset_label Baseline \
  --input_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
  --input_regions //quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p05/Filtered_Regions.rds

# Script 04
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/04_align_other_datasets_to_reference_cpgs_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --reference_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/01_reference_filter_regions/Baseline/cov3_75pct/Filtered_BSseq.rds \
  --dataset1_label 36-38wk \
  --dataset1_meta /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
  --dataset1_cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/36-38wk/cytosine_reports \
  --dataset2_label Postpartum \
  --dataset2_meta /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/sample_info.xlsx \
  --dataset2_cyto_dir /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_minidata/Postpartum/cytosine_reports \
  --cov 1 \
  --per_sample 0.01 \
  --workers 8  

# Script 05
Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/05_build_other_dataset_region_methylation_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
  --reference_regions /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/02_reference_region_filter/Baseline/cov3_75pct/covMin4_methSD0p05/Filtered_Regions.rds \
  --dataset1_label 36-38wk \
  --dataset1_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/04_align_to_reference_cpgs/36-38wk/Filtered_BSseq.rds \
  --dataset2_label Postpartum \
  --dataset2_bs /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/04_align_to_reference_cpgs/Postpartum/Filtered_BSseq.rds

# Script 06

Rscript scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
  --input_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/03_region_methylation/Baseline/cov3_75pct/covMin4_methSD0p08/Region_Methylation.rds \
  --sample_info /path/to/Baseline_sample_info.xlsx \
  --sample_id_col Sample_ID \
  --protected_traits_file /path/to/protected_traits.txt \
  --technical_traits_file /path/to/technical_traits.txt

  Rscript scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
  --input_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/03_region_methylation/TP36_38weeks/cov1_1pct/covMin4_methSD0p08/Region_Methylation.rds \
  --sample_info /path/to/TP36_sample_info.xlsx \
  --sample_id_col Sample_ID \
  --protected_traits_file /path/to/protected_traits.txt \
  --technical_traits_file /path/to/technical_traits.txt

  Rscript scripts/consensus/06_pc_diagnostics_consensus.R \
  --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
  --input_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/03_region_methylation/TP36_38weeks/cov1_1pct/covMin4_methSD0p08/Region_Methylation.rds \
  --sample_info /path/to/TP36_sample_info.xlsx \
  --sample_id_col Sample_ID \
  --protected_traits_file /path/to/protected_traits.txt \
  --technical_traits_file /path/to/technical_traits.txt