#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 06: PC Derivation and PC-Trait Diagnostics for Consensus
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   - Load region-level methylation from consensus Script 03, 05,
#     or 05b (recommended when analyzing multiple timepoints)
#   - Load analysis-ready sample information
#   - Align samples between methylation and sample info
#   - Derive top PCs using sva::num.sv + SVD
#   - Compute PC-trait association statistics (Pearson and bicor)
#   - Generate PC-trait diagnostic heatmaps and trait dendrogram
#   - Save reusable outputs for downstream methylation adjustment
#
# REQUIRED INPUTS
#   --project_root : root directory of the analysis project
#   --input_meth   : path to Region_Methylation.rds from consensus
#                    Script 03, 05, or 05b. When analyzing multiple
#                    timepoints, use the NA-filtered matrices from
#                    Script 05b (e.g. Baseline_Methylation_complete.rds)
#   --sample_info  : path to sample information file (.xlsx, .csv, .tsv, .txt)
#
# OPTIONAL INPUTS
#   --sample_id_col         : column in sample_info containing sample IDs
#                             [default = use existing rownames]
#   --protected_traits_file : text file with one protected trait per line
#   --technical_traits_file : text file with one technical trait per line
#   --threads               : number of WGCNA threads [default = 4]
#   --top_trait_n           : number of top traits for overall heatmap [default = 50]
#   --top_cell_n            : number of top PC-trait cells for association heatmap [default = 50]
#   --primary_cor_method    : primary correlation method for diagnostic plots:
#                             bicor or pearson [default = bicor]
#   --selected_trait_mode   : protected, protected_and_technical, none [default = auto]
#   --plot_width            : plot width in inches [default = 12]
#   --plot_height           : plot height in inches [default = 10]
#   --plot_dpi              : plot DPI [default = 600]
#   --base_size             : base text size for plots [default = 11]
#   --axis_text_size        : axis text size for plots [default = 8]
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/06_pc_diagnostics/<dataset_label>/
#       PCs.rds
#       PC_Metadata.tsv
#       Variance_Explained_First20_PCs.pdf
#       PC_Trait_Correlation_Stats_Pearson.tsv
#       PC_Trait_Correlation_Stats_Bicor.tsv
#       PC_Trait_Bicor_Heatmap_PUBLICATION.pdf   (if selected traits are available)
#       PC_Trait_Bicor_Heatmap_Top<n>Cells.pdf
#       Trait_Dendrogram_From_PC_Associations.pdf
#       sample_info_non_numeric_columns_removed.txt
#       sample_info_all_na_numeric_columns_removed.txt
#       sample_info_zero_variance_numeric_columns_removed.txt
#       traits_requested_protected.txt           (if provided)
#       traits_found_protected.txt               (if provided)
#       traits_missing_protected.txt             (if provided)
#       traits_requested_technical.txt           (if provided)
#       traits_found_technical.txt               (if provided)
#       traits_missing_technical.txt             (if provided)
#       run_parameters.txt
#       sessionInfo.txt
#
# NOTES
#   - sample_info is assumed to be analysis-ready before entering the pipeline.
#   - No project-specific trait recoding is performed here.
#   - If protected traits are provided, they are used in the SVA model matrix
#     for num.sv(). If not provided, an intercept-only model is used.
#   - This script should be run separately for each dataset/timepoint in the
#     consensus workflow.
#   - The methylation matrix passed to --input_meth must be free of NA values.
#     sva::num.sv() and SVD-based PCA cannot handle missing values. Imputation
#     with 0 is not appropriate for WGBS data where NAs reflect insufficient
#     coverage rather than zero methylation.
#   - When analyzing multiple timepoints, run consensus Script 05b BEFORE this
#     script to identify regions that are complete (NA-free) across all
#     timepoints. Use the filtered matrices output by 05b as --input_meth here:
#         <project_root>/comethyl_output/consensus/05b_shared_complete_regions/
#             <dataset_label>_Methylation_complete.rds
#   - If analyzing a single timepoint in isolation, complete cases can be
#     extracted directly from the Region_Methylation.rds produced by Script 03
#     or 05, but using 05b is still recommended for consistency and
#     reproducibility across timepoints.
#   - Output is written to 06_pc_diagnostics/<dataset_label>/ where
#     dataset_label is automatically derived from the input methylation path.
#
# RECOMMENDED PIPELINE ORDER (multi-timepoint)
#   01_reference_filter_regions         ->
#   02_reference_region_filter          ->
#   03_region_methylation               ->  (reference dataset)
#   04_align_to_reference_cpgs          ->  (non-reference datasets)
#   05_build_other_dataset_...          ->  (non-reference datasets)
#   05b_filter_shared_complete_regions  ->  (all timepoints together)
#   06_pc_diagnostics                   ->  (run separately per timepoint)

# Example
# Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --input_meth /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/05b_shared_complete_regions/Baseline_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/metadata/sample_info_baseline.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/technical_traits.txt \
#   --selected_trait_mode   protected

#   Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --input_meth /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/05b_shared_complete_regions/36-38wk_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/metadata/sample_info_36_38wk.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/technical_traits.txt \
#   --selected_trait_mode   protected
  
#   Rscript /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/scripts/comethyl_scripts/consensus/06_pc_diagnostics_consensus.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo \
#   --input_meth /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/comethyl_output/consensus/05b_shared_complete_regions/Postpartum_Methylation_complete.rds \
#   --sample_info /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/metadata/sample_info_postpartum.xlsx \
#   --protected_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/protected_traits.txt \
#   --technical_traits_file /quobyte/lasallegrp/projects/wgbs_growell_comethylation_analysis/tests/demo/config/technical_traits.txt \
#   --selected_trait_mode   protected
# ================================================================

message("Starting Script 6 ✓")
suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(comethyl)
  library(sva)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(gridExtra)
  library(grid)
})

# ------------------------------------------------------------
# Load helper functions
# ------------------------------------------------------------
script_file_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
if (length(script_file_arg) == 0) {
  stop("Could not determine script path from commandArgs().")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_file_arg[1])))
helper_file <- file.path(script_dir, "helper.R")
if (!file.exists(helper_file)) {
  stop("helper.R not found next to this script: ", helper_file)
}
source(helper_file)

# ------------------------------------------------------------
# Parse command-line arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character",
              help = "Root directory of the project"),

  make_option("--input_meth", type = "character",
              help = "Path to Region_Methylation.rds from consensus Script 03, 05, or 05b"),

  make_option("--sample_info", type = "character",
              help = "Path to sample information file (.xlsx, .csv, .tsv, .txt)"),

  make_option("--sample_id_col", type = "character", default = NULL,
              help = "Column in sample_info containing sample IDs [default = rownames]"),

  make_option("--protected_traits_file", type = "character", default = NULL,
              help = "Optional text file with one protected trait per line"),

  make_option("--technical_traits_file", type = "character", default = NULL,
              help = "Optional text file with one technical trait per line"),

  make_option("--threads", type = "integer", default = 4,
              help = "Number of WGCNA threads [default = 4]"),

  make_option("--top_trait_n", type = "integer", default = 50,
              help = "Number of top traits for overall heatmap [default = 50]"),

  make_option("--top_cell_n", type = "integer", default = 50,
              help = "Number of top PC-trait cells for association heatmap [default = 50]"),

  make_option("--primary_cor_method", type = "character", default = "bicor",
              help = "Primary correlation method for diagnostics: bicor or pearson [default = bicor]"),

  make_option("--selected_trait_mode", type = "character", default = NULL,
              help = "Trait mode for selected-trait heatmap: protected, protected_and_technical, or none [default = auto]"),

  make_option("--plot_width", type = "double", default = 12,
              help = "Plot width in inches [default = 12]"),

  make_option("--plot_height", type = "double", default = 10,
              help = "Plot height in inches [default = 10]"),

  make_option("--plot_dpi", type = "integer", default = 600,
              help = "Plot DPI [default = 600]"),

  make_option("--base_size", type = "double", default = 11,
              help = "Base plot text size [default = 11]"),

  make_option("--axis_text_size", type = "double", default = 8,
              help = "Axis text size [default = 8]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$input_meth))   stop("--input_meth is required")
if (is.null(opt$sample_info))  stop("--sample_info is required")

if (!dir.exists(opt$project_root))  stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$input_meth))   stop("Input methylation file not found: ", opt$input_meth)
if (!file.exists(opt$sample_info))  stop("Sample info file not found: ", opt$sample_info)

if (!is.null(opt$protected_traits_file) && !file.exists(opt$protected_traits_file)) {
  stop("Protected traits file not found: ", opt$protected_traits_file)
}
if (!is.null(opt$technical_traits_file) && !file.exists(opt$technical_traits_file)) {
  stop("Technical traits file not found: ", opt$technical_traits_file)
}
if (opt$threads < 1)      stop("--threads must be >= 1")
if (opt$top_trait_n < 1)  stop("--top_trait_n must be >= 1")
if (opt$top_cell_n < 1)   stop("--top_cell_n must be >= 1")
if (!tolower(opt$primary_cor_method) %in% c("bicor", "pearson")) {
  stop("--primary_cor_method must be 'bicor' or 'pearson'")
}
if (!is.null(opt$selected_trait_mode) &&
    !tolower(opt$selected_trait_mode) %in% c("protected", "protected_and_technical", "none")) {
  stop("--selected_trait_mode must be one of: protected, protected_and_technical, none")
}

# ------------------------------------------------------------
# Configure threads
# ------------------------------------------------------------
options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------------
# Derive output directory
# ------------------------------------------------------------
input_filename <- basename(opt$input_meth)
input_dir      <- dirname(opt$input_meth)

if (grepl("05b_shared_complete_regions", input_dir)) {
  dataset_label <- sub("_Methylation_complete\\.rds$", "", input_filename)
  dataset_label <- sub("\\.rds$", "", dataset_label)
  cpg_label     <- "05b_shared_complete"
  region_label  <- "shared"
} else {
  region_label  <- basename(input_dir)
  cpg_label     <- basename(dirname(input_dir))
  dataset_label <- basename(dirname(dirname(input_dir)))
}

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "06_pc_diagnostics")
out_dir       <- file.path(step_dir, dataset_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Dataset label: ", dataset_label)
message("Output directory: ", out_dir)

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------
pcs_rds                  <- file.path(out_dir, "PCs.rds")
pc_meta_tsv              <- file.path(out_dir, "PC_Metadata.tsv")
var_plot_pdf             <- file.path(out_dir, "Variance_Explained_First20_PCs.pdf")
pc_trait_pearson_tsv     <- file.path(out_dir, "PC_Trait_Correlation_Stats_Pearson.tsv")
pc_trait_bicor_tsv       <- file.path(out_dir, "PC_Trait_Correlation_Stats_Bicor.tsv")
heatmap_pub_pdf          <- file.path(out_dir, "PC_Trait_Bicor_Heatmap_PUBLICATION.pdf")
heatmap_topcells_pdf     <- file.path(out_dir, paste0("PC_Trait_Bicor_Heatmap_Top", opt$top_cell_n, "Cells.pdf"))
trait_dendro_pdf         <- file.path(out_dir, "Trait_Dendrogram_From_PC_Associations.pdf")
removed_non_numeric_file <- file.path(out_dir, "sample_info_non_numeric_columns_removed.txt")
removed_all_na_file      <- file.path(out_dir, "sample_info_all_na_numeric_columns_removed.txt")
removed_zero_var_file    <- file.path(out_dir, "sample_info_zero_variance_numeric_columns_removed.txt")

# ------------------------------------------------------------
# Helper: write vector safely
# ------------------------------------------------------------
write_vector_file <- function(x, file) {
  x <- unique(as.character(x))
  if (length(x) == 0) x <- character(0)
  writeLines(x, con = file)
}

# ------------------------------------------------------------
# Helper: derive PCs from methylation
# ------------------------------------------------------------
getPCs_generic <- function(meth,
                           mod = matrix(1, nrow = ncol(meth), ncol = 1),
                           save = TRUE,
                           file = "PCs.rds",
                           verbose = TRUE) {
  if (verbose) message("[getPCs_generic] Determining number of PCs using sva::num.sv")

  mod_noNA <- mod

  if (ncol(mod_noNA) > 1) {
    intercept_col <- apply(mod_noNA, 2, function(x) length(unique(x)) == 1)
    trait_cols    <- mod_noNA[, !intercept_col, drop = FALSE]

    col_var    <- apply(trait_cols, 2, function(x) var(x, na.rm = TRUE))
    trait_cols <- trait_cols[, !is.na(col_var) & col_var > 0, drop = FALSE]

    if (ncol(trait_cols) > 0) {
      trait_cols <- scale(trait_cols, center = TRUE, scale = TRUE)
      trait_cols[is.na(trait_cols) | is.infinite(trait_cols)] <- 0
      mod_noNA <- cbind(mod_noNA[, intercept_col, drop = FALSE], trait_cols)
    } else {
      mod_noNA <- mod_noNA[, intercept_col, drop = FALSE]
    }
  }

  mod_noNA[is.na(mod_noNA) | is.infinite(mod_noNA)] <- 0

  if (any(!is.finite(mod_noNA))) {
    warning("[getPCs_generic] Non-finite values remain in model matrix; falling back to intercept-only")
    mod_noNA <- matrix(1, nrow = ncol(meth), ncol = 1)
  }

  n_pc <- sva::num.sv(meth, mod = mod_noNA, method = "be", seed = 5)

  n_samples  <- ncol(meth)
  n_features <- nrow(meth)

  if (verbose) {
    message("------------------------------------------------------------")
    message("[SVA / PCA DIAGNOSTICS]")
    message("------------------------------------------------------------")
    message("Input matrix: ", n_features, " features (regions) x ", n_samples, " samples")

    if (ncol(mod_noNA) > 1) {
      message("Model matrix: includes protected variables (signal preserved during SV estimation)")
    } else {
      message("Model matrix: intercept-only (all variation treated as unwanted)")
    }

    message("")
    message("num.sv() result:")
    message("  \u2192 Estimated latent factors (SV-like PCs): ", n_pc)
    message("")
    message("Interpretation:")
    message("  \u2192 These PCs represent major axes of variation NOT explained by the model")
    message("  \u2192 In a protected model, they primarily capture technical/unwanted variation")
    message("  \u2192 They are used as adjustment covariates in downstream models")
    message("------------------------------------------------------------")
  }

  meth_t          <- t(meth)
  meth_t_centered <- scale(meth_t, center = TRUE, scale = FALSE)
  ss              <- svd(meth_t_centered)

  PCs <- ss$u[, seq_len(n_pc), drop = FALSE]
  dimnames(PCs) <- list(rownames(meth_t), paste0("PC_", seq_len(n_pc)))

  var_explained <- (ss$d^2) / sum(ss$d^2)

  if (save) saveRDS(PCs, file = file)

  list(
    PCs             = PCs,
    singular_values = ss$d,
    var_explained   = var_explained,
    n_pc            = n_pc
  )
}

# ------------------------------------------------------------
# Load methylation
# ------------------------------------------------------------
meth <- readRDS(opt$input_meth)

if (!(is.matrix(meth) || is.data.frame(meth))) {
  stop("input_meth must contain a matrix-like object (matrix/data.frame).")
}
meth <- as.matrix(meth)

if (!is.numeric(meth))          stop("Methylation matrix must be numeric.")
if (nrow(meth) < 1)             stop("Methylation matrix has zero rows.")
if (ncol(meth) < 2)             stop("Methylation matrix must have at least 2 samples (columns).")
if (is.null(colnames(meth)))    stop("Methylation matrix must have sample IDs in colnames.")
if (anyDuplicated(colnames(meth))) {
  dup_ids <- unique(colnames(meth)[duplicated(colnames(meth))])
  stop("Duplicate methylation sample IDs found. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}

n_na <- sum(is.na(meth))
if (n_na > 0) {
  stop(
    "Methylation matrix contains ", n_na, " NA values. ",
    "Run Script 05b first to obtain a complete (NA-free) matrix, ",
    "or use complete.cases() to remove incomplete regions before running this script."
  )
}

message("Loaded methylation matrix: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ------------------------------------------------------------
# Load sample info
# ------------------------------------------------------------
sample_info <- readSampleInfo(
  file          = opt$sample_info,
  sample_id_col = opt$sample_id_col,
  verbose       = TRUE
)

if (nrow(sample_info) < 2)      stop("Sample info must contain at least 2 samples after loading.")
if (is.null(rownames(sample_info))) stop("Sample info must have rownames after ID resolution.")
if (anyDuplicated(rownames(sample_info))) {
  dup_ids <- unique(rownames(sample_info)[duplicated(rownames(sample_info))])
  stop("Duplicate sample IDs found in sample info. Example: ", paste(head(dup_ids, 10), collapse = ", "))
}

# ------------------------------------------------------------
# Align samples
# ------------------------------------------------------------
common_samples <- intersect(colnames(meth), rownames(sample_info))
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between methylation matrix and sample info.")
}
if (length(common_samples) < 5) {
  warning("Only ", length(common_samples), " overlapping samples found between methylation matrix and sample info.")
}

meth        <- meth[, common_samples, drop = FALSE]
sample_info <- sample_info[common_samples, , drop = FALSE]

message("After alignment: ", nrow(meth), " regions x ", ncol(meth), " samples")

# ------------------------------------------------------------
# Keep numeric traits only for PC-trait correlations
# ------------------------------------------------------------
is_numeric_atomic   <- vapply(sample_info, function(x) is.numeric(x) && is.atomic(x), logical(1))
removed_non_numeric <- names(sample_info)[!is_numeric_atomic]
sample_info_num     <- sample_info[, is_numeric_atomic, drop = FALSE]

write_vector_file(removed_non_numeric, removed_non_numeric_file)

if (ncol(sample_info_num) == 0) {
  stop("No numeric trait columns available in sample_info after filtering.")
}

all_na   <- vapply(sample_info_num, function(x) all(is.na(x)), logical(1))
zero_var <- vapply(sample_info_num, function(x) {
  vals <- x[!is.na(x)]
  if (length(vals) <= 1) return(TRUE)
  length(unique(vals)) == 1
}, logical(1))

removed_all_na   <- names(sample_info_num)[all_na]
removed_zero_var <- names(sample_info_num)[zero_var]

sample_info_num <- sample_info_num[, !(all_na | zero_var), drop = FALSE]

write_vector_file(removed_all_na,   removed_all_na_file)
write_vector_file(removed_zero_var, removed_zero_var_file)

if (ncol(sample_info_num) == 0) {
  stop("No usable numeric traits remain after removing all-NA and zero-variance columns.")
}

message("Numeric trait columns retained for diagnostics: ", ncol(sample_info_num))

# ------------------------------------------------------------
# Load and resolve protected / technical traits
# ------------------------------------------------------------
protected_traits_requested <- readTraitFile(opt$protected_traits_file, verbose = TRUE)
technical_traits_requested <- readTraitFile(opt$technical_traits_file, verbose = TRUE)

protected_resolved <- resolveTraits(
  requested_traits = protected_traits_requested,
  available_traits = colnames(sample_info_num),
  label            = "protected_traits",
  verbose          = TRUE
)

technical_resolved <- resolveTraits(
  requested_traits = technical_traits_requested,
  available_traits = colnames(sample_info_num),
  label            = "technical_traits",
  verbose          = TRUE
)

if (!is.null(opt$protected_traits_file)) {
  write_vector_file(protected_resolved$requested, file.path(out_dir, "traits_requested_protected.txt"))
  write_vector_file(protected_resolved$found,     file.path(out_dir, "traits_found_protected.txt"))
  write_vector_file(protected_resolved$missing,   file.path(out_dir, "traits_missing_protected.txt"))
}

if (!is.null(opt$technical_traits_file)) {
  write_vector_file(technical_resolved$requested, file.path(out_dir, "traits_requested_technical.txt"))
  write_vector_file(technical_resolved$found,     file.path(out_dir, "traits_found_technical.txt"))
  write_vector_file(technical_resolved$missing,   file.path(out_dir, "traits_missing_technical.txt"))
}

# ------------------------------------------------------------
# Determine selected trait mode
# ------------------------------------------------------------
selected_trait_mode <- opt$selected_trait_mode
if (is.null(selected_trait_mode)) {
  selected_trait_mode <- if (length(protected_resolved$found) > 0) "protected" else "none"
}
selected_trait_mode <- tolower(selected_trait_mode)

selected_traits <- character(0)
if (selected_trait_mode == "protected") {
  selected_traits <- protected_resolved$found
} else if (selected_trait_mode == "protected_and_technical") {
  selected_traits <- unique(c(protected_resolved$found, technical_resolved$found))
}

# ------------------------------------------------------------
# Build model matrix for num.sv
# ------------------------------------------------------------
old_na_action <- getOption("na.action")
options(na.action = "na.pass")
on.exit(options(na.action = old_na_action), add = TRUE)

if (length(protected_resolved$found) > 0) {
  protected_df <- sample_info_num[, protected_resolved$found, drop = FALSE]
  mod_formula  <- paste("~", paste(colnames(protected_df), collapse = " + "))
  mod          <- model.matrix(as.formula(mod_formula), data = protected_df)
  message("Using protected traits in SVA model matrix: ",
          paste(colnames(protected_df), collapse = ", "))
} else {
  mod_formula  <- "~ 1"
  mod          <- matrix(1, nrow = ncol(meth), ncol = 1)
  rownames(mod) <- colnames(meth)
  colnames(mod) <- "(Intercept)"
  message("No protected traits provided/found; using intercept-only model for num.sv")
}

if (nrow(mod) != ncol(meth)) {
  stop("Model matrix row count (", nrow(mod),
       ") does not match methylation sample count (", ncol(meth), ").")
}

# ------------------------------------------------------------
# Derive PCs
# ------------------------------------------------------------
pc_out <- getPCs_generic(
  meth    = meth,
  mod     = mod,
  save    = TRUE,
  file    = pcs_rds,
  verbose = TRUE
)

PCs             <- pc_out$PCs
var_explained   <- pc_out$var_explained
singular_values <- pc_out$singular_values
n_pc            <- pc_out$n_pc
n_total_pcs     <- length(var_explained)

# ------------------------------------------------------------
# Save PC metadata
# ------------------------------------------------------------
pc_meta <- data.frame(
  PC                 = paste0("PC_", seq_along(singular_values)),
  singular_value     = singular_values,
  variance_explained = var_explained,
  rank               = seq_along(singular_values),
  stringsAsFactors   = FALSE
)

readr::write_tsv(pc_meta, pc_meta_tsv)

# ------------------------------------------------------------
# Variance explained plot (first 20 PCs only)
# ------------------------------------------------------------
n_plot_pcs <- min(20, n_total_pcs)

var_df <- data.frame(
  PC = seq_len(n_plot_pcs),
  variance_explained = var_explained[seq_len(n_plot_pcs)]
)

vline_x <- if (n_pc <= n_plot_pcs) n_pc else NA_real_
label_x <- if (n_pc <= n_plot_pcs) min(n_pc + 0.5, n_plot_pcs) else n_plot_pcs - 2

p_var <- ggplot(var_df, aes(x = PC, y = variance_explained)) +
  geom_line(color = "grey40") +
  geom_point(size = 2, color = "steelblue") +
  {
    if (!is.na(vline_x)) geom_vline(xintercept = vline_x, linetype = "dashed", color = "red", linewidth = 0.8)
  } +
  {
    if (!is.na(vline_x)) annotate(
      "text",
      x = label_x,
      y = max(var_df$variance_explained) * 0.95,
      label = paste0("num.sv = ", n_pc),
      hjust = 0,
      size = 4,
      color = "red"
    )
  } +
  scale_x_continuous(breaks = seq_len(n_plot_pcs)) +
  theme_bw(base_size = opt$base_size) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title = paste0("Variance Explained by First ", n_plot_pcs, " PCs: ", dataset_label),
    x     = "Principal Component",
    y     = "Proportion of Variance Explained"
  )

ggsave(
  filename  = var_plot_pdf,
  plot      = p_var,
  width     = opt$plot_width,
  height    = opt$plot_height,
  dpi       = opt$plot_dpi,
  units     = "in",
  limitsize = FALSE
)

message("Variance explained plot saved for first ", n_plot_pcs, " PCs")
if (n_pc > n_plot_pcs) {
  message("num.sv selected ", n_pc, " PCs, which is beyond the first ", n_plot_pcs, " shown in the plot")
} else {
  message("Vertical line indicates num.sv-selected PCs (n = ", n_pc, ")")
}

# ------------------------------------------------------------
# PC-trait correlations
# ------------------------------------------------------------
PCtraitCor_pearson <- getMEtraitCor(
  PCs,
  colData = sample_info_num,
  corType = "pearson",
  file    = pc_trait_pearson_tsv
)

PCtraitCor_bicor <- getMEtraitCor(
  PCs,
  colData = sample_info_num,
  corType = "bicor",
  file    = pc_trait_bicor_tsv
)

primary_method <- tolower(opt$primary_cor_method)
primary_stats  <- if (primary_method == "bicor") PCtraitCor_bicor else PCtraitCor_pearson

# ------------------------------------------------------------
# Publication heatmap: selected traits
# ------------------------------------------------------------
if (length(selected_traits) > 0 && selected_trait_mode != "none") {
  message("Saving publication heatmap (selected traits): ", heatmap_pub_pdf)
  plotPCTrait(
    MEtraitCor     = primary_stats,
    subsetTraits   = selected_traits,
    label_mode     = "star",
    autoSize       = TRUE,
    cell_in        = 0.22,
    max_width      = 14,
    max_height     = 8,
    force_x_labels = TRUE,
    axis.text.size = opt$axis_text_size,
    base_size      = opt$base_size,
    dpi            = opt$plot_dpi,
    file           = heatmap_pub_pdf,
    verbose        = TRUE
  )
} else {
  message("Skipping publication heatmap: no selected traits under selected_trait_mode = ",
          selected_trait_mode)
}

# ------------------------------------------------------------
# Top cells heatmap
# ------------------------------------------------------------
message("Saving top cells heatmap (top ", opt$top_cell_n, " associations): ", heatmap_topcells_pdf)
plotPCTrait(
  MEtraitCor     = primary_stats,
  topCells       = opt$top_cell_n,
  label_mode     = "star",
  autoSize       = FALSE,
  width          = 10,
  height         = 16,
  force_x_labels = TRUE,
  force_y_labels = TRUE,
  axis.text.size = 9,
  wrapY          = 40,
  base_size      = opt$base_size,
  dpi            = opt$plot_dpi,
  file           = heatmap_topcells_pdf,
  verbose        = TRUE
)

# ------------------------------------------------------------
# Trait dendrogram from PC associations
# Increased width so trait labels are easier to read
# ------------------------------------------------------------
plotTraitDendrogramFromPC(
  pc_trait_stats = primary_stats,
  output_file    = trait_dendro_pdf,
  cor_column     = primary_method,
  p_column       = "p",
  width          = 20,
  height         = 12,
  verbose        = TRUE
)

# ------------------------------------------------------------
# Save run parameters
# ------------------------------------------------------------
run_params <- c(
  paste("project_root:", opt$project_root),
  paste("input_meth:", opt$input_meth),
  paste("sample_info:", opt$sample_info),
  paste("sample_id_col:", ifelse(is.null(opt$sample_id_col), "NULL", opt$sample_id_col)),
  paste("protected_traits_file:", ifelse(is.null(opt$protected_traits_file), "NULL", opt$protected_traits_file)),
  paste("technical_traits_file:", ifelse(is.null(opt$technical_traits_file), "NULL", opt$technical_traits_file)),
  paste("threads:", opt$threads),
  paste("top_trait_n:", opt$top_trait_n),
  paste("top_cell_n:", opt$top_cell_n),
  paste("primary_cor_method:", primary_method),
  paste("selected_trait_mode:", selected_trait_mode),
  paste("plot_width:", opt$plot_width),
  paste("plot_height:", opt$plot_height),
  paste("plot_dpi:", opt$plot_dpi),
  paste("base_size:", opt$base_size),
  paste("axis_text_size:", opt$axis_text_size),
  paste("dataset_label:", dataset_label),
  paste("cpg_label:", cpg_label),
  paste("region_label:", region_label),
  paste("model_formula:", mod_formula),
  paste("n_samples_after_alignment:", ncol(meth)),
  paste("n_regions:", nrow(meth)),
  paste("n_total_pcs_from_svd:", n_total_pcs),
  paste("n_pcs_selected_by_numsv:", n_pc),
  paste("n_numeric_traits_used:", ncol(sample_info_num)),
  paste("n_protected_traits_found:", length(protected_resolved$found)),
  paste("n_technical_traits_found:", length(technical_resolved$found)),
  paste("output_dir:", out_dir),
  paste("date:", as.character(Sys.time()))
)

writeLines(run_params, con = file.path(out_dir, "run_parameters.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(out_dir, "sessionInfo.txt"))

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------
message("\u2713 Consensus Script 06 complete: PC derivation and diagnostics finished")
message("Output folder: ", out_dir)

message("------------------------------------------------------------")
message("[LATENT STRUCTURE SUMMARY]")
message("------------------------------------------------------------")
message("Samples: ", nrow(PCs))
message("Total latent dimensions (full PCA space): ", n_total_pcs)
message("Selected latent factors for adjustment (num.sv): ", n_pc)
message("")
message("Returned adjustment matrix: ", nrow(PCs), " x ", ncol(PCs), " (samples x selected PCs)")
message("")
message("Variance captured by selected PCs: ",
        round(sum(var_explained[seq_len(n_pc)]) * 100, 2), "%")
message("")
message("Interpretation:")
message("  \u2192 num.sv estimates how many hidden sources of variation remain after accounting for protected traits")
message("  \u2192 These selected PCs act as surrogate variables (SV-like) for unwanted variation")
message("  \u2192 They are used for downstream methylation adjustment")
message("------------------------------------------------------------")