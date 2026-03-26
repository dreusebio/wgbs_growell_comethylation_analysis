#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 11: Consensus Module Membership
#
# PURPOSE
#   - Load Consensus_Modules.rds
#   - Load dataset methylation matrices used in consensus analysis
#   - Compute region-module membership per dataset
#   - Save membership tables for each dataset
# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(openxlsx)
})

options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_meth", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_meth", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_meth", type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--membership_cor", type = "character", default = "bicor"),
  make_option("--max_p_outliers", type = "double", default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) stop(label, " must be matrix-like.")
  x <- as.matrix(x)
  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (is.null(rownames(x))) stop(label, " must have region IDs in rownames.")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  x
}

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$consensus_modules_rds)) stop("--consensus_modules_rds is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth)) stop("--dataset1_meth is required")

if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found")
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found")

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth)) stop("Provide both dataset2_label and dataset2_meth")
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found")
}
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) stop("Provide both dataset3_label and dataset3_meth")
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found")
}

consensusMods <- readRDS(opt$consensus_modules_rds)

dataset_inputs <- list(
  list(label = opt$dataset1_label, meth_file = opt$dataset1_meth)
)
if (dataset2_provided) dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset2_label, meth_file = opt$dataset2_meth)
if (dataset3_provided) dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset3_label, meth_file = opt$dataset3_meth)

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "11_consensus_membership", opt$adjustment_version)
dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

for (ds in dataset_inputs) {
  ds_label <- ds$label
  ds_dir <- file.path(step_dir, ds_label)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  meth <- validate_meth_matrix(readRDS(ds$meth_file), ds_label)

  # WGCNA membership expects samples x features
  expr_df <- as.data.frame(t(meth))
  MEs <- consensusMods$multiMEs[[ds_label]]$data

  if (tolower(opt$membership_cor) == "bicor") {
    mm <- bicorAndPvalue(
      expr_df,
      MEs,
      use = "pairwise.complete.obs",
      maxPOutliers = opt$max_p_outliers
    )
    membership_df <- as.data.frame(mm$bicor)
    p_df <- as.data.frame(mm$p)
  } else {
    mm <- corAndPvalue(expr_df, MEs, use = "pairwise.complete.obs")
    membership_df <- as.data.frame(mm$cor)
    p_df <- as.data.frame(mm$p)
  }

  colnames(membership_df) <- gsub("^ME", "", colnames(membership_df))
  colnames(p_df) <- paste0(gsub("^ME", "", colnames(p_df)), "_p")

  membership_df$RegionID <- rownames(membership_df)
  membership_df$Module <- consensusMods$colors

  out_df <- cbind(membership_df, p_df)

  write.table(
    out_df,
    file = file.path(ds_dir, paste0(ds_label, "_Consensus_Probe_Module_Membership.txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  openxlsx::write.xlsx(
    out_df,
    file = file.path(ds_dir, paste0(ds_label, "_Consensus_Probe_Module_Membership.xlsx")),
    rowNames = FALSE
  )

  writeLines(
    c(
      paste("project_root:", opt$project_root),
      paste("consensus_modules_rds:", opt$consensus_modules_rds),
      paste("dataset_label:", ds_label),
      paste("dataset_meth:", ds$meth_file),
      paste("adjustment_version:", opt$adjustment_version),
      paste("membership_cor:", opt$membership_cor),
      paste("date:", as.character(Sys.time()))
    ),
    con = file.path(ds_dir, "run_parameters.txt")
  )
}

writeLines(capture.output(sessionInfo()), con = file.path(step_dir, "sessionInfo.txt"))

message("✓ Script 11 complete: consensus membership finished")