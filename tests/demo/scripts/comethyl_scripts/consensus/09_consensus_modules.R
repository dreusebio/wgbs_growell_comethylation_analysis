#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 09: Consensus Module Detection
#
# PURPOSE
#   - Load 2-3 dataset methylation matrices for the same adjustment version
#   - Ensure region sets are shared and ordered identically
#   - Build multiExpr object for WGCNA consensus analysis
#   - Run blockwiseConsensusModules()
#   - Save shared consensus module object
#   - Save dataset-specific eigengenes and module tables
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/09_consensus_modules/<adjustment_version>/
#       shared/
#           Consensus_Modules.rds
#           Consensus_Module_Distribution.csv
#           Consensus_Region_Assignments.tsv
#           Consensus_Dendrogram_with_Colors.pdf
#           run_parameters.txt
#           sessionInfo.txt
#
#       <dataset_label>/
#           <dataset_label>_Consensus_MEs.csv
#           <dataset_label>_Consensus_Module_Eigengene_Dendrogram.pdf
# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Parse args
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--regions_file", type = "character"),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_meth", type = "character"),
  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_meth", type = "character", default = NULL),
  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_meth", type = "character", default = NULL),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--power", type = "integer", default = NULL),
  make_option("--chosen_power_file", type = "character", default = NULL),
  make_option("--consensus_cor", type = "character", default = "bicor"),
  make_option("--network_type", type = "character", default = "signed"),
  make_option("--tom_type", type = "character", default = "signed"),
  make_option("--deep_split", type = "integer", default = 4),
  make_option("--min_module_size", type = "integer", default = 5),
  make_option("--merge_cut_height", type = "double", default = 0.1),
  make_option("--max_block_size", type = "integer", default = 40000),
  make_option("--network_calibration", type = "character", default = "full quantile"),
  make_option("--save_consensus_toms", type = "character", default = "TRUE"),
  make_option("--threads", type = "integer", default = 4),
  make_option("--max_p_outliers", type = "double", default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE")
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) stop(label, " must be matrix-like.")
  x <- as.matrix(x)
  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (nrow(x) < 1) stop(label, " has zero rows.")
  if (ncol(x) < 2) stop(label, " must have at least 2 samples.")
  if (is.null(rownames(x))) stop(label, " must have region IDs in rownames.")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(rownames(x))) stop(label, " has duplicated region IDs.")
  if (anyDuplicated(colnames(x))) stop(label, " has duplicated sample IDs.")
  x
}

read_chosen_power <- function(file) {
  if (is.null(file) || !nzchar(file) || !file.exists(file)) return(NULL)
  x <- readLines(file, warn = FALSE)
  hit <- grep("^chosen_power:", x, value = TRUE)
  if (length(hit) == 0) return(NULL)
  as.integer(trimws(sub("^chosen_power:\\s*", "", hit[1])))
}

safe_hclust_plot <- function(hc, file, main = "", cex = 0.6, h = NULL) {
  grDevices::pdf(file, width = 10, height = 5)
  par(mar = c(0, 5, 1, 1))
  plot(hc, main = main, xlab = "", sub = "", cex = cex)
  if (!is.null(h)) abline(h = h, col = "red")
  dev.off()
}

# ------------------------------------------------------------
# Validate args
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$regions_file)) stop("--regions_file is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth)) stop("--dataset1_meth is required")

if (!dir.exists(opt$project_root)) stop("project_root not found")
if (!file.exists(opt$regions_file)) stop("regions_file not found")
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found")

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth)) {
    stop("Provide both dataset2_label and dataset2_meth")
  }
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found")
}
if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) {
    stop("Provide both dataset3_label and dataset3_meth")
  }
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found")
}

consensus_cor <- tolower(opt$consensus_cor)
if (!consensus_cor %in% c("bicor", "pearson")) stop("--consensus_cor must be bicor or pearson")

save_consensus_toms <- parse_bool(opt$save_consensus_toms, "--save_consensus_toms")

chosen_power <- opt$power
if (is.null(chosen_power)) {
  chosen_power <- read_chosen_power(opt$chosen_power_file)
}
if (is.null(chosen_power) || is.na(chosen_power)) {
  stop("Provide --power or a valid --chosen_power_file")
}

WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------------
# Output dirs
# ------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "09_consensus_modules", opt$adjustment_version)
shared_dir <- file.path(step_dir, "shared")
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load regions
# ------------------------------------------------------------
regions <- read.delim(opt$regions_file, stringsAsFactors = FALSE, check.names = FALSE)
required_region_cols <- c("RegionID", "chr", "start", "end")
missing_region_cols <- setdiff(required_region_cols, colnames(regions))
if (length(missing_region_cols) > 0) {
  stop("regions_file missing required columns: ", paste(missing_region_cols, collapse = ", "))
}

# ------------------------------------------------------------
# Load datasets
# ------------------------------------------------------------
dataset_inputs <- list(
  list(label = opt$dataset1_label, meth_file = opt$dataset1_meth)
)
if (dataset2_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset2_label, meth_file = opt$dataset2_meth)
}
if (dataset3_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(label = opt$dataset3_label, meth_file = opt$dataset3_meth)
}

meth_list <- list()
for (ds in dataset_inputs) {
  m <- validate_meth_matrix(readRDS(ds$meth_file), ds$label)
  meth_list[[ds$label]] <- m
}

# ------------------------------------------------------------
# Ensure shared region set and ordering
# ------------------------------------------------------------
ref_regions <- rownames(meth_list[[dataset_inputs[[1]]$label]])

for (nm in names(meth_list)) {
  if (!identical(rownames(meth_list[[nm]]), ref_regions)) {
    if (!all(ref_regions %in% rownames(meth_list[[nm]]))) {
      stop("Dataset ", nm, " does not contain all reference regions.")
    }
    meth_list[[nm]] <- meth_list[[nm]][ref_regions, , drop = FALSE]
  }
}

regions_use <- regions[match(ref_regions, regions$RegionID), , drop = FALSE]
if (any(is.na(regions_use$RegionID))) {
  stop("Some methylation rownames are not present in regions_file")
}

# ------------------------------------------------------------
# Build multiExpr
# WGCNA expects samples x features
# ------------------------------------------------------------
multiExpr <- lapply(meth_list, function(m) {
  list(data = as.data.frame(t(m)))
})

exprSize <- checkSets(multiExpr)
nSets <- exprSize$nSets

# ------------------------------------------------------------
# Run consensus modules
# ------------------------------------------------------------
consensusMods <- blockwiseConsensusModules(
  multiExpr,
  checkMissingData = FALSE,
  maxBlockSize = opt$max_block_size,
  corType = consensus_cor,
  maxPOutliers = opt$max_p_outliers,
  power = chosen_power,
  networkType = opt$network_type,
  checkPower = FALSE,
  TOMType = opt$tom_type,
  networkCalibration = opt$network_calibration,
  saveConsensusTOMs = save_consensus_toms,
  deepSplit = opt$deep_split,
  minModuleSize = opt$min_module_size,
  mergeCutHeight = opt$merge_cut_height,
  verbose = 5
)

saveRDS(consensusMods, file = file.path(shared_dir, "Consensus_Modules.rds"))

# ------------------------------------------------------------
# Shared outputs
# ------------------------------------------------------------
module_dist <- as.data.frame(sort(table(consensusMods$colors), decreasing = TRUE))
colnames(module_dist) <- c("Module", "Regions")
write.csv(module_dist, file.path(shared_dir, "Consensus_Module_Distribution.csv"), row.names = FALSE)

consensus_regions <- regions_use
consensus_regions$module <- consensusMods$colors
write.table(
  consensus_regions,
  file = file.path(shared_dir, "Consensus_Region_Assignments.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

grDevices::pdf(file.path(shared_dir, "Consensus_Dendrogram_with_Colors.pdf"), width = 10, height = 5)
plotDendroAndColors(
  dendro = consensusMods$dendrograms[[1]],
  colors = consensusMods$colors,
  groupLabels = "Consensus Modules",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  marAll = c(1, 5, 1, 0),
  main = "",
  cex.colorLabels = 1.0
)
dev.off()

# ------------------------------------------------------------
# Dataset-specific outputs
# ------------------------------------------------------------
for (ds in names(consensusMods$multiMEs)) {
  ds_dir <- file.path(step_dir, ds)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  me_df <- consensusMods$multiMEs[[ds]]$data
  write.csv(me_df, file = file.path(ds_dir, paste0(ds, "_Consensus_MEs.csv")), row.names = TRUE)

  hc <- if (consensus_cor == "bicor") {
    hclust(as.dist(1 - bicor(me_df, maxPOutliers = opt$max_p_outliers)), method = "average")
  } else {
    hclust(as.dist(1 - cor(me_df, use = "pairwise.complete.obs")), method = "average")
  }

  safe_hclust_plot(
    hc = hc,
    file = file.path(ds_dir, paste0(ds, "_Consensus_Module_Eigengene_Dendrogram.pdf")),
    h = opt$merge_cut_height
  )
}

# ------------------------------------------------------------
# Logs
# ------------------------------------------------------------
writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("regions_file:", opt$regions_file),
    paste("adjustment_version:", opt$adjustment_version),
    paste("datasets:", paste(names(meth_list), collapse = ", ")),
    paste("consensus_cor:", consensus_cor),
    paste("power:", chosen_power),
    paste("network_type:", opt$network_type),
    paste("tom_type:", opt$tom_type),
    paste("deep_split:", opt$deep_split),
    paste("min_module_size:", opt$min_module_size),
    paste("merge_cut_height:", opt$merge_cut_height),
    paste("max_block_size:", opt$max_block_size),
    paste("network_calibration:", opt$network_calibration),
    paste("save_consensus_toms:", save_consensus_toms),
    paste("threads:", opt$threads),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(shared_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()), con = file.path(shared_dir, "sessionInfo.txt"))

message("✓ Script 09 complete: consensus modules built")