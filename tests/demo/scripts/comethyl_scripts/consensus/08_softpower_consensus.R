#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 08: Consensus Soft-Threshold Power Analysis
#
# Pipeline: comethyl WGBS consensus analysis
#
# PURPOSE
#   - Load one adjustment version across multiple datasets
#   - Run soft-threshold analysis separately for each dataset
#   - Generate a combined multi-dataset soft-power plot
#   - Automatically choose one shared power across datasets
#
# DESIGN
#   This script is organized by adjustment version, not by dataset.
#   Example:
#     run all v1 matrices together across Baseline / TP36_38weeks / Postpartum
#     run all v2 matrices together only if they exist
#     run all v3 matrices together only if they exist
#
# REQUIRED INPUTS
#   --project_root
#   --dataset1_label
#   --dataset1_meth
#
# OPTIONAL INPUTS
#   --dataset2_label
#   --dataset2_meth
#   --dataset3_label
#   --dataset3_meth
#
#   --adjustment_version   label for this consensus run
#                          [default = unadjusted]
#   --softpower_cor        bicor or pearson [default = pearson]
#   --power_min            minimum power [default = 1]
#   --power_max            maximum power [default = 20]
#   --block_size           [default = 20000]
#   --gc_interval          [default = 19999]
#   --threads              [default = 4]
#   --plot_sample_dendro   TRUE/FALSE [default = TRUE]
#   --scale_free_cutoff    [default = 0.8]
#
# OUTPUTS
#   <project_root>/comethyl_output/consensus/06_soft_power/<adjustment_version>/
#       chosen_power.txt
#       combined_softpower_summary.tsv
#       Consensus_SoftPower_Combined.pdf
#       run_parameters.txt
#       sessionInfo.txt
#
#       <dataset_label>/
#           SoftPower_<cor>.rds
#           SoftPower_<cor>_table.tsv
#           Sample_Dendrogram_<cor>.pdf   (optional)

# Examples

# Rscript scripts/consensus/08_consensus_softpower.R \
#   --project_root /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis \
#   --dataset1_label Baseline \
#   --dataset1_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/05_methylation_adjustment/Baseline/cov3_75pct/covMin4_methSD0p08/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --dataset2_label TP36_38weeks \
#   --dataset2_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/05_methylation_adjustment/TP36_38weeks/cov1_1pct/covMin4_methSD0p08/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --dataset3_label Postpartum \
#   --dataset3_meth /quobyte/lasallegrp/projects/wgbs_growell_dmr_analysis/comethyl_output/consensus/05_methylation_adjustment/Postpartum/cov1_1pct/covMin4_methSD0p08/v1_all_pcs/Adjusted_Region_Methylation_allPCs.rds \
#   --adjustment_version v1_all_pcs \
#   --softpower_cor bicor \
#   --power_min 1 \
#   --power_max 20 \
#   --scale_free_cutoff 0.8
# ================================================================
message("Starting ✓")

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(comethyl)
})

# ------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------
option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--dataset1_label", type = "character"),
  make_option("--dataset1_meth", type = "character"),

  make_option("--dataset2_label", type = "character", default = NULL),
  make_option("--dataset2_meth", type = "character", default = NULL),

  make_option("--dataset3_label", type = "character", default = NULL),
  make_option("--dataset3_meth", type = "character", default = NULL),

  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--softpower_cor", type = "character", default = "pearson"),
  make_option("--power_min", type = "integer", default = 1),
  make_option("--power_max", type = "integer", default = 20),
  make_option("--block_size", type = "integer", default = 20000),
  make_option("--gc_interval", type = "integer", default = 19999),
  make_option("--threads", type = "integer", default = 4),
  make_option("--plot_sample_dendro", type = "character", default = "TRUE"),
  make_option("--scale_free_cutoff", type = "double", default = 0.8)
)

opt <- parse_args(OptionParser(option_list = option_list))

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_bool <- function(x, arg_name) {
  x2 <- tolower(trimws(as.character(x)))
  if (x2 %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x2 %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(arg_name, " must be TRUE or FALSE. Got: ", x)
}

validate_meth_matrix <- function(x, label) {
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop(label, " must be matrix-like.")
  }
  x <- as.matrix(x)
  if (!is.numeric(x)) stop(label, " must be numeric.")
  if (nrow(x) < 1) stop(label, " has zero rows.")
  if (ncol(x) < 2) stop(label, " must have at least 2 samples.")
  if (is.null(colnames(x))) stop(label, " must have sample IDs in colnames.")
  if (anyDuplicated(colnames(x))) {
    dup_ids <- unique(colnames(x)[duplicated(colnames(x))])
    stop(label, " has duplicated sample IDs. Example: ",
         paste(head(dup_ids, 10), collapse = ", "))
  }
  x
}

extract_sft_table <- function(sft) {
  if (!is.list(sft) || is.null(sft$fitIndices)) {
    stop("Soft-power object does not contain fitIndices.")
  }
  as.data.frame(sft$fitIndices)
}

choose_consensus_power <- function(summary_df, cutoff = 0.8) {
  if (nrow(summary_df) == 0) stop("summary_df is empty.")

  power_summary <- summary_df %>%
    group_by(Power) %>%
    summarise(
      n_datasets = n(),
      n_meet_cutoff = sum(SFT_R2 >= cutoff, na.rm = TRUE),
      min_SFT_R2 = min(SFT_R2, na.rm = TRUE),
      median_SFT_R2 = median(SFT_R2, na.rm = TRUE),
      mean_connectivity = mean(Mean_Connectivity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Power)

  # Best case: all datasets meet cutoff
  all_ok <- power_summary %>% filter(n_meet_cutoff == n_datasets)
  if (nrow(all_ok) > 0) {
    chosen <- all_ok$Power[1]
    reason <- "smallest power where all datasets meet scale-free cutoff"
    return(list(power = chosen, reason = reason, summary = power_summary))
  }

  # Fallback: maximize number of datasets meeting cutoff, then smallest power
  fallback <- power_summary %>%
    arrange(desc(n_meet_cutoff), Power)

  chosen <- fallback$Power[1]
  reason <- "fallback: smallest power maximizing number of datasets meeting cutoff"

  list(power = chosen, reason = reason, summary = power_summary)
}

make_combined_softpower_plot <- function(summary_df, out_file, cutoff = 0.8,
                                         chosen_power = NULL) {
  p1 <- ggplot(summary_df, aes(x = Power, y = SFT_R2,
                                color = Dataset, group = Dataset)) +
    geom_line() +
    geom_point() +
    geom_text(aes(label = Power), vjust = -0.7, size = 2.8, show.legend = FALSE) +
    geom_hline(yintercept = cutoff, linetype = "dashed", color = "black") +
    { if (!is.null(chosen_power))
        geom_vline(xintercept = chosen_power, linetype = "dotted",
                   color = "black", linewidth = 0.8) } +
    { if (!is.null(chosen_power))
        annotate("text", x = chosen_power + 0.3,
                 y = 0.05, label = paste0("power = ", chosen_power),
                 hjust = 0, size = 3.5) } +
    scale_x_continuous(breaks = unique(summary_df$Power)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Consensus Soft-Threshold Selection",
      x     = "Soft Threshold (Power)",
      y     = "Scale-Free Topology Model Fit (R\u00b2)"
    )

  p2 <- ggplot(summary_df, aes(x = Power, y = Mean_Connectivity,
                                color = Dataset, group = Dataset)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = unique(summary_df$Power)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Mean Connectivity",
      x     = "Soft Threshold (Power)",
      y     = "Mean Connectivity"
    )

  p3 <- ggplot(summary_df, aes(x = Power, y = Median_Connectivity,
                                color = Dataset, group = Dataset)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = unique(summary_df$Power)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Median Connectivity",
      x     = "Soft Threshold (Power)",
      y     = "Median Connectivity"
    )

  p4 <- ggplot(summary_df, aes(x = Power, y = Max_Connectivity,
                                color = Dataset, group = Dataset)) +
    geom_line() +
    geom_point() +
    scale_x_continuous(breaks = unique(summary_df$Power)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank()) +
    labs(
      title = "Max Connectivity",
      x     = "Soft Threshold (Power)",
      y     = "Max Connectivity"
    )

  grDevices::pdf(out_file, width = 14, height = 10)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
  dev.off()
}

# ------------------------------------------------------------
# Validate arguments
# ------------------------------------------------------------
if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$dataset1_label)) stop("--dataset1_label is required")
if (is.null(opt$dataset1_meth)) stop("--dataset1_meth is required")

if (!dir.exists(opt$project_root)) stop("Project root does not exist: ", opt$project_root)
if (!file.exists(opt$dataset1_meth)) stop("dataset1_meth not found: ", opt$dataset1_meth)

dataset2_provided <- !is.null(opt$dataset2_label) || !is.null(opt$dataset2_meth)
dataset3_provided <- !is.null(opt$dataset3_label) || !is.null(opt$dataset3_meth)

if (dataset2_provided) {
  if (is.null(opt$dataset2_label) || is.null(opt$dataset2_meth)) {
    stop("Provide both --dataset2_label and --dataset2_meth.")
  }
  if (!file.exists(opt$dataset2_meth)) stop("dataset2_meth not found: ", opt$dataset2_meth)
}

if (dataset3_provided) {
  if (is.null(opt$dataset3_label) || is.null(opt$dataset3_meth)) {
    stop("Provide both --dataset3_label and --dataset3_meth.")
  }
  if (!file.exists(opt$dataset3_meth)) stop("dataset3_meth not found: ", opt$dataset3_meth)
}

softpower_cor <- tolower(opt$softpower_cor)
if (!softpower_cor %in% c("bicor", "pearson")) {
  stop("--softpower_cor must be bicor or pearson")
}
if (opt$power_min < 1) stop("--power_min must be >= 1")
if (opt$power_max < opt$power_min) stop("--power_max must be >= power_min")
if (opt$power_max > 20) stop("For consensus soft power, power_max should not exceed 20.")
if (opt$block_size < 1) stop("--block_size must be >= 1")
if (opt$gc_interval < 0) stop("--gc_interval must be >= 0")
if (opt$threads < 1) stop("--threads must be >= 1")
if (opt$scale_free_cutoff <= 0 || opt$scale_free_cutoff > 1) stop("--scale_free_cutoff must be > 0 and <= 1")

plot_sample_dendro <- parse_bool(opt$plot_sample_dendro, "--plot_sample_dendro")
powerVector <- seq(opt$power_min, opt$power_max)

# ------------------------------------------------------------
# Threads
# ------------------------------------------------------------
WGCNA::enableWGCNAThreads(nThreads = opt$threads)

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "08_soft_power")
out_dir <- file.path(step_dir, opt$adjustment_version)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", out_dir)

# ------------------------------------------------------------
# Build dataset input list
# ------------------------------------------------------------
dataset_inputs <- list(
  list(label = opt$dataset1_label, file = opt$dataset1_meth)
)

if (dataset2_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(
    label = opt$dataset2_label,
    file = opt$dataset2_meth
  )
}

if (dataset3_provided) {
  dataset_inputs[[length(dataset_inputs) + 1]] <- list(
    label = opt$dataset3_label,
    file = opt$dataset3_meth
  )
}

# ------------------------------------------------------------
# Run soft power per dataset
# ------------------------------------------------------------
all_summary <- list()

for (ds in dataset_inputs) {
  ds_label <- ds$label
  ds_file <- ds$file

  message("\n==============================")
  message("Running soft power for dataset: ", ds_label)
  message("==============================")

  meth <- validate_meth_matrix(readRDS(ds_file), ds_label)

  ds_out_dir <- file.path(out_dir, ds_label)
  dir.create(ds_out_dir, recursive = TRUE, showWarnings = FALSE)

  sft_rds <- file.path(ds_out_dir, paste0("SoftPower_", softpower_cor, ".rds"))
  sft_table_tsv <- file.path(ds_out_dir, paste0("SoftPower_", softpower_cor, "_table.tsv"))
  dendro_pdf <- file.path(ds_out_dir, paste0("Sample_Dendrogram_", softpower_cor, ".pdf"))

  sft <- getSoftPower(
    meth,
    powerVector = powerVector,
    corType = softpower_cor,
    file = sft_rds,
    blockSize = opt$block_size,
    gcInterval = opt$gc_interval
  )

  sft_table <- extract_sft_table(sft)

  # Standardize columns
  colnames_lower <- tolower(colnames(sft_table))

  power_col <- colnames(sft_table)[match("power", colnames_lower)]
  fit_col <- colnames(sft_table)[match("sft.r.sq", colnames_lower)]
  if (is.na(fit_col)) fit_col <- colnames(sft_table)[match("sft.rsq", colnames_lower)]
  mean_col <- colnames(sft_table)[match("mean.k.", colnames_lower)]
  if (is.na(mean_col)) mean_col <- colnames(sft_table)[match("mean.k", colnames_lower)]
  median_col <- colnames(sft_table)[match("median.k.", colnames_lower)]
  if (is.na(median_col)) median_col <- colnames(sft_table)[match("median.k", colnames_lower)]
  max_col <- colnames(sft_table)[match("max.k.", colnames_lower)]
  if (is.na(max_col)) max_col <- colnames(sft_table)[match("max.k", colnames_lower)]

  if (any(is.na(c(power_col, fit_col, mean_col, median_col, max_col)))) {
    stop("Could not identify expected columns in soft-power table for dataset: ", ds_label)
  }

  summary_df <- tibble::tibble(
    Dataset = ds_label,
    Power = sft_table[[power_col]],
    SFT_R2 = abs(sft_table[[fit_col]]),
    Mean_Connectivity = sft_table[[mean_col]],
    Median_Connectivity = sft_table[[median_col]],
    Max_Connectivity = sft_table[[max_col]]
  )

  readr::write_tsv(summary_df, sft_table_tsv)
  all_summary[[ds_label]] <- summary_df

  if (isTRUE(plot_sample_dendro)) {
    tryCatch({
      d <- getDendro(meth, distance = "euclidean")
      plotDendro(d, file = dendro_pdf, expandY = c(0.25, 0.08))
    }, error = function(e) {
      warning("Sample dendrogram failed for ", ds_label, ": ", conditionMessage(e))
    })
  }
}

combined_summary <- bind_rows(all_summary)
readr::write_tsv(combined_summary, file.path(out_dir, "combined_softpower_summary.tsv"))

# ------------------------------------------------------------
# Choose one shared power
# ------------------------------------------------------------
chosen <- choose_consensus_power(combined_summary, cutoff = opt$scale_free_cutoff)

writeLines(
  c(
    paste("adjustment_version:", opt$adjustment_version),
    paste("scale_free_cutoff:", opt$scale_free_cutoff),
    paste("chosen_power:", chosen$power),
    paste("reason:", chosen$reason)
  ),
  con = file.path(out_dir, "chosen_power.txt")
)

# ------------------------------------------------------------
# Combined plot
# ------------------------------------------------------------
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  warning("gridExtra not installed; skipping combined softpower PDF.")
} else {
  make_combined_softpower_plot(
  summary_df   = combined_summary,
  out_file     = file.path(out_dir, paste0("Consensus_SoftPower_Combined_", softpower_cor, ".pdf")),
  cutoff       = opt$scale_free_cutoff,
  chosen_power = chosen$power
  )
}

# ------------------------------------------------------------
# Save top-level run parameters
# ------------------------------------------------------------
writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("adjustment_version:", opt$adjustment_version),
    paste("dataset1_label:", opt$dataset1_label),
    paste("dataset1_meth:", opt$dataset1_meth),
    paste("dataset2_label:", ifelse(is.null(opt$dataset2_label), "NULL", opt$dataset2_label)),
    paste("dataset2_meth:", ifelse(is.null(opt$dataset2_meth), "NULL", opt$dataset2_meth)),
    paste("dataset3_label:", ifelse(is.null(opt$dataset3_label), "NULL", opt$dataset3_label)),
    paste("dataset3_meth:", ifelse(is.null(opt$dataset3_meth), "NULL", opt$dataset3_meth)),
    paste("softpower_cor:", softpower_cor),
    paste("power_min:", opt$power_min),
    paste("power_max:", opt$power_max),
    paste("power_vector:", paste(powerVector, collapse = ",")),
    paste("block_size:", opt$block_size),
    paste("gc_interval:", opt$gc_interval),
    paste("threads:", opt$threads),
    paste("plot_sample_dendro:", plot_sample_dendro),
    paste("scale_free_cutoff:", opt$scale_free_cutoff),
    paste("chosen_power:", chosen$power),
    paste("chosen_reason:", chosen$reason),
    paste("n_datasets:", length(dataset_inputs)),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(out_dir, "run_parameters.txt")
)

writeLines(
  capture.output(sessionInfo()),
  con = file.path(out_dir, "sessionInfo.txt")
)

message("\nCONSENSUS SOFT-POWER COMPLETE ✓")
message("Chosen shared power: ", chosen$power)
message("Outputs saved under:\n  ", out_dir)