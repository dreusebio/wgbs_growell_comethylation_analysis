#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 10: Consensus Diagnostics
#
# PURPOSE
#   - Load Consensus_Modules.rds
#   - Save consensus eigengene network plots
#   - Save per-dataset ME correlation heatmaps
#   - Save per-dataset sample dendrograms from consensus MEs
# ================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(WGCNA)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

option_list <- list(
  make_option("--project_root", type = "character"),
  make_option("--consensus_modules_rds", type = "character"),
  make_option("--adjustment_version", type = "character", default = "unadjusted"),
  make_option("--consensus_cor", type = "character", default = "bicor"),
  make_option("--max_p_outliers", type = "double", default = 0.1)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$project_root)) stop("--project_root is required")
if (is.null(opt$consensus_modules_rds)) stop("--consensus_modules_rds is required")
if (!file.exists(opt$consensus_modules_rds)) stop("consensus_modules_rds not found")

consensusMods <- readRDS(opt$consensus_modules_rds)

pipeline_root <- file.path(opt$project_root, "comethyl_output", "consensus")
step_dir <- file.path(pipeline_root, "10_consensus_diagnostics", opt$adjustment_version)
shared_dir <- file.path(step_dir, "shared")
dir.create(shared_dir, recursive = TRUE, showWarnings = FALSE)

# Shared eigengene networks
consensusMEs <- consensusOrderMEs(consensusMods$multiMEs)

grDevices::pdf(file.path(shared_dir, "Consensus_Eigengene_Networks.pdf"), width = 8, height = 7)
par(cex = 0.8)
plotEigengeneNetworks(
  consensusMEs,
  setLabels = names(consensusMods$multiMEs),
  plotDendrograms = FALSE,
  marHeatmap = c(3, 3, 2, 1),
  zlimPreservation = c(0.5, 1),
  xLabelsAngle = 90
)
dev.off()

# Per dataset
for (ds in names(consensusMods$multiMEs)) {
  ds_dir <- file.path(step_dir, ds)
  dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)

  MEs <- consensusMods$multiMEs[[ds]]$data

  hc_samples <- hclust(as.dist(1 - cor(t(MEs), use = "pairwise.complete.obs")), method = "average")
  grDevices::pdf(file.path(ds_dir, paste0(ds, "_Consensus_ME_Sample_Dendrogram.pdf")), width = 10, height = 5)
  par(mar = c(0, 5, 1, 1))
  plot(hc_samples, main = "", xlab = "", sub = "", cex = 0.6)
  dev.off()

  me_cor <- if (tolower(opt$consensus_cor) == "bicor") {
    bicor(MEs, maxPOutliers = opt$max_p_outliers)
  } else {
    cor(MEs, use = "pairwise.complete.obs")
  }

  grDevices::pdf(file.path(ds_dir, paste0(ds, "_Consensus_ME_Correlation_Heatmap.pdf")), width = 8, height = 7)
  labeledHeatmap(
    Matrix = me_cor,
    xLabels = colnames(MEs),
    yLabels = colnames(MEs),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = NULL,
    setStdMargins = FALSE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = paste0(ds, " Consensus MEs")
  )
  dev.off()
}

writeLines(
  c(
    paste("project_root:", opt$project_root),
    paste("consensus_modules_rds:", opt$consensus_modules_rds),
    paste("adjustment_version:", opt$adjustment_version),
    paste("consensus_cor:", opt$consensus_cor),
    paste("date:", as.character(Sys.time()))
  ),
  con = file.path(shared_dir, "run_parameters.txt")
)

writeLines(capture.output(sessionInfo()), con = file.path(shared_dir, "sessionInfo.txt"))

message("✓ Script 10 complete: consensus diagnostics finished")
