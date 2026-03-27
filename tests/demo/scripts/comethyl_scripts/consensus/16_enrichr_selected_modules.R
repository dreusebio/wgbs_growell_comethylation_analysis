#!/usr/bin/env Rscript

# ================================================================
# SCRIPT 16: Consensus Module Enrichment
#
# PURPOSE
#   - Load annotated regions from script 15
#   - Determine which modules to enrich via:
#       (a) --module_list_file  (explicit list), OR
#       (b) --stats_file from one or more datasets (12A) +
#           --trait_list_file + --p_thresh
#           (selects modules significant in ANY dataset)
#   - Run KEGG and/or Enrichr per selected module
#   - Write per-module enrichment tables and plots
#   - Write summary tables and run logs
#
# NOTE
#   Module assignments are shared across datasets — enrichment
#   runs once per adjustment version. When using stats-based
#   selection, modules significant in ANY dataset are included.
#
# OUTPUT STRUCTURE
#   comethyl_output/consensus/16_module_enrichment/<adjustment_version>/
#       selected_modules.txt
#       selected_modules_source.tsv   (which dataset each came from)
#       selected_modules_summary.tsv
#       <MODULE>_KEGG.tsv
#       <MODULE>_KEGG.xlsx
#       <MODULE>_KEGG_dotplot.pdf
#       <MODULE>_enrichr.xlsx
#       run_parameters.txt
#       run_log.txt
#       sessionInfo.txt
# ================================================================
message("Starting Script 16 ✓")

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# ============================================================
# Helpers
# ============================================================
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[idx + 1])
  default
}

trim_or_null <- function(x) {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x)) return(NULL)
  x
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",")[[1]])
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(default)
  tolower(trimws(x)) %in% c("true", "t", "1", "yes", "y")
}

stop_if_missing <- function(x, label) {
  if (is.null(x) || !nzchar(x))
    stop("Missing required argument: ", label, call. = FALSE)
}

validate_file_exists <- function(path, label) {
  if (!file.exists(path))
    stop(label, " not found: ", path, call. = FALSE)
}

validate_dir_exists <- function(path, label) {
  if (!dir.exists(path))
    stop(label, " not found: ", path, call. = FALSE)
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_log <- function(logfile, ...) {
  txt <- paste0("[", timestamp_now(), "] ", paste0(..., collapse = ""))
  message(txt)
  cat(txt, "\n", file = logfile, append = TRUE)
}

write_lines_safe <- function(x, file) {
  writeLines(as.character(x), con = file, useBytes = TRUE)
}

# ============================================================
# Input readers
# ============================================================
read_list_file <- function(path) {
  x <- readLines(path, warn = FALSE)
  unique(trimws(x[nzchar(trimws(x))]))
}

load_annotation_files <- function(annotation_dir) {
  annotated_file <- file.path(annotation_dir, "Annotated_Regions.tsv")
  gene_list_file <- file.path(annotation_dir, "Module_Gene_List.tsv")
  summary_file   <- file.path(annotation_dir, "Module_Gene_Summary.tsv")

  validate_file_exists(annotated_file, "Annotated_Regions.tsv")
  validate_file_exists(gene_list_file, "Module_Gene_List.tsv")
  validate_file_exists(summary_file,   "Module_Gene_Summary.tsv")

  annotated <- read.delim(annotated_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
  gene_list  <- read.delim(gene_list_file,  stringsAsFactors = FALSE,
                            check.names = FALSE)
  summary    <- read.delim(summary_file,    stringsAsFactors = FALSE,
                            check.names = FALSE)

  # Validate
  req_a <- c("RegionID", "chr", "start", "end", "module")
  miss_a <- setdiff(req_a, colnames(annotated))
  if (length(miss_a) > 0)
    stop("Annotated_Regions.tsv missing columns: ",
         paste(miss_a, collapse = ", "), call. = FALSE)

  req_g <- c("module", "gene_symbol")
  miss_g <- setdiff(req_g, colnames(gene_list))
  if (length(miss_g) > 0)
    stop("Module_Gene_List.tsv missing columns: ",
         paste(miss_g, collapse = ", "), call. = FALSE)

  list(annotated = annotated, gene_list = gene_list, summary = summary)
}

# ============================================================
# Stats-based module selection
# ============================================================
standardize_stats_columns <- function(df, label) {
  cn <- tolower(colnames(df))
  for (col in c("module", "trait")) {
    if (!col %in% cn)
      stop(label, " must contain '", col, "' column.", call. = FALSE)
    names(df)[match(col, cn)] <- col
  }
  p_candidates <- c("p", "pvalue", "p.value", "p_value")
  p_hit <- intersect(p_candidates, cn)
  if (length(p_hit) == 0)
    stop(label, " must contain a p-value column (p, pvalue, p.value, or p_value).",
         call. = FALSE)
  names(df)[match(p_hit[1], cn)] <- "p"
  df
}

read_stats_file <- function(path, label) {
  ext <- tolower(tools::file_ext(path))
  df <- if (ext %in% c("xlsx", "xls")) {
    openxlsx::read.xlsx(path)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext == "csv") {
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported stats file format: ", path, call. = FALSE)
  }
  standardize_stats_columns(df, label)
}

select_modules_from_stats <- function(stats_df, trait_list, p_thresh) {
  stats_df %>%
    dplyr::mutate(p = as.numeric(p)) %>%
    dplyr::filter(!is.na(p), trait %in% trait_list, p <= p_thresh) %>%
    dplyr::pull(module) %>%
    unique() %>%
    as.character()
}

# ============================================================
# Gene list builder
# ============================================================
build_gene_vector <- function(gene_list_df, module_name) {
  gene_list_df %>%
    dplyr::filter(module == module_name) %>%
    dplyr::pull(gene_symbol) %>%
    as.character() %>%
    trimws() %>%
    unique() %>%
    .[!is.na(.) & . != ""]
}

# ============================================================
# KEGG enrichment
# ============================================================
run_kegg <- function(gene_symbols, out_prefix, organism = "hsa") {
  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]
  if (length(gene_symbols) == 0) {
    message("  KEGG: no genes — skipping")
    return(invisible(NULL))
  }

  gene_df <- suppressMessages(
    clusterProfiler::bitr(gene_symbols,
                          fromType = "SYMBOL",
                          toType   = "ENTREZID",
                          OrgDb    = org.Hs.eg.db)
  )
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    message("  KEGG: no Entrez IDs mapped — skipping")
    return(invisible(NULL))
  }

  kegg_res <- suppressMessages(
    clusterProfiler::enrichKEGG(
      gene         = unique(gene_df$ENTREZID),
      organism     = organism,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  )
  kegg_df <- as.data.frame(kegg_res)
  if (is.null(kegg_df) || nrow(kegg_df) == 0) {
    message("  KEGG: no results returned")
    return(invisible(NULL))
  }

  write.table(kegg_df,
              paste0(out_prefix, "_KEGG.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  openxlsx::write.xlsx(kegg_df,
                        paste0(out_prefix, "_KEGG.xlsx"),
                        rowNames = FALSE)

  if (all(c("Description", "p.adjust") %in% colnames(kegg_df))) {
    top <- kegg_df %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = min(25, nrow(kegg_df)))

    if (nrow(top) > 0) {
      top$Description <- factor(top$Description,
                                 levels = rev(top$Description))
      p <- ggplot(top, aes(x = -log10(p.adjust), y = Description)) +
        geom_point() +
        theme_bw(base_size = 12) +
        labs(title = "Top KEGG pathways",
             x = "-log10(adj p)", y = NULL)
      ggsave(paste0(out_prefix, "_KEGG_dotplot.pdf"),
             plot = p, width = 9, height = 6)
    }
  }

  invisible(kegg_df)
}

# ============================================================
# Enrichr
# ============================================================
run_enrichr <- function(gene_symbols, out_prefix, dbs) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    warning("enrichR not installed; skipping Enrichr.")
    return(invisible(NULL))
  }
  gene_symbols <- unique(na.omit(as.character(gene_symbols)))
  gene_symbols <- gene_symbols[nzchar(gene_symbols)]
  if (length(gene_symbols) == 0 || length(dbs) == 0)
    return(invisible(NULL))

  enrichR::setEnrichrSite("Enrichr")
  res <- tryCatch(
    enrichR::enrichr(gene_symbols, dbs),
    error = function(e) {
      warning("Enrichr failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(res) || length(res) == 0) return(invisible(NULL))

  wb <- openxlsx::createWorkbook()
  wrote_any <- FALSE
  for (nm in names(res)) {
    df <- res[[nm]]
    if (is.null(df) || nrow(df) == 0) next
    sheet <- substr(gsub("[^A-Za-z0-9]", "_", nm), 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, df)
    wrote_any <- TRUE
  }
  if (wrote_any)
    openxlsx::saveWorkbook(wb, paste0(out_prefix, "_enrichr.xlsx"),
                            overwrite = TRUE)

  invisible(res)
}

# ============================================================
# Read arguments
# ============================================================
project_root       <- trim_or_null(get_arg("--project_root"))
annotation_dir     <- trim_or_null(get_arg("--annotation_dir"))
adjustment_version <- trim_or_null(get_arg("--adjustment_version", "unadjusted"))

module_list_file   <- trim_or_null(get_arg("--module_list_file"))

# Per-dataset stats files from 12A for stats-based selection
# Modules significant in ANY dataset will be included
dataset1_label      <- trim_or_null(get_arg("--dataset1_label"))
dataset1_stats_file <- trim_or_null(get_arg("--dataset1_stats_file"))
dataset2_label      <- trim_or_null(get_arg("--dataset2_label"))
dataset2_stats_file <- trim_or_null(get_arg("--dataset2_stats_file"))
dataset3_label      <- trim_or_null(get_arg("--dataset3_label"))
dataset3_stats_file <- trim_or_null(get_arg("--dataset3_stats_file"))

trait_list_file    <- trim_or_null(get_arg("--trait_list_file"))
p_thresh           <- as.numeric(get_arg("--p_thresh", "0.05"))
do_kegg            <- to_bool(get_arg("--do_kegg",    "true"),  default = TRUE)
do_enrichr         <- to_bool(get_arg("--do_enrichr", "false"), default = FALSE)
enrichr_dbs        <- split_csv(get_arg("--enrichr_dbs",
  "GO_Biological_Process_2025,GO_Cellular_Component_2025,GO_Molecular_Function_2025,KEGG_2021_Human"))
organism           <- trim_or_null(get_arg("--organism", "hsa"))

# ============================================================
# Validate
# ============================================================
stop_if_missing(project_root,   "--project_root")
stop_if_missing(annotation_dir, "--annotation_dir")

validate_dir_exists(project_root,   "project_root")
validate_dir_exists(annotation_dir, "annotation_dir")

# Module selection validation
if (!is.null(module_list_file)) {
  validate_file_exists(module_list_file, "--module_list_file")
} else {
  if (is.null(trait_list_file))
    stop("Provide either --module_list_file OR --trait_list_file with at least one stats file.",
         call. = FALSE)
  validate_file_exists(trait_list_file, "--trait_list_file")

  # At least one stats file must be provided
  any_stats <- !is.null(dataset1_stats_file) ||
               !is.null(dataset2_stats_file) ||
               !is.null(dataset3_stats_file)
  if (!any_stats)
    stop("Provide at least one --dataset*_stats_file when using --trait_list_file.",
         call. = FALSE)

  if (!is.null(dataset1_stats_file))
    validate_file_exists(dataset1_stats_file, "--dataset1_stats_file")
  if (!is.null(dataset2_stats_file))
    validate_file_exists(dataset2_stats_file, "--dataset2_stats_file")
  if (!is.null(dataset3_stats_file))
    validate_file_exists(dataset3_stats_file, "--dataset3_stats_file")
}

if (do_enrichr && length(enrichr_dbs) == 0)
  stop("do_enrichr=TRUE but no --enrichr_dbs provided.", call. = FALSE)

# ============================================================
# Output directory
# ============================================================
pipeline_root <- file.path(project_root, "comethyl_output", "consensus")
step_dir      <- file.path(pipeline_root, "16_module_enrichment", adjustment_version)
safe_dir_create(step_dir)
message("Output directory: ", step_dir)

log_file    <- file.path(step_dir, "run_log.txt")
params_file <- file.path(step_dir, "run_parameters.txt")

# ============================================================
# Load annotation files from script 15
# ============================================================
append_log(log_file, "Loading annotation files from: ", annotation_dir)
anno_files   <- load_annotation_files(annotation_dir)
gene_list_df <- anno_files$gene_list
module_summary <- anno_files$summary
available_modules <- unique(as.character(gene_list_df$module))
append_log(log_file, "Available annotated modules: ",
           paste(sort(available_modules), collapse = ", "))

# ============================================================
# Module selection
# ============================================================
selected_modules <- character(0)
selection_mode   <- NULL

# Track which dataset each module came from (for source table)
module_source_rows <- list()

if (!is.null(module_list_file)) {
  # (a) Explicit list
  selected_modules <- read_list_file(module_list_file)
  selection_mode   <- "module_list_file"
  append_log(log_file, "Module selection: explicit list (",
             length(selected_modules), " modules)")

} else {
  # (b) Stats-based — union across all provided datasets
  selection_mode <- "stats_based_union_across_datasets"
  trait_list     <- read_list_file(trait_list_file)
  append_log(log_file, "Trait list: ", paste(trait_list, collapse = ", "))

  dataset_stats <- list()
  if (!is.null(dataset1_stats_file))
    dataset_stats[[dataset1_label]] <- dataset1_stats_file
  if (!is.null(dataset2_stats_file))
    dataset_stats[[dataset2_label]] <- dataset2_stats_file
  if (!is.null(dataset3_stats_file))
    dataset_stats[[dataset3_label]] <- dataset3_stats_file

  for (ds_label in names(dataset_stats)) {
    sf       <- dataset_stats[[ds_label]]
    stats_df <- read_stats_file(sf, paste0(ds_label, " stats_file"))
    mods     <- select_modules_from_stats(stats_df,
                                          trait_list = trait_list,
                                          p_thresh   = p_thresh)
    append_log(log_file, ds_label, ": ", length(mods),
               " modules selected (p <= ", p_thresh, ")")

    for (m in mods) {
      module_source_rows[[length(module_source_rows) + 1]] <-
        data.frame(module = m, dataset = ds_label,
                   stringsAsFactors = FALSE)
    }
    selected_modules <- union(selected_modules, mods)
  }
}

# Filter to modules that actually have annotation
selected_modules <- unique(as.character(selected_modules))
selected_modules <- selected_modules[nzchar(selected_modules)]
not_found        <- setdiff(selected_modules, available_modules)
selected_modules <- intersect(selected_modules, available_modules)

if (length(not_found) > 0)
  append_log(log_file, "WARNING: modules not found in annotation, skipped: ",
             paste(not_found, collapse = ", "))

append_log(log_file, "Final selected modules (", length(selected_modules), "): ",
           paste(sort(selected_modules), collapse = ", "))

# Save selection outputs
write_lines_safe(selected_modules,
                 file.path(step_dir, "selected_modules.txt"))
message("Saved: selected_modules.txt")

# Module source table (only relevant for stats-based selection)
if (length(module_source_rows) > 0) {
  source_df <- dplyr::bind_rows(module_source_rows) %>%
    dplyr::filter(module %in% selected_modules) %>%
    dplyr::arrange(module, dataset)
  write.table(source_df,
              file.path(step_dir, "selected_modules_source.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved: selected_modules_source.tsv")
}

# Selected module summary
sel_summary <- module_summary %>%
  dplyr::filter(module %in% selected_modules)
write.table(sel_summary,
            file.path(step_dir, "selected_modules_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: selected_modules_summary.tsv")

if (length(selected_modules) == 0) {
  append_log(log_file, "No modules selected — no enrichment to run.")
  message("WARNING: No modules selected after filtering. Exiting.")
} else {

  # ============================================================
  # Per-module enrichment
  # ============================================================
  for (m in selected_modules) {
    gene_vec   <- build_gene_vector(gene_list_df, m)
    out_prefix <- file.path(step_dir, m)

    append_log(log_file, "Module ", m, ": ", length(gene_vec), " genes")

    if (do_kegg) {
      message("  Running KEGG for module: ", m)
      tryCatch(
        run_kegg(gene_vec, out_prefix, organism = organism),
        error = function(e)
          append_log(log_file, "KEGG failed for ", m, ": ", conditionMessage(e))
      )
    }

    if (do_enrichr) {
      message("  Running Enrichr for module: ", m)
      tryCatch(
        run_enrichr(gene_vec, out_prefix, dbs = enrichr_dbs),
        error = function(e)
          append_log(log_file, "Enrichr failed for ", m, ": ", conditionMessage(e))
      )
    }
  }
}

# ============================================================
# Logs
# ============================================================
write_lines_safe(
  c(
    paste0("timestamp\t",          timestamp_now()),
    paste0("project_root\t",       project_root),
    paste0("annotation_dir\t",     annotation_dir),
    paste0("adjustment_version\t", adjustment_version),
    paste0("selection_mode\t",     selection_mode),
    paste0("module_list_file\t",
           ifelse(is.null(module_list_file), "", module_list_file)),
    paste0("trait_list_file\t",
           ifelse(is.null(trait_list_file), "", trait_list_file)),
    paste0("dataset1_label\t",
           ifelse(is.null(dataset1_label), "", dataset1_label)),
    paste0("dataset1_stats_file\t",
           ifelse(is.null(dataset1_stats_file), "", dataset1_stats_file)),
    paste0("dataset2_label\t",
           ifelse(is.null(dataset2_label), "", dataset2_label)),
    paste0("dataset2_stats_file\t",
           ifelse(is.null(dataset2_stats_file), "", dataset2_stats_file)),
    paste0("dataset3_label\t",
           ifelse(is.null(dataset3_label), "", dataset3_label)),
    paste0("dataset3_stats_file\t",
           ifelse(is.null(dataset3_stats_file), "", dataset3_stats_file)),
    paste0("p_thresh\t",           p_thresh),
    paste0("do_kegg\t",            do_kegg),
    paste0("do_enrichr\t",         do_enrichr),
    paste0("organism\t",           organism),
    paste0("enrichr_dbs\t",        paste(enrichr_dbs, collapse = ",")),
    paste0("n_selected_modules\t", length(selected_modules)),
    paste0("selected_modules\t",   paste(sort(selected_modules), collapse = ", "))
  ),
  params_file
)

writeLines(capture.output(sessionInfo()),
           con = file.path(step_dir, "sessionInfo.txt"))

append_log(log_file, "Script 16 complete")
message("\n✓ Script 16 complete: consensus module enrichment finished")
message("Outputs saved under: ", step_dir)