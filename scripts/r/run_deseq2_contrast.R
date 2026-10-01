#!/usr/bin/env Rscript

# Headless R execution bridge for Bulk Treatment-versus-Matched-Control Differential Expression v1
# Executes official Bioconductor DESeq2 negative-binomial Wald test

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    counts = NULL,
    coldata = NULL,
    sizefactors = NULL,
    design = NULL,
    output = NULL
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--counts" && i < length(args)) {
      params$counts <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--coldata" && i < length(args)) {
      params$coldata <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--sizefactors" && i < length(args)) {
      params$sizefactors <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--design" && i < length(args)) {
      params$design <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      params$output <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

if (is.null(params$counts) || is.null(params$coldata) || is.null(params$sizefactors) || is.null(params$design) || is.null(params$output)) {
  cat("Usage: Rscript run_deseq2_contrast.R --counts <counts.tsv> --coldata <coldata.tsv> --sizefactors <sizefactors.tsv> --design <formula> --output <output.json>\n", file = stderr())
  quit(status = 1)
}

timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
r_ver <- as.character(getRversion())
os_info <- as.character(Sys.info()["sysname"])

# Version resolution
deseq2_available <- suppressWarnings(requireNamespace("DESeq2", quietly = TRUE))
deseq2_ver <- if (deseq2_available) as.character(packageVersion("DESeq2")) else "unavailable"

bioc_ver <- "unavailable"
if (suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
  bioc_ver <- as.character(BiocManager::version())
} else if (suppressWarnings(requireNamespace("BiocVersion", quietly = TRUE))) {
  bv <- as.character(packageVersion("BiocVersion"))
  bioc_ver <- sub("^([0-9]+\\.[0-9]+).*", "\\1", bv)
}

write_json <- function(obj, out_path) {
  json_str <- jsonlite::toJSON(obj, auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 8)
  cat(json_str, file = out_path)
}

if (!deseq2_available) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "R package 'DESeq2' is not available in the active environment."
  ), params$output)
  quit(status = 2)
}

# 1. Read count matrix
counts_mat <- tryCatch({
  df <- read.delim(params$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "integer"
  mat
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to read count matrix:", e$message)
  ), params$output)
  quit(status = 3)
})

# 2. Read coldata
coldata_df <- tryCatch({
  read.delim(params$coldata, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to read coldata:", e$message)
  ), params$output)
  quit(status = 3)
})

# 3. Read size factors
sf_df <- tryCatch({
  read.delim(params$sizefactors, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to read size factors:", e$message)
  ), params$output)
  quit(status = 3)
})

# Align sample ordering
sample_names <- colnames(counts_mat)
if (!all(sample_names %in% rownames(coldata_df))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Sample IDs in count matrix do not match coldata rows."
  ), params$output)
  quit(status = 4)
}
coldata_df <- coldata_df[sample_names, , drop = FALSE]

if (!all(sample_names %in% sf_df$sample_id)) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Sample IDs in count matrix do not match size factors table."
  ), params$output)
  quit(status = 4)
}

sf_named <- setNames(as.numeric(sf_df$size_factor), as.character(sf_df$sample_id))
sf_vector <- sf_named[sample_names]

if (any(is.na(sf_vector)) || any(!is.finite(sf_vector)) || any(sf_vector <= 0)) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Invalid size factors: must be non-null, finite, and strictly positive."
  ), params$output)
  quit(status = 4)
}

# Factor configuration: condition MUST have 'control' as reference level
if (!("condition" %in% colnames(coldata_df))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "coldata missing required column 'condition'."
  ), params$output)
  quit(status = 4)
}

coldata_df$condition <- factor(coldata_df$condition, levels = c("control", "treatment"))
if (any(is.na(coldata_df$condition))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "condition values must be either 'control' or 'treatment'."
  ), params$output)
  quit(status = 4)
}

if ("biological_source" %in% colnames(coldata_df)) {
  coldata_df$biological_source <- factor(coldata_df$biological_source)
}

# Check design formula and model matrix
formula_obj <- as.formula(params$design)
design_mat <- tryCatch({
  model.matrix(formula_obj, data = coldata_df)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to build model matrix for design:", e$message)
  ), params$output)
  quit(status = 5)
})

rank <- qr(design_mat)$rank
resid_df <- nrow(coldata_df) - rank

if (rank < ncol(design_mat)) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    design_matrix_rank = rank,
    residual_degrees_of_freedom = resid_df,
    error_message = sprintf("Model matrix is rank deficient: rank %d < columns %d.", rank, ncol(design_mat))
  ), params$output)
  quit(status = 5)
}

if (resid_df <= 0) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    design_matrix_rank = rank,
    residual_degrees_of_freedom = resid_df,
    error_message = sprintf("Residual degrees of freedom must be > 0, got %d.", resid_df)
  ), params$output)
  quit(status = 5)
}

# Construct DESeqDataSet
dds <- tryCatch({
  DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata_df, design = formula_obj)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("DESeqDataSetFromMatrix failed:", e$message)
  ), params$output)
  quit(status = 6)
})

# Assign locked size factors
DESeq2::sizeFactors(dds) <- sf_vector

# Fit DESeq2 model
# Wald test, betaPrior=FALSE, minReplicatesForReplace=Inf (never silently replace outliers)
fit_res <- tryCatch({
  d <- DESeq2::DESeq(
    dds,
    test = "Wald",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf,
    quiet = TRUE
  )
  list(success = TRUE, dds = d)
}, error = function(e) {
  list(success = FALSE, message = e$message)
})

if (!fit_res$success) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    design_matrix_rank = rank,
    residual_degrees_of_freedom = resid_df,
    error_message = paste("DESeq2 model fitting failed:", fit_res$message)
  ), params$output)
  quit(status = 7)
}

# Extract results with explicit contrast: treatment vs control
# Positive log2FC = Treatment > Control
# Negative log2FC = Treatment < Control
# independentFiltering = FALSE to retain complete evaluable gene universe
res <- tryCatch({
  DESeq2::results(
    fit_res$dds,
    contrast = c("condition", "treatment", "control"),
    independentFiltering = FALSE
  )
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    deseq2_version = deseq2_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    design_matrix_rank = rank,
    residual_degrees_of_freedom = resid_df,
    error_message = paste("DESeq2 results extraction failed:", e$message)
  ), params$output)
  quit(status = 8)
})

res_df <- as.data.frame(res)
gene_ids <- rownames(res_df)

gene_results <- vector("list", length(gene_ids))
for (i in seq_along(gene_ids)) {
  gid <- gene_ids[i]
  bm <- res_df$baseMean[i]
  lfc <- res_df$log2FoldChange[i]
  se <- res_df$lfcSE[i]
  st <- res_df$stat[i]
  pv <- res_df$pvalue[i]
  padj <- res_df$padj[i]

  # Determine status and notes
  if (is.na(bm) || bm == 0) {
    status_str <- "ZERO_COUNTS"
    notes_str <- "All samples have zero counts; statistics undefined."
  } else if (is.na(pv)) {
    status_str <- "OUTLIER"
    notes_str <- "p-value is NA due to Cook's distance outlier detection."
  } else {
    status_str <- "OK"
    notes_str <- NA_character_
  }

  gene_results[[i]] <- list(
    canonical_gene_id = gid,
    base_mean = if (is.na(bm)) 0.0 else as.numeric(bm),
    log2_fold_change = if (is.na(lfc)) NA_real_ else as.numeric(lfc),
    lfc_standard_error = if (is.na(se)) NA_real_ else as.numeric(se),
    wald_statistic = if (is.na(st)) NA_real_ else as.numeric(st),
    p_value = if (is.na(pv)) NA_real_ else as.numeric(pv),
    adjusted_p_value_bh = if (is.na(padj)) NA_real_ else as.numeric(padj),
    result_status = status_str,
    result_notes = notes_str
  )
}

out_data <- list(
  status = "SUCCESS",
  design_formula = params$design,
  design_matrix_rank = rank,
  residual_degrees_of_freedom = resid_df,
  sample_ids = sample_names,
  size_factors_used = as.list(sf_vector),
  r_version = r_ver,
  bioc_version = bioc_ver,
  deseq2_version = deseq2_ver,
  operating_system = os_info,
  execution_timestamp = timestamp,
  gene_results = gene_results,
  error_message = NA_character_
)

write_json(out_data, params$output)
quit(status = 0)
