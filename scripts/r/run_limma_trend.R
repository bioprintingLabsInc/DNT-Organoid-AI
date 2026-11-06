#!/usr/bin/env Rscript

# Headless R execution bridge for Bulk Processed Expression Differential Inference v1
# Executes official Bioconductor limma-trend empirical Bayes linear modeling for continuous data

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    expression = NULL,
    coldata = NULL,
    contrast = NULL,
    design = "~ 0 + condition",
    output = NULL
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--expression" && i < length(args)) {
      params$expression <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--coldata" && i < length(args)) {
      params$coldata <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--contrast" && i < length(args)) {
      params$contrast <- args[i + 1]
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

if (is.null(params$expression) || is.null(params$coldata) || is.null(params$contrast) || is.null(params$output)) {
  cat("Usage: Rscript run_limma_trend.R --expression <expr.tsv> --coldata <coldata.tsv> --contrast <contrast_str> [--design <formula>] --output <out.json>\n", file = stderr())
  quit(status = 1)
}

timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
r_ver <- as.character(getRversion())
os_info <- as.character(Sys.info()["sysname"])

# Version resolution
limma_available <- suppressWarnings(requireNamespace("limma", quietly = TRUE))
limma_ver <- if (limma_available) as.character(packageVersion("limma")) else "unavailable"

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

if (!limma_available) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "R package 'limma' is not available in the active environment."
  ), params$output)
  quit(status = 2)
}

# 1. Read continuous expression matrix
expr_mat <- tryCatch({
  df <- read.delim(params$expression, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  mat
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to read expression matrix:", e$message)
  ), params$output)
  quit(status = 3)
})

# Verify finite values
if (any(!is.finite(expr_mat))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Expression matrix contains non-finite (NaN, NA, Inf) values."
  ), params$output)
  quit(status = 3)
}

# 2. Read coldata
coldata_df <- tryCatch({
  read.delim(params$coldata, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to read coldata:", e$message)
  ), params$output)
  quit(status = 3)
})

# Align sample ordering
sample_names <- colnames(expr_mat)
if (!all(sample_names %in% rownames(coldata_df))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Sample IDs in expression matrix do not match coldata rows."
  ), params$output)
  quit(status = 4)
}
coldata_df <- coldata_df[sample_names, , drop = FALSE]

if (!("condition" %in% colnames(coldata_df))) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Coldata must contain a 'condition' column."
  ), params$output)
  quit(status = 4)
}

# Biological replication verification: each condition in contrast must have >= 2 samples
condition_counts <- table(coldata_df$condition)
if (any(condition_counts < 2)) {
  under_rep <- names(condition_counts)[condition_counts < 2]
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste0("Insufficient biological replication: condition(s) ", paste(under_rep, collapse = ", "), " have < 2 biological replicates.")
  ), params$output)
  quit(status = 5)
}

# Ensure condition is factor
coldata_df$condition <- factor(coldata_df$condition)

# 3. Design matrix construction
design_formula <- tryCatch({
  as.formula(params$design)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Invalid design formula:", e$message)
  ), params$output)
  quit(status = 6)
})

design_mat <- tryCatch({
  model.matrix(design_formula, data = coldata_df)
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("Failed to build model.matrix:", e$message)
  ), params$output)
  quit(status = 6)
})

# Clean colnames of design_mat
if (grepl("0\\s*\\+\\s*condition", params$design) || grepl("condition\\s*-\\s*1", params$design)) {
  colnames(design_mat) <- make.names(levels(coldata_df$condition))
} else {
  colnames(design_mat) <- make.names(colnames(design_mat))
}

# Rank check
qr_d <- qr(design_mat)
rank_val <- qr_d$rank
p_cols <- ncol(design_mat)
n_samples <- nrow(design_mat)

if (rank_val < p_cols) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste0("Design matrix is rank-deficient (rank = ", rank_val, ", columns = ", p_cols, "). Model cannot be uniquely fit.")
  ), params$output)
  quit(status = 7)
}

if (n_samples - rank_val < 1) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = "Zero residual degrees of freedom: cannot estimate variance for empirical Bayes moderation."
  ), params$output)
  quit(status = 7)
}

# 4. Fit linear model and limma-trend
fit_result <- tryCatch({
  fit <- limma::lmFit(expr_mat, design_mat)
  
  # Format contrast string using clean make.names if necessary
  contrast_clean <- params$contrast
  contrast_mat <- limma::makeContrasts(contrasts = contrast_clean, levels = design_mat)
  fit_cont <- limma::contrasts.fit(fit, contrast_mat)
  fit_eb <- limma::eBayes(fit_cont, trend = TRUE)
  top_tab <- limma::topTable(fit_eb, number = Inf, adjust.method = "BH", sort.by = "none")
  
  list(
    fit = fit_eb,
    table = top_tab,
    df_residual = fit_eb$df.residual[1]
  )
}, error = function(e) {
  write_json(list(
    status = "ERROR",
    r_version = r_ver,
    bioc_version = bioc_ver,
    limma_version = limma_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = paste("limma-trend model fitting failed:", e$message)
  ), params$output)
  quit(status = 8)
})

top_df <- fit_result$table
gene_names <- rownames(top_df)

gene_results <- vector("list", nrow(top_df))
for (idx in seq_len(nrow(top_df))) {
  gene_results[[idx]] <- list(
    canonical_gene_id = gene_names[idx],
    base_mean = as.numeric(top_df$AveExpr[idx]),
    log2_fold_change = as.numeric(top_df$logFC[idx]),
    lfc_standard_error = NULL, # limma topTable does not provide Wald standard error; uninvented
    moderated_t = as.numeric(top_df$t[idx]),
    p_value = as.numeric(top_df$P.Value[idx]),
    adjusted_p_value_bh = as.numeric(top_df$adj.P.Val[idx]),
    b_statistic = as.numeric(top_df$B[idx]),
    result_status = "OK"
  )
}

output_payload <- list(
  status = "SUCCESS",
  analysis_method = "limma_trend",
  test_statistic_name = "moderated_t",
  r_version = r_ver,
  bioc_version = bioc_ver,
  limma_version = limma_ver,
  operating_system = os_info,
  execution_timestamp = timestamp,
  design_formula = params$design,
  contrast = params$contrast,
  design_matrix_rank = rank_val,
  residual_degrees_of_freedom = as.numeric(fit_result$df_residual),
  input_gene_count = nrow(expr_mat),
  input_sample_count = ncol(expr_mat),
  gene_results = gene_results
)

write_json(output_payload, params$output)
quit(status = 0)
