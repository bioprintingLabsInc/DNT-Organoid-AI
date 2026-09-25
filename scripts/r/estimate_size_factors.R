#!/usr/bin/env Rscript

# Headless R execution bridge for Bulk RNA-seq Normalization v1
# Executes official Bioconductor DESeq2::estimateSizeFactorsForMatrix

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(input = NULL, method = "ratio", output = NULL)
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--input" && i < length(args)) {
      params$input <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--method" && i < length(args)) {
      params$method <- args[i + 1]
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

if (is.null(params$input) || is.null(params$output)) {
  cat("Usage: Rscript estimate_size_factors.R --input <counts.tsv> --method <ratio|poscounts> --output <result.json>\n", file = stderr())
  quit(status = 1)
}

timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
r_ver <- as.character(getRversion())
os_info <- as.character(Sys.info()["sysname"])

write_output_json <- function(status, method, size_factors, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, error_message, out_path) {
  # Build JSON string deterministically
  sf_json <- "null"
  if (!is.null(size_factors)) {
    sf_entries <- sapply(names(size_factors), function(nm) {
      val <- size_factors[[nm]]
      if (is.na(val) || is.nan(val)) "null" else sprintf("\"%s\": %s", nm, as.character(val))
    })
    sf_json <- paste0("{", paste(sf_entries, collapse = ", "), "}")
  }
  
  err_json <- if (is.null(error_message)) "null" else sprintf("\"%s\"", gsub("\"", "\\\\\"", error_message))
  
  json_str <- sprintf(
    '{\n  "status": "%s",\n  "method": "%s",\n  "size_factors": %s,\n  "r_version": "%s",\n  "bioc_version": "%s",\n  "deseq2_version": "%s",\n  "operating_system": "%s",\n  "execution_timestamp": "%s",\n  "error_message": %s\n}\n',
    status, method, sf_json, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, err_json
  )
  cat(json_str, file = out_path)
}

# Check DESeq2 availability
deseq2_available <- suppressWarnings(requireNamespace("DESeq2", quietly = TRUE))
deseq2_ver <- if (deseq2_available) as.character(packageVersion("DESeq2")) else "unavailable"

bioc_ver <- "unavailable"
if (suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
  bioc_ver <- as.character(BiocManager::version())
} else if (deseq2_available) {
  # Infer Bioconductor release from DESeq2 1.46.0 -> Bioc 3.20
  bioc_ver <- "3.20"
}

if (!deseq2_available) {
  write_output_json("ERROR", params$method, NULL, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, "R package 'DESeq2' is not available in the active environment.", params$output)
  quit(status = 2)
}

# Read count matrix
counts <- tryCatch({
  df <- read.delim(params$input, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  as.matrix(df)
}, error = function(e) {
  write_output_json("ERROR", params$method, NULL, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, paste("Failed to read input matrix:", e$message), params$output)
  quit(status = 3)
})

# Run DESeq2::estimateSizeFactorsForMatrix
res <- tryCatch({
  sf <- DESeq2::estimateSizeFactorsForMatrix(counts, type = params$method)
  list(success = TRUE, sf = sf)
}, error = function(e) {
  list(success = FALSE, message = e$message)
})

if (res$success) {
  write_output_json("SUCCESS", params$method, res$sf, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, NULL, params$output)
  quit(status = 0)
} else {
  write_output_json("ERROR", params$method, NULL, r_ver, bioc_ver, deseq2_ver, os_info, timestamp, res$message, params$output)
  quit(status = 4)
}
