#!/usr/bin/env Rscript

# Headless R execution bridge for Stage 8C Pathway and Signaling Analysis
# Executes preranked GSEA using clusterProfiler, org.Hs.eg.db, fgsea, and KEGG

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  params <- list(
    ranks = NULL,
    output = NULL,
    min_size = 15,
    max_size = 500,
    pvalue_cutoff = 1.0,
    kegg_rds = "data/references/kegg_hsa_pathways.rds",
    seed = 42,
    gsea_param = 1.0,
    score_type = "std",
    nperm_simple = 1000,
    eps = 0,
    padjust_method = "BH"
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--ranks" && i < length(args)) {
      params$ranks <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      params$output <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--min_size" && i < length(args)) {
      params$min_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--max_size" && i < length(args)) {
      params$max_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--pvalue_cutoff" && i < length(args)) {
      params$pvalue_cutoff <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--kegg_rds" && i < length(args)) {
      params$kegg_rds <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--seed" && i < length(args)) {
      params$seed <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--gsea_param" && i < length(args)) {
      params$gsea_param <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--score_type" && i < length(args)) {
      params$score_type <- as.character(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--nperm_simple" && i < length(args)) {
      params$nperm_simple <- as.integer(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--eps" && i < length(args)) {
      params$eps <- as.numeric(args[i + 1])
      i <- i + 2
    } else if (args[i] == "--padjust_method" && i < length(args)) {
      params$padjust_method <- as.character(args[i + 1])
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

params <- parse_args(args)

if (is.null(params$ranks) || is.null(params$output)) {
  cat("Usage: Rscript run_preranked_gsea.R --ranks <ranks.tsv> --output <output.json> [--min_size 15] [--max_size 500] [--pvalue_cutoff 1.0] [--kegg_rds <path>] [--seed 42] [--gsea_param 1] [--score_type std] [--nperm_simple 1000]\n", file = stderr())
  quit(status = 1)
}

timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
r_ver <- as.character(getRversion())
os_info <- as.character(Sys.info()["sysname"])

write_error <- function(msg, code = 2) {
  err_obj <- list(
    status = "ERROR",
    r_version = r_ver,
    operating_system = os_info,
    execution_timestamp = timestamp,
    error_message = msg
  )
  jsonlite::write_json(err_obj, params$output, auto_unbox = TRUE, pretty = TRUE)
  cat(msg, "\n", file = stderr())
  quit(status = code)
}

# Check package dependencies
req_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "GO.db", "fgsea", "KEGGREST", "AnnotationDbi")
missing_pkgs <- req_pkgs[!sapply(req_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  write_error(paste("Missing required R packages:", paste(missing_pkgs, collapse = ", ")), code = 2)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(GO.db)
  library(fgsea)
  library(KEGGREST)
  library(AnnotationDbi)
})

pkg_versions <- list(
  clusterProfiler = as.character(packageVersion("clusterProfiler")),
  org_Hs_eg_db = as.character(packageVersion("org.Hs.eg.db")),
  GO_db = as.character(packageVersion("GO.db")),
  fgsea = as.character(packageVersion("fgsea")),
  KEGGREST = as.character(packageVersion("KEGGREST"))
)

# Read ranks
ranks_df <- tryCatch({
  read.delim(params$ranks, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}, error = function(e) {
  write_error(paste("Failed to read ranks file:", e$message), code = 3)
})

if (!all(c("entrez_gene_id", "ranking_value") %in% colnames(ranks_df))) {
  write_error("Ranks file must contain 'entrez_gene_id' and 'ranking_value' columns", code = 3)
}

# Clean and sort geneList
valid_idx <- !is.na(ranks_df$entrez_gene_id) & ranks_df$entrez_gene_id != "" & !is.na(ranks_df$ranking_value)
ranks_df <- ranks_df[valid_idx, ]

geneList <- as.numeric(ranks_df$ranking_value)
names(geneList) <- as.character(ranks_df$entrez_gene_id)

# Order decreasingly
ord <- order(geneList, decreasing = TRUE)
geneList <- geneList[ord]
ranks_df <- ranks_df[ord, ]

# Symbol lookup map
approved_symbols <- if ("approved_symbol" %in% colnames(ranks_df)) ranks_df$approved_symbol else rep("", nrow(ranks_df))
entrez2symbol <- setNames(approved_symbols, as.character(ranks_df$entrez_gene_id))

set.seed(params$seed)

# Helper function to parse clusterProfiler result to structured list
parse_gsea_results <- function(res_obj, db_name, e2s_map) {
  if (is.null(res_obj)) return(list())
  df <- tryCatch(as.data.frame(res_obj), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(list())

  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    core_genes <- if (!is.null(row$core_enrichment) && !is.na(row$core_enrichment) && row$core_enrichment != "") {
      unlist(strsplit(as.character(row$core_enrichment), "/"))
    } else {
      character(0)
    }
    core_symbols <- unname(e2s_map[core_genes])
    core_symbols[is.na(core_symbols)] <- ""

    nes_val <- if (!is.na(row$NES)) as.numeric(row$NES) else 0.0
    direction <- if (nes_val >= 0) "POSITIVE" else "NEGATIVE"
    padj_val <- if (!is.na(row$p.adjust)) as.numeric(row$p.adjust) else 1.0
    pval_val <- if (!is.na(row$pvalue)) as.numeric(row$pvalue) else 1.0
    es_val <- if (!is.na(row$enrichmentScore)) as.numeric(row$enrichmentScore) else 0.0
    psize_val <- if (!is.na(row$setSize)) as.integer(row$setSize) else 0

    out[[i]] <- list(
      source_database = db_name,
      pathway_id = as.character(row$ID),
      pathway_name = as.character(row$Description),
      enrichment_score = es_val,
      normalized_enrichment_score = nes_val,
      p_value = pval_val,
      adjusted_p_value = padj_val,
      pathway_size = psize_val,
      direction = direction,
      leading_edge_gene_ids = as.list(core_genes),
      leading_edge_symbols = as.list(core_symbols),
      ranking_metric = "wald_statistic",
      analysis_version = paste0("clusterProfiler_", pkg_versions$clusterProfiler),
      is_significant = (padj_val < 0.05)
    )
  }
  return(out)
}

# 1. Run GO Biological Process GSEA with explicit parameters
cat("Executing GO Biological Process preranked GSEA...\n")
go_res <- tryCatch({
  clusterProfiler::gseGO(
    geneList = geneList,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    exponent = params$gsea_param,
    minGSSize = params$min_size,
    maxGSSize = params$max_size,
    pvalueCutoff = params$pvalue_cutoff,
    pAdjustMethod = params$padjust_method,
    eps = params$eps,
    scoreType = params$score_type,
    nPermSimple = params$nperm_simple,
    seed = TRUE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Warning: GO GSEA encountered an error:", e$message, "\n", file = stderr())
  NULL
})

go_list <- parse_gsea_results(go_res, "GO_BP", entrez2symbol)
cat("GO BP completed. Evaluated pathways:", length(go_list), "\n")

# 2. Run KEGG GSEA with explicit parameters
cat("Executing KEGG preranked GSEA...\n")
kegg_res <- NULL
if (!is.null(params$kegg_rds) && file.exists(params$kegg_rds)) {
  cat("Loading governed KEGG annotations from:", params$kegg_rds, "\n")
  kegg_data <- tryCatch(readRDS(params$kegg_rds), error = function(e) NULL)
  if (!is.null(kegg_data) && !is.null(kegg_data$KEGGPATHID2EXTID) && !is.null(kegg_data$KEGGPATHID2NAME)) {
    kegg_res <- tryCatch({
      clusterProfiler::GSEA(
        geneList = geneList,
        TERM2GENE = kegg_data$KEGGPATHID2EXTID,
        TERM2NAME = kegg_data$KEGGPATHID2NAME,
        exponent = params$gsea_param,
        minGSSize = params$min_size,
        maxGSSize = params$max_size,
        pvalueCutoff = params$pvalue_cutoff,
        pAdjustMethod = params$padjust_method,
        eps = params$eps,
        scoreType = params$score_type,
        nPermSimple = params$nperm_simple,
        seed = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      cat("Warning: KEGG GSEA via TERM2GENE encountered an error:", e$message, "\n", file = stderr())
      NULL
    })
  }
}

if (is.null(kegg_res)) {
  cat("Falling back to clusterProfiler::gseKEGG online...\n")
  kegg_res <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList = geneList,
      organism = "hsa",
      keyType = "kegg",
      exponent = params$gsea_param,
      minGSSize = params$min_size,
      maxGSSize = params$max_size,
      pvalueCutoff = params$pvalue_cutoff,
      pAdjustMethod = params$padjust_method,
      eps = params$eps,
      scoreType = params$score_type,
      nPermSimple = params$nperm_simple,
      seed = TRUE,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("Warning: KEGG gseKEGG encountered an error:", e$message, "\n", file = stderr())
    NULL
  })
}

kegg_list <- parse_gsea_results(kegg_res, "KEGG", entrez2symbol)
cat("KEGG completed. Evaluated pathways:", length(kegg_list), "\n")

# Write output JSON
out_obj <- list(
  status = "OK",
  r_version = r_ver,
  operating_system = os_info,
  execution_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  packages = pkg_versions,
  gsea_parameters = list(
    exponent = params$gsea_param,
    gsea_param = params$gsea_param,
    score_type = params$score_type,
    min_size = params$min_size,
    max_size = params$max_size,
    pvalue_cutoff = params$pvalue_cutoff,
    padjust_method = params$padjust_method,
    eps = params$eps,
    seed = params$seed,
    nperm_simple = params$nperm_simple
  ),
  go_results = go_list,
  kegg_results = kegg_list
)

jsonlite::write_json(out_obj, params$output, auto_unbox = TRUE, pretty = TRUE)
cat("Wrote GSEA results to:", params$output, "\n")
