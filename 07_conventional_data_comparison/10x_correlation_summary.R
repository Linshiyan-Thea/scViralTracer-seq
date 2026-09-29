#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Droplet based single-cell reference data as a platform control
# ------------------------------------------------------------------------------
# Rationale
#   scViralTracer-seq follows one identified cell over time. A droplet based
#   experiment profiles each cell once and therefore cannot link the state of a
#   cell before infection to the viral load it later supports. The two platforms
#   nevertheless share one measurement, namely how strongly each host gene
#   tracks the viral load across the infected population. Reporting that number
#   on a conventional droplet dataset puts the within-cell correlation obtained
#   by scViralTracer-seq into context.
#
# Method
#   Two modes are supported.
#
#   --correlations FILE
#       Read a table of per-gene correlations that an upstream droplet pipeline
#       has already produced. Column names are detected automatically, so the
#       common spellings of the gene, the coefficient and the p-value columns
#       all work. Benjamini-Hochberg q-values are computed unless the table
#       already carries them.
#
#   --matrix FILE --viral-gene NAME
#       Compute the correlations directly from a gene-by-cell expression matrix.
#       Every cell is one observation, so this is a across-cell correlation
#       rather than the within-cell correlation of script 05.
#
#   In both modes the retained genes are those above the correlation and FDR
#   thresholds. The correlation of a set of query genes is reported separately so
#   that the two platforms can be compared gene by gene.
#
# Input format for --matrix
#   Genes are rows and cells are columns. Supported formats are .xlsx, .csv and
#   .tsv. Cell names are arbitrary because each cell contributes exactly one
#   observation.
#
# Output (written to --output-dir)
#   tenx_all_gene_correlations.csv   one row per gene with r, p-value, q-value
#   tenx_significant_genes.csv       the retained genes only
#   tenx_top_correlated_genes.csv    the strongest positive and negative genes
#   tenx_query_genes.csv             the genes named in --query-genes
#   tenx_correlation_report.txt      settings and headline numbers
#
# Requirements
#   Base R only. Optional: readxl for .xlsx input.
#
# Usage
#   Rscript tenx_correlation_summary.R \
#           --correlations /path/to/all_genes_correlation_with_NP.csv \
#           --output-dir ./results_tenx
#
#   Rscript tenx_correlation_summary.R \
#           --matrix /path/to/tenx_fpkm.xlsx \
#           --viral-gene NP \
#           --output-dir ./results_tenx
#
# Optional flags
#   --min-abs-corr 0.5     minimum absolute Pearson r for a gene to be retained
#   --max-fdr 0.05         maximum Benjamini-Hochberg q-value
#   --rank-limit 20        how many genes to list per direction in the top table
#   --query-genes IFRD2    comma separated genes to report individually
#   --viral-pattern REGEX  used to find the viral row when --viral-gene is not
#                          given (default _NP$ or ^NP$)
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript tenx_correlation_summary.R --output-dir DIR
       ( --correlations FILE | --matrix FILE [--viral-gene NAME] )
       [--min-abs-corr 0.5] [--max-fdr 0.05] [--rank-limit 20]
       [--query-genes IFRD2]"

opt <- parse_cli(list(
  correlations   = "",
  matrix         = "",
  `viral-gene`   = "",
  `viral-pattern` = "(_NP$|^NP$)",
  `output-dir`   = NULL,
  `min-abs-corr` = 0.5,
  `max-fdr`      = 0.05,
  `rank-limit`   = 20L,
  `query-genes`  = "IFRD2"
), USAGE)

if (!xor(nzchar(opt$correlations), nzchar(opt$matrix))) {
  stop("Supply either --correlations or --matrix, not both and not neither",
       call. = FALSE)
}

ensure_dir(opt[["output-dir"]])

# ---------------------------------------------------- column name detection ---
# Upstream droplet pipelines are not consistent about how they name these three
# columns, and a silent mismatch would produce a table full of NA. The first
# match in each candidate list wins and the choice is reported.
pick_column <- function(names, candidates, what) {
  for (cand in candidates) {
    hit <- which(tolower(names) == tolower(cand))
    if (length(hit) > 0) {
      note("Column for ", what, ": '", names[hit[[1]]], "'")
      return(names[hit[[1]]])
    }
  }
  stop("No column of this table looks like a ", what, ". Found: ",
       paste(names, collapse = ", "), call. = FALSE)
}

# -------------------------------------------- mode A. pre-computed table ------
# -------------------------------------------- mode B. compute from a matrix ---
if (nzchar(opt$correlations)) {
  step("Reading a pre-computed correlation table: ", opt$correlations)
  raw <- read_matrix(opt$correlations, first_col_is_name = FALSE)

  gene_col <- pick_column(colnames(raw),
                          c("gene_name", "gene_symbol", "gene", "symbol", "id"),
                          "the gene name")
  cor_col <- pick_column(colnames(raw),
                         c("correlation_with_NP", "correlation_with_np",
                           "correlation", "pearson_r", "pearson", "r", "cor",
                           "estimate"),
                         "the correlation coefficient")
  p_col <- tryCatch(
    pick_column(colnames(raw),
                c("p_value", "p.value", "pvalue", "p", "p_val"),
                "the p-value"),
    error = function(e) "")
  fdr_col <- tryCatch(
    pick_column(colnames(raw), c("fdr", "q_value", "qvalue", "padj", "p_adjust",
                                 "p.adj", "adjusted_p_value"),
                "the adjusted p-value"),
    error = function(e) "")

  results <- data.frame(
    gene_name   = as.character(raw[[gene_col]]),
    correlation = suppressWarnings(as.numeric(raw[[cor_col]])),
    stringsAsFactors = FALSE
  )
  results$p.value <- if (nzchar(p_col)) {
    suppressWarnings(as.numeric(raw[[p_col]]))
  } else {
    rep(NA_real_, nrow(results))
  }

  if (nzchar(fdr_col)) {
    results$fdr <- suppressWarnings(as.numeric(raw[[fdr_col]]))
    note("Adjusted p-values taken from the input table")
  } else if (all(is.na(results$p.value))) {
    stop("The table carries neither p-values nor adjusted p-values, so no ",
         "multiple testing correction is possible", call. = FALSE)
  } else {
    results$fdr <- bh_fdr(results$p.value)
    note("Benjamini-Hochberg q-values computed from the supplied p-values")
  }
  results$n_cells <- NA_integer_

} else {
  step("Reading expression matrix: ", opt$matrix)
  expr <- as_numeric_matrix(read_matrix(opt$matrix))
  note(nrow(expr), " genes x ", ncol(expr), " cells")

  if (ncol(expr) < 5) {
    stop("An across-cell correlation needs more cells than this", call. = FALSE)
  }

  viral_gene <- opt[["viral-gene"]]
  if (!nzchar(viral_gene)) {
    hits <- grep(opt[["viral-pattern"]], rownames(expr), value = TRUE)
    if (length(hits) == 0) {
      stop("No row matched --viral-pattern '", opt[["viral-pattern"]],
           "'. Pass --viral-gene explicitly", call. = FALSE)
    }
    if (length(hits) > 1) {
      note("Multiple rows matched the viral pattern: ",
           paste(utils::head(hits, 5), collapse = ", "))
      note("Using the first one. Pass --viral-gene to be explicit")
    }
    viral_gene <- hits[[1]]
  }
  if (!viral_gene %in% rownames(expr)) {
    stop("Viral gene '", viral_gene, "' is not a row of the matrix", call. = FALSE)
  }
  step("Viral load proxy: ", viral_gene)

  v <- as.numeric(expr[viral_gene, ])
  G <- as.matrix(expr[setdiff(rownames(expr), viral_gene), , drop = FALSE])
  G[!is.finite(G)] <- NA_real_

  # Both the gene and the viral proxy must be present in a cell for that cell to
  # contribute. Dropping cells per gene keeps genes that are detected in most
  # cells instead of discarding them because of one empty droplet.
  compute_one <- function(i) {
    ok <- is.finite(G[i, ]) & is.finite(v)
    if (sum(ok) < 3L) return(c(r = NA_real_, p = NA_real_, n = sum(ok)))
    ct <- safe_pearson(G[i, ok], v[ok])
    c(r = unname(ct[["r"]]), p = unname(ct[["p"]]), n = sum(ok))
  }

  step("Correlating ", nrow(G), " genes across cells")
  res <- t(vapply(seq_len(nrow(G)), compute_one, numeric(3)))
  results <- data.frame(
    gene_name   = rownames(G),
    correlation = res[, "r"],
    p.value     = res[, "p"],
    n_cells     = as.integer(res[, "n"]),
    stringsAsFactors = FALSE
  )
  results$fdr <- bh_fdr(results$p.value)
}

# --------------------------------------------------------------- 3. filtering -
tested <- !is.na(results$correlation) & !is.na(results$fdr)
if (sum(!tested) > 0) {
  note(sum(!tested), " gene(s) without a usable coefficient are excluded")
}
results <- results[tested, ]
results <- results[order(-abs(results$correlation)), ]
rownames(results) <- NULL

significant <- results[abs(results$correlation) > opt[["min-abs-corr"]] &
                       results$fdr < opt[["max-fdr"]], ]

step("Genes above the threshold: ", nrow(significant), " of ", nrow(results))
note("Criteria: |r| > ", opt[["min-abs-corr"]], ", FDR < ", opt[["max-fdr"]])
note("Correlation range: ", fmt_num(min(results$correlation)), " to ",
     fmt_num(max(results$correlation)))

write_table(results,
            file.path(opt[["output-dir"]], "tenx_all_gene_correlations.csv"))
write_table(significant,
            file.path(opt[["output-dir"]], "tenx_significant_genes.csv"))

# ----------------------------------------------------------- 4. ranked table ---
rank_limit <- as.integer(opt[["rank-limit"]])
pos <- results[results$correlation > 0, ]
neg <- results[results$correlation < 0, ]
pos <- utils::head(pos[order(-pos$correlation), ], rank_limit)
neg <- utils::head(neg[order(neg$correlation), ], rank_limit)
pos$direction <- "Positive"
neg$direction <- "Negative"
top <- rbind(pos, neg)
top <- top[, c("direction", setdiff(colnames(top), "direction")), drop = FALSE]

write_table(top,
            file.path(opt[["output-dir"]], "tenx_top_correlated_genes.csv"))

if (nrow(pos) > 0) {
  step("Strongest positive correlations")
  for (i in seq_len(min(10L, nrow(pos)))) {
    note(sprintf("  %-16s r = %6.3f   FDR = %9.2e",
                 pos$gene_name[i], pos$correlation[i], pos$fdr[i]))
  }
}
if (nrow(neg) > 0) {
  step("Strongest negative correlations")
  for (i in seq_len(min(10L, nrow(neg)))) {
    note(sprintf("  %-16s r = %6.3f   FDR = %9.2e",
                 neg$gene_name[i], neg$correlation[i], neg$fdr[i]))
  }
}

# -------------------------------------------------------- 5. query gene table -
query <- split_csv(opt[["query-genes"]])
query_df <- NULL
if (length(query) > 0) {
  idx <- match(query, results$gene_name)
  found <- !is.na(idx)
  if (any(!found)) {
    note("Query genes absent from this dataset: ",
         paste(query[!found], collapse = ", "))
  }
  query_df <- results[idx[found], , drop = FALSE]
  query_df$significant <- query_df$gene_name %in% significant$gene_name
  write_table(query_df,
              file.path(opt[["output-dir"]], "tenx_query_genes.csv"))

  if (nrow(query_df) > 0) {
    step("Query genes")
    for (i in seq_len(nrow(query_df))) {
      note(sprintf("  %-16s r = %6.3f   p = %9.2e   FDR = %9.2e",
                   query_df$gene_name[i], query_df$correlation[i],
                   query_df$p.value[i], query_df$fdr[i]))
    }
  }
}

# ------------------------------------------------------------------ 6. report -
source_desc <- if (nzchar(opt$correlations)) opt$correlations else opt$matrix
mode_desc <- if (nzchar(opt$correlations)) {
  "pre-computed correlation table"
} else {
  paste0("correlations computed here against ", viral_gene)
}

lines <- c(
  "Droplet based single-cell reference data",
  paste(rep("=", 60), collapse = ""),
  paste0("Input               : ", source_desc),
  paste0("Mode                : ", mode_desc),
  "",
  "Thresholds",
  paste0("  correlation       : |r| > ", opt[["min-abs-corr"]]),
  paste0("  adjusted p-value  : Benjamini-Hochberg q < ", opt[["max-fdr"]]),
  "",
  "Counts",
  paste0("  genes tested      : ", nrow(results)),
  paste0("  genes retained    : ", nrow(significant)),
  paste0("  positive          : ", sum(results$correlation > 0)),
  paste0("  negative          : ", sum(results$correlation < 0)),
  "",
  "Distribution of r over tested genes",
  paste0("  minimum           : ", fmt_num(min(results$correlation))),
  paste0("  median            : ", fmt_num(stats::median(results$correlation))),
  paste0("  mean              : ", fmt_num(mean(results$correlation))),
  paste0("  maximum           : ", fmt_num(max(results$correlation)))
)

if (!is.null(query_df) && nrow(query_df) > 0) {
  lines <- c(lines, "", "Query genes",
             sprintf("  %-16s %8s %12s %12s", "gene", "r", "p", "FDR"),
             vapply(seq_len(nrow(query_df)), function(i)
               sprintf("  %-16s %8.3f %12.2e %12.2e",
                       query_df$gene_name[i], query_df$correlation[i],
                       query_df$p.value[i], query_df$fdr[i]),
               character(1)))
}

if (nrow(significant) > 0) {
  top_sig <- utils::head(significant, 20)
  lines <- c(lines, "", paste0("Top ", nrow(top_sig), " retained genes"),
             vapply(seq_len(nrow(top_sig)), function(i)
               sprintf("  %-16s r = %6.3f   FDR = %9.2e",
                       top_sig$gene_name[i], top_sig$correlation[i],
                       top_sig$fdr[i]),
               character(1)))
}

report_path <- file.path(opt[["output-dir"]], "tenx_correlation_report.txt")
writeLines(lines, report_path)
message("  written: ", report_path)

step("Done")
