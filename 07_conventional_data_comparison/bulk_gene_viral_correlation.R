#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Bulk time course correlation of host genes with viral load
# ------------------------------------------------------------------------------
# Rationale
#   scViralTracer-seq profiles one cell at a time, so the number of cells is
#   small. A conventional bulk RNA-seq time course covers the same infection
#   with many more replicates per time point and therefore provides an
#   independent readout. If a host gene tracks the viral load in the bulk time
#   course as well, the within-cell correlation is not an artefact of the
#   single-cell platform.
#
# Method
#   1. Read a gene-by-sample FPKM matrix of a bulk infection time course.
#   2. Parse the hours post infection out of each sample name and order the
#      samples chronologically.
#   3. Take the viral NP transcript as the quantitative proxy for viral load.
#   4. Compute the Pearson correlation between every host gene and NP across
#      the time course, and correct the p-values with the Benjamini-Hochberg
#      procedure.
#   5. Retain genes above the correlation and FDR thresholds, and report the
#      correlation of a set of query genes for a direct comparison with the
#      single-cell results.
#
#   The correlation is computed as one matrix product rather than with a loop
#   over genes. For a gene-by-time matrix G and a viral vector v, centring both
#   and taking the inner product yields every Pearson coefficient at once, with
#   the same numbers a per-gene cor.test loop would give. No parallel backend is
#   therefore needed.
#
# Input
#   A gene-by-sample matrix in FPKM. The first column holds the gene names and
#   every remaining column is one bulk library. Sample names must contain the
#   hours post infection, for example A0hpi, A3hpi, A4hpi, A5hpi, A6hpi.
#   Supported formats are .xlsx, .csv and .tsv.
#
# Output (written to --output-dir)
#   bulk_all_gene_correlations.csv   one row per gene with r, p-value, q-value
#                                    and the number of time points used
#   bulk_significant_genes.csv       the retained genes only
#   bulk_query_genes.csv             the genes named in --query-genes
#   bulk_correlation_summary.txt     settings and headline numbers
#
# Requirements
#   Base R only. Optional: readxl for .xlsx input.
#
# Usage
#   Rscript bulk_gene_viral_correlation.R \
#           --fpkm /path/to/bulk_fpkm.xlsx \
#           --output-dir ./results_bulk
#
# Optional flags
#   --viral-gene NAME        row name of the viral load proxy (default WSN_NP)
#   --time-pattern REGEX     Perl regular expression whose first capture group
#                            is the hours post infection
#                            (default (?i)([0-9]+)\s*hpi)
#   --time-values LIST       comma separated hours, used instead of parsing the
#                            sample names. The order must match the column order
#   --min-corr 0.95          minimum Pearson r for a gene to be retained
#   --direction positive     positive, negative or both
#   --max-fdr 0.05           maximum Benjamini-Hochberg q-value
#   --query-genes LIST       comma separated genes to report individually
#                            (default IFRD2)
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript bulk_gene_viral_correlation.R --fpkm FILE --output-dir DIR
       [--viral-gene WSN_NP] [--time-pattern REGEX] [--time-values 0,3,4,5,6]
       [--min-corr 0.95] [--direction positive|negative|both] [--max-fdr 0.05]
       [--query-genes IFRD2]"

opt <- parse_cli(list(
  fpkm           = NULL,
  `output-dir`   = NULL,
  `viral-gene`   = "WSN_NP",
  `time-pattern` = "(?i)([0-9]+)\\s*hpi",
  `time-values`  = "",
  `min-corr`     = 0.95,
  direction      = "positive",
  `max-fdr`      = 0.05,
  `query-genes`  = "IFRD2"
), USAGE)

if (!opt$direction %in% c("positive", "negative", "both")) {
  stop("--direction must be one of positive, negative, both", call. = FALSE)
}

ensure_dir(opt[["output-dir"]])

# --------------------------------------------------------------- 1. read in ---
step("Reading bulk matrix: ", opt$fpkm)
expr <- as_numeric_matrix(read_matrix(opt$fpkm))
samples <- colnames(expr)
note(nrow(expr), " genes x ", length(samples), " libraries")

if (ncol(expr) < 3) {
  stop("A correlation across a time course needs at least three libraries",
       call. = FALSE)
}

# ------------------------------------------------------- 2. order by time -----
# The hours post infection are read out of the sample name so that the script
# does not depend on the column order of the spreadsheet. A gene and the viral
# proxy are then compared over exactly the same ordered libraries.
if (nzchar(opt[["time-values"]])) {
  times <- as.numeric(split_csv(opt[["time-values"]]))
  if (length(times) != length(samples)) {
    stop("--time-values supplies ", length(times), " values but the matrix has ",
         length(samples), " libraries", call. = FALSE)
  }
} else {
  m <- regmatches(samples, regexec(opt[["time-pattern"]], samples, perl = TRUE))
  hit <- lengths(m) > 0L
  if (!all(hit)) {
    stop("These sample names carry no time point for --time-pattern '",
         opt[["time-pattern"]], "': ",
         paste(utils::head(samples[!hit], 5), collapse = ", "),
         "\nPass --time-values instead", call. = FALSE)
  }
  times <- as.numeric(vapply(m, function(x) x[[2]], character(1)))
}

ord <- order(times)
samples <- samples[ord]
times <- times[ord]
expr <- expr[, samples, drop = FALSE]

step("Time course: ", paste(sprintf("%shpi", times), collapse = " -> "))
note("Libraries: ", paste(samples, collapse = ", "))

if (any(duplicated(times))) {
  note("Some time points are replicated. Each replicate contributes one point ",
       "to the correlation")
}

# --------------------------------------------------- 3. locate the viral proxy -
viral_gene <- opt[["viral-gene"]]
if (!viral_gene %in% rownames(expr)) {
  stop("Viral gene '", viral_gene, "' is not a row of the matrix", call. = FALSE)
}
v <- as.numeric(expr[viral_gene, ])
if (!all(is.finite(v))) {
  stop("The viral load proxy contains missing values", call. = FALSE)
}
if (stats::sd(v) == 0) {
  stop("The viral load proxy is constant across the time course, so no ",
       "correlation can be computed", call. = FALSE)
}
step("Viral load proxy: ", viral_gene)

host_genes <- setdiff(rownames(expr), viral_gene)
G <- as.matrix(expr[host_genes, , drop = FALSE])
G[!is.finite(G)] <- NA_real_

# ------------------------------------------------------------- 4. correlation -
step("Correlating ", length(host_genes), " host genes with ", viral_gene)

# Genes with any missing value are dropped, because a partial series would be
# correlated over fewer points than the viral proxy and the resulting p-values
# would not be comparable. On a gene-by-time matrix, complete.cases returns one
# flag per gene.
complete <- stats::complete.cases(G)
n_dropped <- sum(!complete)
if (n_dropped > 0) {
  note(n_dropped, " gene(s) dropped for missing values")
}
G <- G[complete, , drop = FALSE]

vc <- v - mean(v)
denom_v <- sqrt(sum(vc^2))

Gc <- G - rowMeans(G)
denom_g <- sqrt(rowSums(Gc^2))
r <- as.numeric(Gc %*% vc) / (denom_g * denom_v)

# A gene that is flat across the whole time course has zero variance, so its
# correlation is undefined rather than zero.
flat <- denom_g == 0
r[flat | !is.finite(r)] <- NA_real_
r[r > 1] <- 1
r[r < -1] <- -1
if (sum(flat) > 0) {
  note(sum(flat), " gene(s) are constant across the time course and receive NA")
}

n_t <- length(v)
df_res <- max(n_t - 2L, 1L)
t_stat <- r * sqrt(df_res / pmax(1 - r^2, .Machine$double.eps))
p_val <- 2 * stats::pt(-abs(t_stat), df = df_res)
p_val[!is.finite(r)] <- NA_real_

results <- data.frame(
  gene_name    = rownames(G),
  correlation  = r,
  p.value      = p_val,
  n_timepoints = rep.int(as.integer(n_t), nrow(G)),
  stringsAsFactors = FALSE
)
results$fdr <- bh_fdr(results$p.value)

# --------------------------------------------------------------- 5. filtering -
tested <- !is.na(results$correlation) & !is.na(results$fdr) &
          results$n_timepoints == n_t

keep_dir <- switch(opt$direction,
  positive = results$correlation >  opt[["min-corr"]],
  negative = results$correlation < -opt[["min-corr"]],
  both     = abs(results$correlation) > opt[["min-corr"]])

significant <- results[tested & keep_dir & results$fdr < opt[["max-fdr"]], ]
significant <- significant[order(-abs(significant$correlation)), ]

step("Genes above the threshold: ", nrow(significant), " of ",
     sum(tested), " tested")
note("Criteria: ", opt$direction, " r ",
     if (opt$direction == "negative") "< -" else "> ", opt[["min-corr"]],
     ", FDR < ", opt[["max-fdr"]])

if (nrow(significant) > 0) {
  note("Top genes:")
  top <- utils::head(significant, 10)
  for (i in seq_len(nrow(top))) {
    note(sprintf("  %-16s r = %6.3f   FDR = %9.2e",
                 top$gene_name[i], top$correlation[i], top$fdr[i]))
  }
}

write_table(results[order(-results$correlation), ],
            file.path(opt[["output-dir"]], "bulk_all_gene_correlations.csv"))
write_table(significant,
            file.path(opt[["output-dir"]], "bulk_significant_genes.csv"))

# -------------------------------------------------------- 6. query gene table -
# This is the table used to compare platforms directly. A gene that correlates
# with the viral load inside single cells is looked up here to see whether the
# bulk time course agrees.
query <- split_csv(opt[["query-genes"]])
if (length(query) > 0) {
  idx <- match(query, results$gene_name)
  found <- !is.na(idx)
  if (any(!found)) {
    note("Query genes absent from the matrix: ",
         paste(query[!found], collapse = ", "))
  }
  query_df <- results[idx[found], , drop = FALSE]
  query_df$significant <- query_df$gene_name %in% significant$gene_name
  write_table(query_df,
              file.path(opt[["output-dir"]], "bulk_query_genes.csv"))

  if (nrow(query_df) > 0) {
    step("Query genes")
    for (i in seq_len(nrow(query_df))) {
      note(sprintf("  %-16s r = %6.3f   p = %9.2e   FDR = %9.2e",
                   query_df$gene_name[i], query_df$correlation[i],
                   query_df$p.value[i], query_df$fdr[i]))
    }
  }
}

# ------------------------------------------------------------------ 7. report -
summary_path <- file.path(opt[["output-dir"]], "bulk_correlation_summary.txt")
lines <- c(
  "Bulk time course correlation with viral load",
  paste(rep("=", 60), collapse = ""),
  paste0("Input matrix        : ", opt$fpkm),
  paste0("Viral load proxy    : ", viral_gene),
  paste0("Time points (hpi)   : ", paste(times, collapse = ", ")),
  paste0("Libraries           : ", paste(samples, collapse = ", ")),
  "",
  "Thresholds",
  paste0("  correlation       : ", opt$direction, " r ",
         if (opt$direction == "negative") "< -" else "> ", opt[["min-corr"]]),
  paste0("  adjusted p-value  : Benjamini-Hochberg q < ", opt[["max-fdr"]]),
  "",
  "Counts",
  paste0("  genes in matrix   : ", nrow(expr)),
  paste0("  genes tested      : ", sum(tested)),
  paste0("  genes retained    : ", nrow(significant)),
  "",
  "Distribution of r over tested genes",
  paste0("  minimum           : ", fmt_num(min(results$correlation[tested]))),
  paste0("  median            : ", fmt_num(stats::median(results$correlation[tested]))),
  paste0("  mean              : ", fmt_num(mean(results$correlation[tested]))),
  paste0("  maximum           : ", fmt_num(max(results$correlation[tested])))
)

if (exists("query_df") && nrow(query_df) > 0) {
  lines <- c(lines, "", "Query genes",
             sprintf("  %-16s %8s %12s %12s", "gene", "r", "p", "FDR"),
             vapply(seq_len(nrow(query_df)), function(i)
               sprintf("  %-16s %8.3f %12.2e %12.2e",
                       query_df$gene_name[i], query_df$correlation[i],
                       query_df$p.value[i], query_df$fdr[i]),
               character(1)))
}

if (nrow(significant) > 0) {
  top <- utils::head(significant, 20)
  lines <- c(lines, "", paste0("Top ", nrow(top), " retained genes"),
             vapply(seq_len(nrow(top)), function(i)
               sprintf("  %-16s r = %6.3f   FDR = %9.2e",
                       top$gene_name[i], top$correlation[i], top$fdr[i]),
               character(1)))
}

writeLines(lines, summary_path)
message("  written: ", summary_path)

step("Done")
