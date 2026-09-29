#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Within-cell coupling of host gene expression to viral load
# ------------------------------------------------------------------------------
# Rationale
#   Comparing populations at each time point cannot reveal whether a host gene
#   rises and falls together with the virus inside the same cell. Because
#   scViralTracer-seq samples one cell repeatedly, the correlation can be
#   computed within each cell along its own time course and then tested for
#   consistency across cells. A gene that tracks the virus in most cells, in
#   the same direction, is a candidate regulator rather than a bystander of the
#   average response.
#
# Method
#   1. Take the viral NP transcript as a quantitative proxy for viral load.
#   2. Within each cell, and using only that cell's own time points, compute
#      the Pearson correlation between the expression trajectory of every host
#      gene and the trajectory of NP.
#   3. Summarise the per-cell coefficients for each gene as a mean and a median
#      across cells.
#   4. Test, for each gene, whether the per-cell coefficients are greater than
#      zero using a one-sided one-sample t-test across cells, and correct the
#      resulting p-values with the Benjamini-Hochberg procedure.
#   5. Retain genes that pass the mean correlation, the FDR and the minimum
#      number of informative cells.
#
# Input
#   A gene-by-sample expression matrix in FPKM. Column names must follow the
#   <CellID>.<TimePoint> convention, for example WSN01.0, WSN01.3, WSN01.6.
#   All time points belonging to one cell must share the same <CellID>.
#
# Output (written to --output-dir)
#   per_cell_gene_correlations.csv    long table, one row per gene per cell
#   gene_correlation_summary.csv      one row per gene with the across-cell
#                                     summary and the corrected p-value
#   significant_genes.csv             the retained genes only
#   significant_gene_expression.xlsx  one sheet per retained gene holding the
#                                     paired gene and NP values per cell
#   coupling_summary.txt              settings and headline numbers
#
# Requirements
#   R packages: data.table
#   Optional: furrr and future, needed only when --workers is greater than 1
#   Optional: openxlsx, needed only to write the .xlsx workbook
#
# Usage
#   Rscript per_cell_gene_viral_coupling.R \
#           --fpkm /path/to/merged_fpkm.xlsx \
#           --output-dir ./results_coupling
#
# Optional flags
#   --viral-gene NAME      row name of the viral load proxy (default: the row
#                          matching --viral-pattern)
#   --viral-pattern REGEX  used to find the NP row automatically (default _NP$)
#   --min-mean-corr 0.7    minimum mean correlation across cells
#   --max-fdr 0.01         maximum Benjamini-Hochberg q-value
#   --min-cells 3          minimum number of cells with a usable coefficient
#   --min-timepoints 2     minimum time points a cell must contribute
#   --workers 1            parallel workers. One means sequential execution
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript per_cell_gene_viral_coupling.R --fpkm FILE --output-dir DIR
       [--viral-gene NAME | --viral-pattern REGEX] [--min-mean-corr 0.7]
       [--max-fdr 0.01] [--min-cells 3] [--min-timepoints 2] [--workers 1]"

opt <- parse_cli(list(
  fpkm             = NULL,
  `output-dir`     = NULL,
  `viral-gene`     = "",
  `viral-pattern`  = "_NP$",
  `min-mean-corr`  = 0.7,
  `max-fdr`        = 0.01,
  `min-cells`      = 3L,
  `min-timepoints` = 2L,
  workers          = 1L
), USAGE)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required", call. = FALSE)
}

ensure_dir(opt[["output-dir"]])

# --------------------------------------------------------------- 1. read in ---
step("Reading expression matrix: ", opt$fpkm)
expr <- as_numeric_matrix(read_matrix(opt$fpkm))
meta <- parse_samples(colnames(expr))
expr <- expr[, meta$sample, drop = FALSE]
meta <- meta[match(colnames(expr), meta$sample), , drop = FALSE]

genes <- rownames(expr)
note(length(genes), " genes x ", ncol(expr), " libraries x ",
     length(unique(meta$cell_id)), " cells")
note("Time points: ", paste(sort(unique(meta$time)), collapse = ", "))

# --------------------------------------------------- 2. locate the viral proxy -
viral_gene <- opt[["viral-gene"]]
if (!nzchar(viral_gene)) {
  hits <- grep(opt[["viral-pattern"]], genes, value = TRUE)
  if (length(hits) == 0) {
    stop("No row matched --viral-pattern '", opt[["viral-pattern"]],
         "'. Pass --viral-gene explicitly", call. = FALSE)
  }
  if (length(hits) > 1) {
    note("Multiple rows matched the viral pattern: ",
         paste(utils::head(hits, 5), collapse = ", "))
    note("Using the first one. Pass --viral-gene to be explicit")
  }
  viral_gene <- hits[1]
}
if (!viral_gene %in% genes) {
  stop("Viral gene '", viral_gene, "' is not a row of the matrix", call. = FALSE)
}
step("Viral load proxy: ", viral_gene)

host_genes <- setdiff(genes, viral_gene)

# ------------------------------------------------------------- 3. per-cell r ---
# The correlation is computed as a single matrix product per cell rather than
# with a loop over genes. For a gene matrix G of dimension genes x time points
# and a viral vector v, centring both and taking the inner product yields every
# Pearson coefficient at once. Genes with zero variance across the time points
# of that cell carry no information and are set to NA.
np_row <- as.numeric(expr[viral_gene, ])
names(np_row) <- colnames(expr)

cell_ids <- unique(meta$cell_id)

compute_one_cell <- function(cell_id) {
  cols <- meta$sample[meta$cell_id == cell_id]
  cols <- cols[order(meta$time[match(cols, meta$sample)])]
  if (length(cols) < opt[["min-timepoints"]]) return(NULL)

  # Drop time points where the viral proxy itself is missing.
  v <- np_row[cols]
  keep_t <- is.finite(v)
  if (sum(keep_t) < opt[["min-timepoints"]]) return(NULL)
  cols <- cols[keep_t]
  v <- v[keep_t]

  G <- as.matrix(expr[host_genes, cols, drop = FALSE])
  G[!is.finite(G)] <- NA_real_
  # complete.cases on a genes-by-time matrix returns one flag per gene.
  keep_g <- stats::complete.cases(G)
  G <- G[keep_g, , drop = FALSE]
  if (nrow(G) == 0) return(NULL)

  Gc <- G - rowMeans(G)
  vc <- v - mean(v)
  denom_g <- sqrt(rowSums(Gc^2))
  denom_v <- sqrt(sum(vc^2))
  r <- as.numeric(Gc %*% vc) / (denom_g * denom_v)
  # A flat trajectory has zero variance, so its correlation is undefined.
  r[denom_g == 0 | denom_v == 0 | !is.finite(r) | abs(r) > 1] <- NA_real_

  data.frame(gene_name = rownames(G),
             correlation = r,
             cell_id = cell_id,
             n_timepoints = length(cols),
             stringsAsFactors = FALSE)
}

step("Computing within-cell correlations across ", length(cell_ids), " cells")

workers <- as.integer(opt$workers)
if (workers > 1L) {
  if (!requireNamespace("furrr", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE)) {
    stop("--workers greater than 1 needs the packages 'furrr' and 'future'",
         call. = FALSE)
  }
  n_cores <- min(workers, parallel::detectCores())
  future::plan(future::multisession, workers = n_cores)
  note("Using ", n_cores, " parallel workers")
  on.exit(future::plan(future::sequential), add = TRUE)
  per_cell <- furrr::future_map(cell_ids, compute_one_cell,
                                .progress = TRUE,
                                .options = furrr::furrr_options(seed = TRUE))
} else {
  per_cell <- lapply(cell_ids, compute_one_cell)
}

per_cell <- per_cell[!vapply(per_cell, is.null, logical(1))]
if (length(per_cell) == 0) {
  stop("No cell contributed enough time points. Check the sample names and ",
       "--min-timepoints", call. = FALSE)
}
all_cors <- data.table::rbindlist(per_cell)
note(length(per_cell), " cells contributed, ", nrow(all_cors),
     " gene-cell coefficients")

write_table(as.data.frame(all_cors),
            file.path(opt[["output-dir"]], "per_cell_gene_correlations.csv"))

# ----------------------------------------------------- 4. across-cell summary ---
step("Summarising across cells")

gene_stats <- all_cors[, {
  v <- correlation[is.finite(correlation)]
  n <- length(v)
  if (n >= 2) {
    # One-sided test that the mean per-cell correlation is greater than zero.
    tt <- tryCatch(stats::t.test(v, alternative = "greater"),
                   error = function(e) NULL)
    p <- if (is.null(tt)) NA_real_ else tt$p.value
  } else if (n == 1) {
    p <- NA_real_
  } else {
    p <- NA_real_
  }
  list(
    mean_corr   = if (n > 0) mean(v) else NA_real_,
    median_corr = if (n > 0) stats::median(v) else NA_real_,
    sd_corr     = if (n >= 2) stats::sd(v) else NA_real_,
    min_corr    = if (n > 0) min(v) else NA_real_,
    max_corr    = if (n > 0) max(v) else NA_real_,
    n_cells     = as.integer(n),
    p_value     = p
  )
}, by = gene_name]

gene_stats <- gene_stats[!is.na(mean_corr) & !is.na(p_value) & n_cells >= 1]
gene_stats[, fdr := bh_fdr(p_value)]

significant <- gene_stats[
  mean_corr > opt[["min-mean-corr"]] &
  fdr < opt[["max-fdr"]] &
  n_cells >= opt[["min-cells"]]
]
significant <- significant[order(-mean_corr)]

n_total <- uniqueN(all_cors$gene_name)
step("Retained ", nrow(significant), " of ", n_total, " host genes")
note("Criteria: mean r > ", opt[["min-mean-corr"]],
     ", FDR < ", opt[["max-fdr"]],
     ", informative cells >= ", opt[["min-cells"]])

if (nrow(significant) > 0) {
  note("Top genes:")
  top <- utils::head(significant, 10)
  for (i in seq_len(nrow(top))) {
    note(sprintf("  %-14s mean r = %6.3f   median r = %6.3f   cells = %2d   FDR = %9.2e",
                 top$gene_name[i], top$mean_corr[i], top$median_corr[i],
                 top$n_cells[i], top$fdr[i]))
  }
}

write_table(as.data.frame(gene_stats[order(-mean_corr)]),
            file.path(opt[["output-dir"]], "gene_correlation_summary.csv"))
write_table(as.data.frame(significant),
            file.path(opt[["output-dir"]], "significant_genes.csv"))

# ------------------------------------------ 5. paired expression of retained ---
# This table is what allows a reader to reproduce the per-cell overlay of a
# host gene against the viral load without rerunning any correlation.
sig_genes <- significant$gene_name

if (length(sig_genes) > 0) {
  step("Exporting the paired expression values of the retained genes")

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    note("Package 'openxlsx' is not installed, writing a long CSV instead")
    long <- do.call(rbind, lapply(sig_genes, function(g) {
      do.call(rbind, lapply(cell_ids, function(cid) {
        cols <- meta$sample[meta$cell_id == cid]
        cols <- cols[order(meta$time[match(cols, meta$sample)])]
        data.frame(gene_name = g,
                   cell_id = cid,
                   time = meta$time[match(cols, meta$sample)],
                   gene_expression = as.numeric(expr[g, cols]),
                   np_expression = np_row[cols],
                   stringsAsFactors = FALSE)
      }))
    }))
    write_table(long,
                file.path(opt[["output-dir"]], "significant_gene_expression.csv"))
  } else {
    wb <- openxlsx::createWorkbook()
    used <- character(0)
    for (g in sig_genes) {
      rows <- lapply(cell_ids, function(cid) {
        cols <- meta$sample[meta$cell_id == cid]
        cols <- cols[order(meta$time[match(cols, meta$sample)])]
        data.frame(cell_id = cid,
                   time = meta$time[match(cols, meta$sample)],
                   gene_expression = as.numeric(expr[g, cols]),
                   np_expression = np_row[cols],
                   stringsAsFactors = FALSE)
      })
      gene_df <- do.call(rbind, rows)

      # Sheet names are limited to 31 characters and must stay unique.
      nm <- gsub("[^A-Za-z0-9_. -]", "_", substr(g, 1, 31))
      base_nm <- nm
      k <- 2L
      while (nm %in% used) {
        nm <- paste0(substr(base_nm, 1, 28), "_", k)
        k <- k + 1L
      }
      used <- c(used, nm)

      openxlsx::addWorksheet(wb, sheetName = nm)
      openxlsx::writeData(wb, sheet = nm, x = gene_df)
    }
    out_xlsx <- file.path(opt[["output-dir"]], "significant_gene_expression.xlsx")
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    message("  written: ", out_xlsx)
  }
}

# ------------------------------------------------------------- 6. report ------
lines <- c(
  strrep("=", 64),
  "Within-cell coupling of host genes to viral load",
  strrep("=", 64),
  "",
  paste0("Input matrix           : ", basename(opt$fpkm)),
  paste0("Viral load proxy       : ", viral_gene),
  paste0("Correlation            : Pearson, computed within each cell"),
  paste0("Across-cell test       : one-sided one-sample t-test, greater than zero"),
  paste0("Multiple testing       : Benjamini-Hochberg"),
  "",
  strrep("-", 64),
  "Data",
  strrep("-", 64),
  paste0("  genes in input       : ", length(genes)),
  paste0("  host genes tested    : ", length(host_genes)),
  paste0("  libraries            : ", ncol(expr)),
  paste0("  cells                : ", length(cell_ids)),
  paste0("  cells contributing   : ", length(per_cell)),
  paste0("  time points          : ",
         paste(sort(unique(meta$time)), collapse = ", ")),
  paste0("  min time points/cell : ", opt[["min-timepoints"]]),
  "",
  strrep("-", 64),
  "Retention criteria",
  strrep("-", 64),
  paste0("  mean r >             : ", opt[["min-mean-corr"]]),
  paste0("  FDR <                : ", opt[["max-fdr"]]),
  paste0("  informative cells >= : ", opt[["min-cells"]]),
  "",
  strrep("-", 64),
  "Outcome",
  strrep("-", 64),
  paste0("  genes with a usable summary : ", nrow(gene_stats)),
  paste0("  genes retained              : ", nrow(significant))
)
if (nrow(significant) > 0) {
  lines <- c(lines, "", "  Top 20 retained genes",
             sprintf("  %-14s %8s %8s %6s %10s",
                     "gene", "mean_r", "median_r", "cells", "FDR"))
  top <- utils::head(significant, 20)
  for (i in seq_len(nrow(top))) {
    lines <- c(lines, sprintf("  %-14s %8.3f %8.3f %6d %10.2e",
                              top$gene_name[i], top$mean_corr[i],
                              top$median_corr[i], top$n_cells[i], top$fdr[i]))
  }
}
out_txt <- file.path(opt[["output-dir"]], "coupling_summary.txt")
writeLines(lines, out_txt)
message("  written: ", out_txt)

step("Done")
