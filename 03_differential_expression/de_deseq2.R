#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Differential expression analysis with DESeq2
# ------------------------------------------------------------------------------
# What it does
#   Runs the standard DESeq2 workflow on a gene-by-sample raw count matrix and
#   writes one table with a call for every gene.
#
#   Significance thresholds used throughout the manuscript
#       |log2 fold change| >= 1   and   adjusted p < 0.05
#   Both are exposed as options so that the same script can be reused with
#   different cut-offs, but the defaults reproduce the reported gene sets.
#
# Input
#   --counts      raw count matrix, genes in rows and samples in columns.
#                 .xlsx, .csv or .tsv are all accepted. Counts must be integers,
#                 not FPKM or otherwise normalised values.
#   --coldata     sample sheet with at least the columns
#                     sample   matching a column name of the count matrix
#                     group    the factor of interest
#                 Additional columns are kept and may be used for blocking.
#
# Output (written to --output-dir)
#   DESeq2_results.csv   every tested gene with baseMean, log2FoldChange, lfcSE,
#                        stat, pvalue, padj and a change column of Up / Down / Not
#   DEG_up.csv           up-regulated genes only, sorted by log2 fold change
#   DEG_down.csv         down-regulated genes only
#   DEG_summary.txt      counts of up and down genes and the settings used
#
# Requirements
#   Bioconductor: DESeq2
#   CRAN: readxl (only if the count matrix is .xlsx)
#
# Usage
#   Rscript de_deseq2.R \
#           --counts  /path/to/merged_count.xlsx \
#           --coldata /path/to/group_all.xlsx \
#           --factor group --alt High --ref Low \
#           --output-dir ./results_de
#
# Optional flags
#   --lfc-cutoff 1         absolute log2 fold change threshold
#   --padj-cutoff 0.05     adjusted p-value threshold
#   --blocking cell        add a blocking factor to the design, giving
#                          ~ blocking + group. Use this when the same cell is
#                          sampled at several time points, which is exactly the
#                          longitudinal design of scViralTracer-seq
#   --min-cells-expressed all|N
#                          keep a gene when it has a non-zero count in all
#                          samples (default) or in at least N samples
#   --min-total-count 20   keep a gene when its total count across samples
#                          exceeds this value
#   --independent-filtering TRUE|FALSE
#                          DESeq2 independent filtering (default TRUE)
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript de_deseq2.R --counts FILE --coldata FILE --factor NAME
       --alt LEVEL --ref LEVEL --output-dir DIR
       [--lfc-cutoff 1] [--padj-cutoff 0.05] [--blocking NAME]
       [--min-cells-expressed all|N] [--min-total-count 20]"

opt <- parse_cli(list(
  counts           = NULL,
  coldata          = NULL,
  factor           = NULL,
  alt              = NULL,
  ref              = NULL,
  `output-dir`     = NULL,
  `lfc-cutoff`     = 1,
  `padj-cutoff`    = 0.05,
  blocking         = "",
  `min-cells-expressed` = "all",
  `min-total-count`     = 20,
  `independent-filtering` = TRUE
), USAGE)

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("Package 'DESeq2' is required. Install it with\n",
       "  BiocManager::install('DESeq2')", call. = FALSE)
}

ensure_dir(opt[["output-dir"]])
lfc_cut  <- opt[["lfc-cutoff"]]
padj_cut <- opt[["padj-cutoff"]]

# --------------------------------------------------------------- 1. read in ---
step("Reading the count matrix: ", opt$counts)
cnt <- read_matrix(opt$counts)

# Duplicated gene names are summed rather than dropped, which is what the
# original analysis did. It keeps the row set identical to the annotation.
if (anyDuplicated(rownames(cnt))) {
  n_dup <- sum(duplicated(rownames(cnt)))
  note(n_dup, " duplicated gene names found and summed")
  cnt <- as_numeric_matrix(cnt)
  cnt <- as.data.frame(stats::aggregate(as.matrix(cnt),
                                        by = list(Gene = rownames(cnt)),
                                        FUN = sum, na.rm = TRUE))
  rownames(cnt) <- cnt$Gene
  cnt$Gene <- NULL
}

cnt <- as_numeric_matrix(cnt)
cnt[is.na(cnt)] <- 0

step("Reading the sample sheet: ", opt$coldata)
coldata <- as.data.frame(read_matrix(opt$coldata, first_col_is_name = FALSE))

for (need in c("sample", opt$factor)) {
  if (!need %in% colnames(coldata)) {
    stop("Column '", need, "' is missing from the sample sheet. Found: ",
         paste(colnames(coldata), collapse = ", "), call. = FALSE)
  }
}

blocking <- opt$blocking
if (nzchar(blocking) && !blocking %in% colnames(coldata)) {
  stop("Blocking factor '", blocking, "' is not a column of the sample sheet",
       call. = FALSE)
}

# --------------------------------------------------- 2. align names and order ---
# Sample identifiers are normalised on both sides so that trivial differences
# in case or in disallowed characters cannot silently drop a library.
clean <- function(x) tolower(make.names(as.character(x)))
colnames(cnt) <- clean(colnames(cnt))
coldata$.key <- clean(coldata$sample)

missing_in_counts <- setdiff(coldata$.key, colnames(cnt))
if (length(missing_in_counts)) {
  stop("These samples are listed in the sample sheet but absent from the count ",
       "matrix: ", paste(utils::head(missing_in_counts, 5), collapse = ", "),
       call. = FALSE)
}
dropped <- setdiff(colnames(cnt), coldata$.key)
if (length(dropped)) {
  note(length(dropped), " column(s) of the count matrix have no sample sheet ",
       "entry and are ignored")
  cnt <- cnt[, coldata$.key, drop = FALSE]
} else {
  cnt <- cnt[, coldata$.key, drop = FALSE]
}

# ------------------------------------------------------- 3. build the object ---
step("Building the DESeq2 data set")
counts_int <- round(as.matrix(cnt))
if (!all(counts_int >= 0)) {
  stop("The count matrix contains negative values. DESeq2 needs raw, ",
       "un-normalised counts", call. = FALSE)
}

cd <- data.frame(row.names = coldata$.key)
cd$group <- factor(coldata[[opt$factor]])
if (nzchar(blocking)) cd$block <- factor(coldata[[blocking]])

for (lvl in c(opt$alt, opt$ref)) {
  if (!lvl %in% levels(cd$group)) {
    stop("Level '", lvl, "' is not present in the '", opt$factor,
         "' column. Found: ", paste(levels(cd$group), collapse = ", "),
         call. = FALSE)
  }
}
# Setting the reference level explicitly makes the sign of the fold change
# unambiguous. Positive values mean higher in --alt.
cd$group <- relevel(cd$group, ref = opt$ref)

design <- if (nzchar(blocking)) ~ block + group else ~ group
note("Design: ", paste(deparse(design), collapse = " "))
note("Contrast: ", opt$alt, " versus ", opt$ref,
     " (positive log2 fold change means higher in ", opt$alt, ")")
note("Samples: ", nrow(cd), "   Genes: ", nrow(counts_int))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts_int,
                                      colData = cd,
                                      design = design)

# --------------------------------------------------------- 4. low expression ---
n_expressed <- rowSums(DESeq2::counts(dds) > 0)
min_cells <- if (identical(tolower(opt[["min-cells-expressed"]]), "all")) {
  ncol(dds)
} else {
  as.integer(opt[["min-cells-expressed"]])
}
if (is.na(min_cells)) {
  stop("--min-cells-expressed must be 'all' or a positive integer", call. = FALSE)
}
row_total <- rowSums(DESeq2::counts(dds))
keep <- (n_expressed >= min_cells) & (row_total > opt[["min-total-count"]])

step("Filtering low-expression genes")
note("Expressed in >= ", min_cells, " sample(s) and total count > ",
     opt[["min-total-count"]])
note("Retained ", sum(keep), " of ", nrow(dds), " genes")

dds <- dds[keep, ]
if (nrow(dds) < 2) {
  stop("Fewer than two genes survived filtering. Loosen --min-cells-expressed ",
       "or --min-total-count", call. = FALSE)
}

# ------------------------------------------------------------- 5. run DESeq2 ---
step("Running DESeq2")
suppressMessages(dds <- DESeq2::DESeq(dds))

res <- DESeq2::results(dds,
                       contrast = c("group", opt$alt, opt$ref),
                       alpha = padj_cut,
                       independentFiltering = opt[["independent-filtering"]])

DEG <- as.data.frame(res)
DEG$gene <- rownames(DEG)

# ------------------------------------------------------------- 6. call genes ---
# A gene is called only when both criteria are met. The direction is decided by
# the sign of the fold change alone, because the magnitude criterion has already
# been applied. Using `> 0` here also covers log2 fold changes of exactly 1,
# which a `> 1` test would have mislabelled as down-regulated.
sig <- !is.na(DEG$padj) & DEG$padj < padj_cut &
       !is.na(DEG$log2FoldChange) & abs(DEG$log2FoldChange) >= lfc_cut
DEG$change <- ifelse(sig, ifelse(DEG$log2FoldChange > 0, "Up", "Down"), "Not")
DEG$change <- factor(DEG$change, levels = c("Up", "Down", "Not"))

DEG <- DEG[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat",
               "pvalue", "padj", "change")]

n_up   <- sum(DEG$change == "Up")
n_down <- sum(DEG$change == "Down")
n_tested <- sum(!is.na(DEG$padj))

step("Result")
note("Genes tested        : ", n_tested)
note("Up-regulated        : ", n_up)
note("Down-regulated      : ", n_down)
note("Not significant     : ", n_tested - n_up - n_down)
print(table(DEG$change))

# --------------------------------------------------------------- 7. outputs ---
write_table(DEG, file.path(opt[["output-dir"]], "DESeq2_results.csv"))

up <- DEG[DEG$change == "Up", ]
up <- up[order(-up$log2FoldChange), ]
down <- DEG[DEG$change == "Down", ]
down <- down[order(down$log2FoldChange), ]
write_table(up,   file.path(opt[["output-dir"]], "DEG_up.csv"))
write_table(down, file.path(opt[["output-dir"]], "DEG_down.csv"))

sh <- DESeq2::sizeFactors(dds)
lines <- c(
  strrep("=", 64),
  "Differential expression summary",
  strrep("=", 64),
  "",
  paste0("Count matrix         : ", basename(opt$counts)),
  paste0("Sample sheet         : ", basename(opt$coldata)),
  paste0("Design               : ", paste(deparse(design), collapse = " ")),
  paste0("Contrast             : ", opt$alt, " versus ", opt$ref),
  paste0("Samples              : ", ncol(dds)),
  "",
  strrep("-", 64),
  "Filtering",
  strrep("-", 64),
  paste0("  genes in input                 : ", length(keep)),
  paste0("  expressed in >= N samples      : ", min_cells),
  paste0("  total count >                  : ", opt[["min-total-count"]]),
  paste0("  genes retained                 : ", nrow(dds)),
  "",
  strrep("-", 64),
  "Significance criteria",
  strrep("-", 64),
  paste0("  |log2 fold change| >=          : ", lfc_cut),
  paste0("  adjusted p-value  <            : ", padj_cut),
  paste0("  multiple testing correction    : Benjamini-Hochberg (DESeq2 padj)"),
  paste0("  independent filtering          : ", opt[["independent-filtering"]]),
  "",
  strrep("-", 64),
  "Outcome",
  strrep("-", 64),
  paste0("  genes tested                   : ", n_tested),
  paste0("  up-regulated                   : ", n_up),
  paste0("  down-regulated                 : ", n_down),
  paste0("  not significant                : ", n_tested - n_up - n_down),
  "",
  strrep("-", 64),
  "Library size factors",
  strrep("-", 64),
  paste0("  ", paste(sprintf("%s = %.4f", names(sh), sh), collapse = "\n  "))
)
out_txt <- file.path(opt[["output-dir"]], "DEG_summary.txt")
writeLines(lines, out_txt)
message("  written: ", out_txt)

step("Done")
