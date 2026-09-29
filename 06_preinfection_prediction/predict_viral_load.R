#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Predicting viral load from the pre-infection transcriptome
# ------------------------------------------------------------------------------
# Rationale
#   Conventional single-cell RNA sequencing destroys the cell, so the state a
#   cell was in before infection can never be linked to the viral load it later
#   supported. scViralTracer-seq records both from the same cell, which makes
#   it possible to ask whether the pre-infection transcriptome predicts how
#   much virus that cell will produce.
#
# Method
#   1. Take the pre-infection expression of every host gene as features and the
#      viral NP expression of the same cell at a late time point as the target.
#   2. Drop genes that are not expressed above a floor in enough cells.
#   3. Screen every gene on its own with a univariate linear regression against
#      the target, which is equivalent to a Pearson correlation test, and
#      correct the p-values with the Benjamini-Hochberg procedure.
#   4. Fit two multivariate models on the same feature set, an elastic net
#      regression and a random forest, and report both the apparent and the
#      leave-one-out performance so that the optimistic in-sample fit is never
#      mistaken for predictive accuracy.
#   5. Report feature importance from both multivariate models.
#
# Input
#   --fpkm        gene-by-cell matrix in FPKM holding the PRE-INFECTION
#                 expression of every cell. Column names must follow the
#                 <CellID>.<TimePoint> convention or be plain cell identifiers.
#   --np-row-6h   row name of that matrix carrying the viral load per cell at
#                 the late time point, for example WSN_NP_6hpi
#   --np-row-3h   optional row name carrying an intermediate viral load, kept
#                 in the output for trajectory inspection but never used as a
#                 feature or as a target
#
# Output (written to --output-dir)
#   single_gene_regression_results.csv  one row per gene with beta, r, r squared,
#                                       p-value, adjusted p-value and a
#                                       significant flag
#   model_performance.csv               apparent and leave-one-out RMSE and
#                                       R squared for both multivariate models
#   predictions.csv                     per-cell observed and predicted viral load
#   elasticnet_coefficients.csv         non-zero coefficients of the elastic net
#   randomforest_importance.csv         permutation importance of the random forest
#   all_results.xlsx                    the same tables in one workbook
#   prediction_summary.txt              settings, performance and top genes
#
# Requirements
#   CRAN: glmnet, ranger
#   Optional: openxlsx, needed only to write the .xlsx workbook
#
# Usage
#   Rscript predict_viral_load.R \
#           --fpkm /path/to/merged_fpkm_preinfection.xlsx \
#           --np-row-6h WSN_NP_6hpi \
#           --output-dir ./results_prediction
#
# Optional flags
#   --np-row-3h NAME       intermediate viral load row, retained but not modelled
#   --min-fpkm 1           expression floor for a gene to count as expressed
#   --min-cells 3          a gene is kept when at least this many cells are
#                          above the floor
#   --min-r 0.65           absolute Pearson r used to flag the top-correlated
#                          genes reported in the manuscript
#   --max-fdr 0.05         adjusted p-value used with --min-r
#   --alpha 0.5            elastic net mixing parameter. 0.5 balances the ridge
#                          and the lasso penalty
#   --n-trees 500          random forest trees
#   --min-node-size 1      random forest leaf size
#   --cv loocv|none        leave-one-out cross validation of the multivariate
#                          models. With a handful of cells this is the only
#                          honest estimate of performance
#   --seed 1               seed for the random forest
# ==============================================================================

# Load the shared helpers from the repository root.
.self <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(.self))), "R", "common.R"))
rm(.self)

USAGE <- "Usage: Rscript predict_viral_load.R --fpkm FILE --np-row-6h NAME
       --output-dir DIR [--np-row-3h NAME] [--min-fpkm 1] [--min-cells 3]
       [--min-r 0.65] [--max-fdr 0.05] [--alpha 0.5] [--n-trees 500]
       [--cv loocv|none] [--seed 1]"

opt <- parse_cli(list(
  fpkm            = NULL,
  `np-row-6h`     = NULL,
  `np-row-3h`     = "",
  `output-dir`    = NULL,
  `min-fpkm`      = 1,
  `min-cells`     = 3L,
  `min-r`         = 0.65,
  `max-fdr`       = 0.05,
  alpha           = 0.5,
  `n-trees`       = 500L,
  `min-node-size` = 1L,
  cv              = "loocv",
  seed            = 1L
), USAGE)

for (pkg in c("glmnet", "ranger")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required", call. = FALSE)
  }
}

ensure_dir(opt[["output-dir"]])
set.seed(opt$seed)

# Root mean squared error. Defined here rather than imported, because the only
# package that ships a function of this name for regression is one that the
# rest of the pipeline does not need.
rmse <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (sum(ok) == 0) return(NA_real_)
  sqrt(mean((observed[ok] - predicted[ok])^2))
}

r_squared <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (sum(ok) < 3) return(NA_real_)
  if (stats::sd(observed[ok]) == 0 || stats::sd(predicted[ok]) == 0) return(NA_real_)
  stats::cor(observed[ok], predicted[ok])^2
}

# --------------------------------------------------------------- 1. read in ---
step("Reading expression matrix: ", opt$fpkm)
raw <- read_matrix(opt$fpkm)
raw <- as_numeric_matrix(raw)

np_late_name <- opt[["np-row-6h"]]
np_mid_name <- opt[["np-row-3h"]]

if (!np_late_name %in% rownames(raw)) {
  stop("Row '", np_late_name, "' is not in the matrix. Pass the correct name ",
       "with --np-row-6h", call. = FALSE)
}

np_late <- as.numeric(raw[np_late_name, ])
np_mid <- if (nzchar(np_mid_name) && np_mid_name %in% rownames(raw)) {
  as.numeric(raw[np_mid_name, ])
} else {
  if (nzchar(np_mid_name)) {
    note("--np-row-3h '", np_mid_name, "' not found and ignored")
  }
  rep(NA_real_, ncol(raw))
}

cells <- colnames(raw)
keep_cell <- is.finite(np_late)
if (sum(!keep_cell) > 0) {
  note(sum(!keep_cell), " cell(s) without a viral load value are excluded: ",
       paste(utils::head(cells[!keep_cell], 5), collapse = ", "))
}
np_late <- np_late[keep_cell]
np_mid <- np_mid[keep_cell]
cells <- cells[keep_cell]

target_rows <- c(np_late_name, np_mid_name)
target_rows <- target_rows[nzchar(target_rows)]
feat <- raw[!rownames(raw) %in% target_rows, keep_cell, drop = FALSE]
feat <- as.data.frame(t(feat))
rownames(feat) <- cells

if (nrow(feat) < 5) {
  stop("Only ", nrow(feat), " cells are available. A predictive model needs ",
       "more observations than this", call. = FALSE)
}

step("Modelling data: ", nrow(feat), " cells x ", ncol(feat), " genes")

# --------------------------------------------------------- 2. gene filtering ---
expressed <- colnames(feat)[colSums(feat > opt[["min-fpkm"]], na.rm = TRUE) >=
                              opt[["min-cells"]]]
if (length(expressed) < 2) {
  stop("Fewer than two genes passed the expression filter. Lower --min-fpkm ",
       "or --min-cells", call. = FALSE)
}
note("Genes kept after filtering (FPKM > ", opt[["min-fpkm"]], " in >= ",
     opt[["min-cells"]], " cells): ", length(expressed), " of ", ncol(feat))

feat <- feat[, expressed, drop = FALSE]
feat[!is.finite(as.matrix(feat))] <- 0

# Gene symbols contain characters that R would silently rewrite in a formula,
# so the models run on sanitized column names and every output table carries
# the original symbol back.
gene_names <- colnames(feat)
features <- make.names(gene_names, unique = TRUE)
colnames(feat) <- features
name_map <- data.frame(feature = features, gene_name = gene_names,
                       stringsAsFactors = FALSE)

y <- np_late

# -------------------------------------------------- 3. univariate screening ---
# For a single predictor, the slope test of a linear regression and the test on
# a Pearson correlation are algebraically identical. Everything is therefore
# computed vectorised from the correlation, which gives the same numbers as a
# per-gene regression loop at a small fraction of the cost.
step("Univariate screening of every gene")

X <- as.matrix(feat)
sx <- apply(X, 2, stats::sd)
sy <- stats::sd(y)
n <- length(y)

Xc <- sweep(X, 2, colMeans(X))
yc <- y - mean(y)

r <- as.numeric(t(Xc) %*% yc) / (sx * sy * (n - 1)) * (n - 1)
# The line above reduces to the standard covariance-over-product form. Rewrite
# it explicitly to keep the intent obvious.
r <- as.numeric(t(Xc) %*% yc) / (sqrt(rowSums(Xc^2)) * sqrt(sum(yc^2)))
r[sx == 0 | sy == 0 | !is.finite(r)] <- NA_real_
r[r > 1] <- 1
r[r < -1] <- -1

df_res <- max(n - 2, 1)
t_stat <- r * sqrt(df_res / pmax(1 - r^2, .Machine$double.eps))
p_val <- 2 * stats::pt(-abs(t_stat), df = df_res)
p_val[!is.finite(r)] <- NA_real_

beta <- r * (sy / sx)
beta[sx == 0] <- NA_real_

univariate <- data.frame(
  feature = features,
  gene_name = gene_names,
  pearson_r = r,
  beta = beta,
  r_squared = r^2,
  t_statistic = t_stat,
  p_value = p_val,
  stringsAsFactors = FALSE
)
univariate$p_adj <- bh_fdr(univariate$p_value)
univariate$significant <- !is.na(univariate$p_adj) &
  univariate$p_adj < opt[["max-fdr"]] &
  abs(univariate$pearson_r) >= opt[["min-r"]]
univariate <- univariate[order(-univariate$pearson_r), ]

n_sig <- sum(univariate$significant, na.rm = TRUE)
note("Genes with |r| >= ", opt[["min-r"]], " and adjusted p < ",
     opt[["max-fdr"]], ": ", n_sig)

write_table(univariate,
            file.path(opt[["output-dir"]], "single_gene_regression_results.csv"))

# --------------------------------------------- 4. multivariate elastic net ----
step("Fitting the elastic net regression")

Xf <- X[, is.finite(sx) & sx > 0, drop = FALSE]
kept_features <- colnames(Xf)
if (ncol(Xf) < 2) {
  stop("Fewer than two informative features remain. The multivariate models ",
       "cannot be fitted", call. = FALSE)
}

fit_enet <- function(Xtrain, ytrain, alpha) {
  cvfit <- glmnet::cv.glmnet(Xtrain, ytrain, alpha = alpha,
                             nfolds = min(10L, length(ytrain)))
  glmnet::glmnet(Xtrain, ytrain, alpha = alpha, lambda = cvfit$lambda.min)
}

enet_model <- tryCatch(fit_enet(Xf, y, opt$alpha),
                       error = function(e) {
                         message("  elastic net failed: ", conditionMessage(e))
                         NULL
                       })
enet_pred <- if (!is.null(enet_model)) {
  as.numeric(stats::predict(enet_model, newx = Xf))
} else {
  rep(NA_real_, n)
}

enet_coef <- if (!is.null(enet_model)) {
  cf <- as.numeric(stats::coef(enet_model))[-1]
  cf_df <- data.frame(feature = rownames(stats::coef(enet_model))[-1],
                      coefficient = cf, stringsAsFactors = FALSE)
  cf_df <- merge(name_map, cf_df, by = "feature", all.x = FALSE)
  cf_df[order(-abs(cf_df$coefficient)), , drop = FALSE]
} else {
  data.frame(gene_name = character(), feature = character(),
             coefficient = numeric(), stringsAsFactors = FALSE)
}

if (!is.null(enet_model)) {
  note("Non-zero coefficients: ", sum(enet_coef$coefficient != 0), " of ",
       nrow(enet_coef))
}

# ------------------------------------------------- 5. multivariate random forest
step("Fitting the random forest")

rf_data <- data.frame(NP = y, Xf, check.names = FALSE)
colnames(rf_data)[-1] <- kept_features
mtry_val <- max(1L, min(floor(sqrt(ncol(Xf))), ncol(Xf)))

rf_model <- tryCatch(
  ranger::ranger(NP ~ ., data = rf_data, mtry = mtry_val,
                 num.trees = opt[["n-trees"]],
                 importance = "permutation",
                 min.node.size = opt[["min-node-size"]],
                 seed = opt$seed),
  error = function(e) {
    message("  random forest failed: ", conditionMessage(e))
    NULL
  }
)
rf_pred <- if (!is.null(rf_model)) {
  as.numeric(stats::predict(rf_model, data = Xf)$predictions)
} else {
  rep(NA_real_, n)
}

rf_imp <- if (!is.null(rf_model)) {
  vi <- rf_model$variable.importance
  data.frame(feature = names(vi), importance = as.numeric(vi),
             stringsAsFactors = FALSE)
} else {
  data.frame(feature = character(), importance = numeric(),
             stringsAsFactors = FALSE)
}
if (nrow(rf_imp)) {
  rf_imp <- merge(name_map, rf_imp, by = "feature", all.x = FALSE)
  rf_imp <- rf_imp[order(-rf_imp$importance), ]
}

# --------------------------------------------- 6. leave-one-out performance ---
# Apparent performance is reported for completeness, but with a handful of
# cells it mostly measures how well the model memorised the data. The
# leave-one-out figures are the ones that should be quoted.
do_cv <- tolower(opt$cv) %in% c("loocv", "loo", "cv")
enet_cv <- rep(NA_real_, n)
rf_cv <- rep(NA_real_, n)

if (do_cv && n >= 4) {
  step("Leave-one-out cross validation across ", n, " cells")
  for (i in seq_len(n)) {
    idx <- setdiff(seq_len(n), i)
    if (stats::sd(y[idx]) == 0) next

    enet_i <- tryCatch(fit_enet(Xf[idx, , drop = FALSE], y[idx], opt$alpha),
                       error = function(e) NULL)
    if (!is.null(enet_i)) {
      enet_cv[i] <- as.numeric(stats::predict(enet_i,
                                              newx = Xf[i, , drop = FALSE]))
    }

    rf_i <- tryCatch(
      ranger::ranger(NP ~ ., data = rf_data[idx, , drop = FALSE],
                     mtry = mtry_val, num.trees = opt[["n-trees"]],
                     min.node.size = opt[["min-node-size"]],
                     seed = opt$seed),
      error = function(e) NULL)
    if (!is.null(rf_i)) {
      rf_cv[i] <- as.numeric(stats::predict(rf_i,
                                            data = Xf[i, , drop = FALSE])$predictions)
    }
  }
} else if (do_cv) {
  note("Too few cells for leave-one-out cross validation, skipped")
}

perf <- data.frame(
  Model = c("Elastic net", "Elastic net", "Random forest", "Random forest"),
  Evaluation = c("apparent", "leave-one-out", "apparent", "leave-one-out"),
  RMSE = c(rmse(y, enet_pred), rmse(y, enet_cv),
           rmse(y, rf_pred), rmse(y, rf_cv)),
  R_squared = c(r_squared(y, enet_pred), r_squared(y, enet_cv),
                r_squared(y, rf_pred), r_squared(y, rf_cv)),
  stringsAsFactors = FALSE
)
if (!do_cv) perf$R_squared[perf$Evaluation == "leave-one-out"] <- NA_real_

step("Model performance")
for (i in seq_len(nrow(perf))) {
  note(sprintf("%-14s %-14s RMSE = %8s   R2 = %8s",
               perf$Model[i], perf$Evaluation[i],
               fmt_num(perf$RMSE[i], 4), fmt_num(perf$R_squared[i], 4)))
}

# --------------------------------------------------------------- 7. outputs ---
predictions <- data.frame(
  cell = cells,
  observed_viral_load = y,
  elasticnet_apparent = enet_pred,
  elasticnet_loo = enet_cv,
  randomforest_apparent = rf_pred,
  randomforest_loo = rf_cv,
  stringsAsFactors = FALSE
)
if (!all(is.na(np_mid))) predictions$intermediate_viral_load <- np_mid

write_table(perf, file.path(opt[["output-dir"]], "model_performance.csv"))
write_table(predictions, file.path(opt[["output-dir"]], "predictions.csv"))
write_table(enet_coef, file.path(opt[["output-dir"]], "elasticnet_coefficients.csv"))
write_table(rf_imp, file.path(opt[["output-dir"]], "randomforest_importance.csv"))

if (requireNamespace("openxlsx", quietly = TRUE)) {
  sheets <- list(Univariate_screening = univariate,
                 Model_performance = perf,
                 Predictions = predictions,
                 ElasticNet_coefficients = enet_coef,
                 RandomForest_importance = rf_imp)
  out_xlsx <- file.path(opt[["output-dir"]], "all_results.xlsx")
  openxlsx::write.xlsx(sheets, out_xlsx)
  message("  written: ", out_xlsx)
} else {
  note("Package 'openxlsx' is not installed, the combined workbook was skipped")
}

# --------------------------------------------------------------- 8. summary ---
lines <- c(
  strrep("=", 64),
  "Pre-infection prediction of viral load",
  strrep("=", 64),
  "",
  paste0("Input matrix           : ", basename(opt$fpkm)),
  paste0("Target                 : ", np_late_name),
  paste0("Features               : pre-infection expression of every host gene"),
  "",
  strrep("-", 64),
  "Data",
  strrep("-", 64),
  paste0("  cells modelled       : ", n),
  paste0("  genes in input       : ", nrow(raw) - length(target_rows)),
  paste0("  genes after filter   : ", length(expressed)),
  paste0("  expression floor     : FPKM > ", opt[["min-fpkm"]],
         " in >= ", opt[["min-cells"]], " cells"),
  paste0("  features with variance : ", ncol(Xf)),
  "",
  strrep("-", 64),
  "Univariate screening",
  strrep("-", 64),
  paste0("  test                 : Pearson correlation, equivalently a"),
  paste0("                         single-predictor linear regression"),
  paste0("  multiple testing     : Benjamini-Hochberg"),
  paste0("  top-gene criteria    : |r| >= ", opt[["min-r"]],
         " and adjusted p < ", opt[["max-fdr"]]),
  paste0("  genes passing        : ", n_sig)
)

if (n_sig > 0) {
  lines <- c(lines, "", "  Top 20 genes by Pearson r",
             sprintf("  %-16s %8s %10s %10s", "gene", "r", "p", "FDR"))
  top <- univariate[order(-univariate$pearson_r), ]
  top <- utils::head(top[!is.na(top$pearson_r), ], 20)
  for (i in seq_len(nrow(top))) {
    lines <- c(lines, sprintf("  %-16s %8.3f %10.2e %10.2e",
                              top$gene_name[i], top$pearson_r[i],
                              top$p_value[i], top$p_adj[i]))
  }
}

lines <- c(lines, "",
  strrep("-", 64),
  "Multivariate models",
  strrep("-", 64),
  paste0("  elastic net alpha    : ", opt$alpha),
  paste0("  lambda selection     : cv.glmnet, lambda.min, ",
         min(10L, n), "-fold"),
  paste0("  random forest trees  : ", opt[["n-trees"]]),
  paste0("  random forest mtry   : ", mtry_val),
  paste0("  min node size        : ", opt[["min-node-size"]]),
  paste0("  cross validation     : ", if (do_cv) "leave-one-out" else "none"),
  paste0("  seed                 : ", opt$seed),
  "")
for (i in seq_len(nrow(perf))) {
  lines <- c(lines, sprintf("  %-14s %-14s RMSE = %10s   R2 = %8s",
                            perf$Model[i], perf$Evaluation[i],
                            fmt_num(perf$RMSE[i], 4),
                            fmt_num(perf$R_squared[i], 4)))
}

if (nrow(enet_coef) > 0) {
  nz <- enet_coef[enet_coef$coefficient != 0, ]
  lines <- c(lines, "",
    paste0("  non-zero elastic net coefficients: ", nrow(nz)))
  if (nrow(nz) > 0) {
    nz <- utils::head(nz[order(-abs(nz$coefficient)), ], 15)
    for (i in seq_len(nrow(nz))) {
      lines <- c(lines, sprintf("    %-16s %+10.4f", nz$gene_name[i],
                                nz$coefficient[i]))
    }
  }
}

if (nrow(rf_imp) > 0) {
  lines <- c(lines, "", "  Top 15 random forest permutation importance")
  ti <- utils::head(rf_imp, 15)
  for (i in seq_len(nrow(ti))) {
    lines <- c(lines, sprintf("    %-16s %12.4f", ti$gene_name[i],
                              ti$importance[i]))
  }
}

lines <- c(lines, "",
  strrep("-", 64),
  "Note on interpretation",
  strrep("-", 64),
  "  The apparent figures are in-sample and are optimistic by construction.",
  "  With a small number of cells the leave-one-out figures are the ones that",
  "  should be quoted. A model whose leave-one-out R squared is negative does",
  "  not generalise, which does not invalidate the univariate screening that",
  "  is reported separately above."
)

out_txt <- file.path(opt[["output-dir"]], "prediction_summary.txt")
writeLines(lines, out_txt)
message("  written: ", out_txt)

step("Done")
