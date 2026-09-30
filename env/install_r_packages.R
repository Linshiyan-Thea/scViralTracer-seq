#!/usr/bin/env Rscript
# ==============================================================================
# scViralTracer-seq analysis pipeline -- R dependency installer
# ------------------------------------------------------------------------------
# Run once before using the pipeline.
#
#     Rscript env/install_r_packages.R
#
# Every package below is listed with the script that needs it. Packages marked
# OPTIONAL are only touched on specific input formats or with specific flags, so
# the pipeline runs without them for the common case.
#
#   CRAN
#     readxl        read .xlsx matrices                 all R scripts, OPTIONAL
#     openxlsx      write .xlsx result tables           05 and 06, OPTIONAL
#     data.table    fast grouping of per-cell results   05
#     mclust        two component Gaussian mixture      02
#     glmnet        elastic net regression              06
#     furrr         parallel map over cells             05, OPTIONAL
#     future        parallel backend for furrr          05, OPTIONAL
#
#   Bioconductor
#     DESeq2        differential expression             03
#
# Tested with R 4.3 and Bioconductor 3.18.
# ==============================================================================

cran_pkgs <- c("readxl", "openxlsx", "data.table", "mclust",
                "glmnet", "furrr", "future")
bioc_pkgs <- c("DESeq2")

# ----------------------------------------------------------------- CRAN part ---
missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
} else {
  message("All CRAN packages are already installed")
}

# --------------------------------------------------------- Bioconductor part ---
missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]
if (length(missing_bioc) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)
} else {
  message("All Bioconductor packages are already installed")
}

# -------------------------------------------------------------------- verify ---
message("\nVersion check")
all_pkgs <- c(cran_pkgs, bioc_pkgs)
for (pkg in all_pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  ver <- if (ok) as.character(utils::packageVersion(pkg)) else "MISSING"
  message(sprintf("  %-12s %s", pkg, ver))
}

message("\nR version: ", R.version.string)
message("Done")
