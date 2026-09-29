# scViralTracer-seq analysis pipeline


scViralTracer-seq withdraws a defined volume of cytoplasm from one living cell with an
electric-field-driven nanopipette, so the same cell can be sequenced before infection and
again at several hours post infection (hpi). The repository holds every analysis step that
turns those libraries into the numbers reported in the paper: read processing and library
QC, separation of infected from bystander cells, differential expression, diffusion
pseudotime, within-cell coupling of host genes to viral load, prediction of viral load from
the pre-infection transcriptome, and the comparison against bulk and droplet-based data.

---

## What is included

All analysis scripts, their parameter defaults, the shared helper library and
the environment specifications.

---

## Requirements

### Command line tools

Used by `01_read_processing_and_qc/`.

| Tool | Purpose | Version used |
|---|---|---|
| FastQC / MultiQC | raw read QC | 0.12.1 |
| Trim Galore | adapter and quality trimming | 0.6.10 |
| STAR | spliced alignment | 2.7.11b |
| samtools | BAM handling and down-sampling | 1.19 |
| Subread (`featureCounts`) | gene-level quantification | 2.0.6 |

### R

Used by steps 02, 03, 05, 06 and 07. Tested with R 4.3 and Bioconductor 3.18.

```bash
Rscript env/install_r_packages.R
```

This installs `readxl`, `openxlsx`, `data.table`, `mclust`, `glmnet`, `ranger`, `furrr`,
`future` from CRAN and `DESeq2` from Bioconductor, then prints the resolved versions.
`readxl`, `openxlsx`, `furrr` and `future` are optional and only touched for `.xlsx`
input, `.xlsx` output, or parallel execution.

### Python

Used by step 04. Tested with Python 3.10.

```bash
conda env create -f env/environment.yml && conda activate scviraltracer
# or, without conda
pip install -r env/requirements.txt
```

`scanpy`, `numpy`, `pandas`, `scipy` and, for `.xlsx` input, `openpyxl`.

---

## Directory layout

```
.
├── README.md
├── .gitignore
├── R/
│   └── common.R                        shared CLI parsing, matrix IO, statistics
├── env/
│   ├── environment.yml                 conda environment for Python and the CLI tools
│   ├── requirements.txt                pip equivalent for the Python part
│   └── install_r_packages.R            one-shot R dependency installer
├── 01_read_processing_and_qc/
│   ├── 01_align_and_quantify.sh        trimming, STAR, featureCounts, down-sampling
│   └── 02_library_statistics.sh        consolidated per-library QC table
├── 02_infected_vs_bystander/
│   └── classify_infected_bystander.R   two-component mixture on total viral load
├── 03_differential_expression/
│   └── de_deseq2.R                     DESeq2, optionally with a per-cell blocking factor
├── 04_pseudotime/
│   └── scanpy_dpt.py                   diffusion pseudotime and its validation
├── 05_within_cell_gene_viral_coupling/
│   └── per_cell_gene_viral_coupling.R  per-cell Pearson r of every gene against NP
├── 06_preinfection_prediction/
│   └── predict_viral_load.R            regression of viral load on the pre-infection state
└── 07_conventional_data_comparison/
    ├── bulk_gene_viral_correlation.R   bulk time course correlation
    └── tenx_correlation_summary.R      droplet-based data summary
```

---

## Input conventions

Read these three conventions before running anything, because every script depends on them.

### 1. Sample naming

```
<CellID>.<TimePoint>
```

`<CellID>` identifies one physical cell and `<TimePoint>` is an integer number of hours
post infection. For example a cell sampled three times gives

```
WSN01.0    WSN01.3    WSN01.6
```

The pre-infection library is `.0`. Scripts 02, 04 and 05 recover the longitudinal pairing
by parsing this name, so a library whose name does not follow the convention is rejected
with an explicit error rather than silently dropped.

Bulk and droplet-based libraries are not longitudinally paired, and those scripts parse
their own naming. Step 07 documents the pattern it expects.

### 2. Matrix orientation

Genes are **rows**, samples are **columns**. The first column holds the gene names.
`.xlsx`, `.csv` and `.tsv` are accepted interchangeably by every R and Python script.

Two kinds of matrix are used and they are not interchangeable.

| Matrix | Content | Consumed by |
|---|---|---|
| count | raw integer fragment counts | 03 (DESeq2), 04 (trajectory) |
| FPKM | length and depth normalised values | 02, 05, 06, 07 |

Step 03 rejects non-integer input on purpose, because feeding normalised values to DESeq2
invalidates its variance model.

### 3. Viral gene names

Viral rows are recognised by name, not by position. The defaults match the annotation used
here, in which the eight IAV WSN segments are prefixed with the strain name.

| Default pattern | Scripts |
|---|---|
| `^WSN_` | 02 |
| `_NP$` | 05 (the viral load proxy) |
| `^WSN` | 04 (the rows zeroed before infection) |
| `WSN_NP` | 07 (bulk) |

Override any of them with the `--viral-genes`, `--viral-gene` or `--viral-pattern` flag of
the corresponding script.

---

## Running the pipeline

The steps are independent once the matrices exist, so run only what you need. The commands
below assume the repository root is the working directory and that inputs live under
`./data`.

### Step 01 — read processing and library QC

```bash
bash 01_read_processing_and_qc/01_align_and_quantify.sh \
     --fastq-dir  ./data/fastq \
     --work-dir   ./work \
     --star-index ./reference/star_index \
     --gtf        ./reference/combined_host_virus.gtf \
     --threads 12

bash 01_read_processing_and_qc/02_library_statistics.sh \
     --work-dir ./work
```

The reference is a combined host and viral annotation. The host part is the GENCODE primary
assembly annotation and the viral part is the eight IAV A/WSN/1933 segments appended as
extra contigs, so that viral reads are assigned to viral genes instead of being counted as
unmapped. A third reference, the canine genome, is used only for the stock-purity control
reported in the supplementary material.

Down-sampling to a common depth is the last stage of step 01 and it matters for this assay.
Libraries from a single cell differ several-fold in size, and a per-cell comparison of
expression would otherwise be confounded by that difference. Every library is therefore
reduced to the same number of properly paired fragments and quantified a second time.

| Flag | Default | Meaning |
|---|---|---|
| `--fastq-dir` | required | directory of `<SAMPLE>_1.fq.gz` / `_2.fq.gz` pairs |
| `--work-dir` | required | all outputs are written below this directory |
| `--star-index` | required | directory of the STAR genome index |
| `--gtf` | required | combined host and viral annotation |
| `--threads` | 12 | threads for trimming, alignment and quantification |
| `--tmp-dir` | `<work-dir>/tmp` | scratch space for `featureCounts` |
| `--target-fragments` | 1000000 | common down-sampling depth |
| `--seed` | 123 | down-sampling seed, fixed so the result is reproducible |
| `--skip-fastqc` | off | skip the raw read QC stage |
| `--skip-trim` | off | align the raw reads directly |
| `--mt-contig` | `chrM` | mitochondrial contig name, step 02 only |
| `--mapq` | 255 | minimum MAPQ for a uniquely mapped read, step 02 only |

Output of step 02 is `<work-dir>/02_results/full_stats_complete.tsv`, one row per library,
carrying the uniquely mapped fraction, the assigned fragment count at the common depth, the
mitochondrial fraction, the assigned fragment count of the full library and the detected
gene number.

### Step 02 — infected versus bystander cells

```bash
Rscript 02_infected_vs_bystander/classify_infected_bystander.R \
        --fpkm ./data/merged_fpkm.xlsx \
        --output-dir ./results/02_classification
```

Residual viral reads from the inoculum are present in every library, so a zero/non-zero
cut-off does not separate the populations. The log-transformed total viral expression is
bimodal instead, and a two-component Gaussian mixture with equal variance is fitted at the
classification time point. The component means are seeded with the background estimated at
the pre-infection time point and with the upper range of the observed values. The component
with the lower fitted mean is reported as bystander and the other as infected, which is
decided by the fit rather than assumed.

| Flag | Default | Meaning |
|---|---|---|
| `--fpkm` | required | gene-by-cell FPKM matrix |
| `--output-dir` | required | where the tables are written |
| `--viral-genes` | empty | explicit viral row names, overrides the pattern |
| `--viral-pattern` | `^WSN_` | regular expression matching viral rows |
| `--background-time` | smallest time point | time point defining the background |
| `--classify-time` | largest time point | time point at which cells are classified |
| `--model-name` | `E` | `mclust` univariate model, `E` equal or `V` varying variance |

### Step 03 — differential expression

```bash
Rscript 03_differential_expression/de_deseq2.R \
        --counts  ./data/merged_count.xlsx \
        --coldata ./data/group_all.xlsx \
        --factor group --alt High --ref Low \
        --blocking cell \
        --output-dir ./results/03_de
```

`--coldata` is a sample sheet with at least a `sample` column matching the count matrix
column names and a column named after `--factor`. `--ref` fixes the reference level so the
sign of the log2 fold change is unambiguous, and `--alt` is the numerator of the contrast.

**Significance thresholds used throughout the manuscript**

```
|log2 fold change| >= 1    and    adjusted p < 0.05
```

These are the defaults of `--lfc-cutoff` and `--padj-cutoff`. Adjusted p-values are the
Benjamini-Hochberg values that DESeq2 returns in its `padj` column.

`--blocking` deserves attention for this design. The post-infection time points are
compared against the pre-infection state **of the same cells**, so each cell contributes
libraries to both arms of the contrast. Passing the cell identifier as a blocking factor
gives the design `~ cell + group`, which removes the between-cell variation that would
otherwise inflate the dispersion. Omit the flag for an unpaired comparison.

| Flag | Default | Meaning |
|---|---|---|
| `--counts` | required | integer count matrix |
| `--coldata` | required | sample sheet |
| `--factor` | required | column of `--coldata` holding the group labels |
| `--alt` | required | numerator level of the contrast |
| `--ref` | required | reference level of the contrast |
| `--output-dir` | required | where the tables are written |
| `--lfc-cutoff` | 1 | absolute log2 fold change threshold |
| `--padj-cutoff` | 0.05 | adjusted p-value threshold |
| `--blocking` | empty | column used as a blocking factor for a paired design |
| `--min-cells-expressed` | `all` | keep a gene detected in all samples, or in at least N |
| `--min-total-count` | 20 | keep a gene whose total count exceeds this value |
| `--independent-filtering` | `TRUE` | DESeq2 independent filtering |

### Step 04 — diffusion pseudotime

```bash
python 04_pseudotime/scanpy_dpt.py \
       --matrix-dir ./data \
       --output-dir ./results/04_pseudotime \
       --time-points 0,3,6
```

`--matrix-dir` must hold one count and one FPKM matrix per time point, named
`merged_count_<T>.<ext>` and `merged_fpkm_<T>.<ext>`.

The ordering is computed from the count matrix only, as counts per million followed by
`log1p`. FPKM values are carried alongside for the gene-level correlations, because they
are the unit in which the rest of the pipeline reports expression.

Root placement is the one place where sampling time is used. The root is restricted to the
libraries taken before infection and, within those, the library with the smallest first
principal component is chosen so that the choice is deterministic. The neighbour graph, the
diffusion map and the pseudotime values themselves are computed without reference to the
time labels, which is what makes the subsequent agreement between pseudotime and recorded
time a genuine test rather than a tautology.

Viral rows are set to zero before infection. Any viral read in a pre-infection library comes
from the inoculum and not from transcription by that cell, so keeping it would place a cell
on the ordering according to a signal it did not produce. Host genes keep their true basal
value, because that basal level is exactly what the assay is designed to record.

| Flag | Default | Meaning |
|---|---|---|
| `--matrix-dir` | required | directory of the per-time-point matrices |
| `--output-dir` | required | where the tables are written |
| `--time-points` | `0,3,6` | integer hours post infection |
| `--count-prefix`, `--fpkm-prefix` | `merged_count`, `merged_fpkm` | file name prefixes |
| `--sample-pattern` | `.*` | keep only columns matching this expression |
| `--n-pcs` | 15 | principal components |
| `--n-neighbors` | 10 | neighbours in the k-nearest-neighbour graph |
| `--diffmap-comps` | 15 | diffusion map components |
| `--scale-max-value` | 10 | clipping value applied when scaling |
| `--filter-min-count` | 1 | drop a gene whose total count is below this |
| `--filter-min-samples` | 3 | keep a gene detected in at least this many libraries |
| `--target-sum` | 1e6 | library size normalisation target |
| `--random-state` | 0 | seed for every stochastic step |
| `--root-cell` | empty | name the root library explicitly |
| `--root-time` | earliest time point | time point the root is drawn from |
| `--viral-pattern` | `^WSN` | rows zeroed before infection |
| `--genes` | ISG and viral panel | genes correlated with pseudotime |
| `--fpkm-offset` | 0.1 | offset added before the log2 transform of FPKM |

The gene filter that actually decides the gene set is the "detected in at least
`--filter-min-samples` libraries" rule, because the preceding `min_counts` filter is applied
after the log transform.

### Step 05 — within-cell coupling of host genes to viral load

```bash
Rscript 05_within_cell_gene_viral_coupling/per_cell_gene_viral_coupling.R \
        --fpkm ./data/merged_fpkm.xlsx \
        --output-dir ./results/05_coupling
```

For every cell, and using only that cell's own time points, the Pearson correlation between
each host gene and the viral NP transcript is computed. The per-cell coefficients are then
summarised per gene and tested for consistency across cells with a **one-sided one-sample
t-test against zero**, followed by Benjamini-Hochberg correction. One-sided is correct here
because the question is whether a gene rises together with the virus, and a gene that moves
in the opposite direction in every cell is a different finding rather than a weaker one.

| Flag | Default | Meaning |
|---|---|---|
| `--fpkm` | required | gene-by-cell FPKM matrix |
| `--output-dir` | required | where the tables are written |
| `--viral-gene` | from pattern | row name of the viral load proxy |
| `--viral-pattern` | `_NP$` | regular expression locating the NP row |
| `--min-mean-corr` | 0.7 | minimum mean correlation across cells |
| `--max-fdr` | 0.01 | maximum Benjamini-Hochberg q-value |
| `--min-cells` | 3 | minimum number of cells with a usable coefficient |
| `--min-timepoints` | 2 | minimum time points a cell must contribute |
| `--workers` | 1 | parallel workers. One means sequential execution |

A gene is retained when it passes all three of `mean r > 0.7`, `FDR < 0.01` and at least
three informative cells.

### Step 06 — predicting viral load from the pre-infection transcriptome

```bash
Rscript 06_preinfection_prediction/predict_viral_load.R \
        --fpkm ./data/merged_fpkm_preinfection.xlsx \
        --np-row-6h WSN_NP_6hpi \
        --np-row-3h WSN_NP_3hpi \
        --output-dir ./results/06_prediction
```

The feature matrix holds the pre-infection expression of every host gene and the target is
the viral NP expression of the same cell at the late time point. Three models are reported.

1. A univariate screen of every gene. For a single predictor the slope test of a linear
   regression and the test on a Pearson correlation are algebraically identical, so the
   screen is computed vectorised from the correlation.
2. An elastic net regression with mixing parameter `alpha = 0.5`.
3. A random forest with permutation importance.

Both multivariate models are reported with their apparent fit **and** with leave-one-out
cross-validation. With a handful of cells the apparent fit is strongly optimistic, and the
cross-validated number is the one to quote.

| Flag | Default | Meaning |
|---|---|---|
| `--fpkm` | required | pre-infection gene-by-cell FPKM matrix |
| `--np-row-6h` | required | row holding the late viral load per cell |
| `--np-row-3h` | empty | row holding an intermediate viral load, retained not modelled |
| `--output-dir` | required | where the tables are written |
| `--min-fpkm` | 1 | expression floor for a gene to count as expressed |
| `--min-cells` | 3 | a gene is kept when this many cells are above the floor |
| `--min-r` | 0.65 | absolute Pearson r flagging the top-correlated genes |
| `--max-fdr` | 0.05 | adjusted p-value used together with `--min-r` |
| `--alpha` | 0.5 | elastic net mixing parameter |
| `--n-trees` | 500 | random forest trees |
| `--min-node-size` | 1 | random forest leaf size |
| `--cv` | `loocv` | `loocv` or `none` |
| `--seed` | 1 | seed for the random forest |

The gene set carried into the enrichment analysis is the one flagged by `--min-r` and
`--max-fdr`.

### Step 07 — comparison with conventional data

Bulk time course:

```bash
Rscript 07_conventional_data_comparison/bulk_gene_viral_correlation.R \
        --fpkm ./data/bulk_fpkm.xlsx \
        --query-genes IFRD2 \
        --output-dir ./results/07_bulk
```

Sample names must contain the hours post infection. The default pattern `(?i)([0-9]+)\s*hpi`
reads it out of names such as `A0hpi`, `A3hpi`, `A4hpi`, `A5hpi` and `A6hpi`. Pass
`--time-values` if your naming differs. Because a bulk time course has few points, the
retention threshold is deliberately strict at `r > 0.95`.

Droplet-based data:

```bash
Rscript 07_conventional_data_comparison/tenx_correlation_summary.R \
        --correlations ./data/all_genes_correlation_with_NP.csv \
        --query-genes IFRD2 \
        --output-dir ./results/07_tenx
```

Column names of the input table are detected automatically, so the usual spellings of the
gene, coefficient and p-value columns all work. Alternatively compute the correlations from
a matrix with `--matrix` and `--viral-gene`. Retention is `|r| > 0.5` with `FDR < 0.05`.

The point of this step is a negative control. `--query-genes` reports the correlation of a
specific gene on the conventional platform next to the same gene's within-cell correlation
from step 05.

---

## Statistical thresholds at a glance

| Analysis | Script | Threshold |
|---|---|---|
| Differential expression | 03 | \|log2 fold change\| >= 1 **and** adjusted *p* < 0.05 |
| Within-cell coupling | 05 | mean *r* > 0.7 **and** FDR < 0.01 **and** >= 3 informative cells |
| Pre-infection prediction | 06 | \|*r*\| >= 0.65 **and** adjusted *p* < 0.05 |
| Bulk time course | 07 | *r* > 0.95 **and** FDR < 0.05 |
| Droplet-based data | 07 | \|*r*\| > 0.5 **and** FDR < 0.05 |

Multiple testing is corrected with the Benjamini-Hochberg procedure everywhere, and the
corrected value is always reported in a column named `fdr`, `p_adj` or `padj`.

Correlations are Pearson unless a Spearman value is reported alongside. All correlation
tests are two-sided with the single exception of step 05, where the across-cell test of the
per-cell coefficients is one-sided against zero for the reason given above.

---

## Output files

Every script writes plain tables and a text summary of the settings it ran with. Nothing is
written outside `--output-dir`.

| Script | Files |
|---|---|
| 02 | `cell_classification.csv`, `classification_summary.txt` |
| 03 | `DESeq2_results.csv`, `DEG_up.csv`, `DEG_down.csv`, `DEG_summary.txt` |
| 04 | `dpt_data.csv`, `per_cell_trajectory.csv`, `dpt_validation.csv`, `pseudotime_by_timepoint.csv`, `gene_pseudotime_correlation.csv`, `pseudotime_summary.txt`, `anndata.h5ad` |
| 05 | `per_cell_gene_correlations.csv`, `gene_correlation_summary.csv`, `significant_genes.csv`, `significant_gene_expression.xlsx`, `coupling_summary.txt` |
| 06 | `single_gene_regression_results.csv`, `model_performance.csv`, `predictions.csv`, `elasticnet_coefficients.csv`, `randomforest_importance.csv`, `all_results.xlsx`, `prediction_summary.txt` |
| 07 | `bulk_all_gene_correlations.csv`, `bulk_significant_genes.csv`, `bulk_query_genes.csv`, `bulk_correlation_summary.txt`, `tenx_all_gene_correlations.csv`, `tenx_significant_genes.csv`, `tenx_top_correlated_genes.csv`, `tenx_query_genes.csv`, `tenx_correlation_report.txt` |

Step 04 additionally writes the UMAP and diffusion coordinates into `dpt_data.csv`. Those
columns are the source data of the trajectory figures, which is why they are exported even
though the script produces no figure.

---

## Reproducibility

* **No hard-coded paths.** Every input and output location is a command line argument. The
  R scripts locate the shared helper library relative to their own file, so they run from
  any working directory.
* **Fixed seeds.** Down-sampling uses `--seed`, the random forest uses `--seed`, and every
  stochastic step of the pseudotime analysis uses `--random-state`.
* **Recorded settings.** Each `*_summary.txt` and `*_report.txt` file restates the
  thresholds and parameters that produced the tables next to it.
* **Explicit failures.** A sample name that does not follow the convention, a non-integer
  count matrix, or a missing viral row stops the script with a message naming the offending
  value. Nothing is silently coerced.

---

## Notes on the released version

These scripts are a cleaned release of the analysis code used for the manuscript. Four
problems in the internal version were fixed on the way, and the fixes change results, so
they are recorded here.

1. **Up and down labels in step 03.** The internal version assigned the direction from
   `log2FoldChange > 1`, which mislabelled a gene with a log2 fold change of exactly 1 as
   down-regulated. The direction now comes from the sign of the estimate and the magnitude
   is used only for the significance call.
2. **Paired design in step 03.** The longitudinal comparison is against the same cells
   before infection, which the internal version ran as an unpaired contrast. `--blocking`
   now expresses that pairing.
3. **Pre-infection zeroing in step 04.** The internal version zeroed the pre-infection
   values of every exported gene, including host genes, which discarded the basal
   expression that the assay exists to measure. Zeroing is now restricted to viral rows.
4. **Undefined helper in step 06.** The internal version called a root-mean-squared-error
   function from a package it deliberately did not load. The metric is now defined locally.

---

## Citation

If you use this code, please cite the manuscript. The diffusion pseudotime method is
described in

> Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. Diffusion pseudotime robustly
> reconstructs lineage branching. *Nat Methods*. 2016;13(10):845-848. doi:10.1038/nmeth.3971

## Licence

Released under the MIT Licence. See `LICENSE`.
