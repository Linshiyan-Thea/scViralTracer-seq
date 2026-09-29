#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
scViralTracer-seq analysis pipeline
Diffusion pseudotime ordering of longitudinally sampled single cells
==============================================================================

What it does
    Because scViralTracer-seq samples the same living cell more than once, the
    true temporal order of the libraries is known. Diffusion pseudotime is
    computed without using that order, so agreement between the two is a
    genuine test of whether the assay resolves infection as a coherent
    transcriptional progression.

    1. Load the per-time-point count and FPKM matrices and merge them into one
       AnnData object.
    2. Library-size normalise to counts per million and apply log1p.
    3. Filter genes, scale, run PCA and build the k-nearest-neighbour graph.
    4. Compute UMAP and the diffusion map.
    5. Place the root cell and compute diffusion pseudotime.
    6. Correlate pseudotime with the recorded sampling time.
    7. Correlate the expression of the eight viral genomic segments with
       pseudotime, with Benjamini-Hochberg correction across genes.

    No figure is produced. The UMAP coordinates are written out as source data
    so that plotting can be done separately.

Input
    A directory holding one count matrix and one FPKM matrix per time point,
    named

        merged_count_<TIME>.<ext>      merged_fpkm_<TIME>.<ext>

    where <TIME> is an integer number of hours, for example 0, 3 and 6, and
    <ext> is .xlsx, .csv or .tsv. Genes are rows, samples are columns.

    Sample names must follow the <CellID>.<TimePoint> convention used
    everywhere else in this pipeline, for example WSN01.0, WSN01.3, WSN01.6.
    Every time point of one cell must share the same <CellID>.

Output (written to --output-dir)
    dpt_data.csv                     one row per library, with time, cell,
                                     pseudotime, UMAP coordinates and the
                                     expression of the requested genes
    per_cell_trajectory.csv          the three libraries of each cell ordered
                                     by sampling time, ready to be connected
    dpt_validation.csv               pseudotime against recorded time
    gene_pseudotime_correlation.csv  per-gene correlation with pseudotime
    pseudotime_summary.txt           parameters and headline statistics
    anndata.h5ad                     the processed object

Reference
    Haghverdi L, Buettner M, Wolf FA, Buettner F, Theis FJ.
    Diffusion pseudotime robustly reconstructs lineage branching.
    Nat Methods. 2016;13(10):845-848. doi:10.1038/nmeth.3971

Requirements
    scanpy, numpy, pandas, scipy, openpyxl (only for .xlsx input)

Usage
    python scanpy_dpt.py --matrix-dir ./matrix --output-dir ./results_pseudotime

    python scanpy_dpt.py \\
        --matrix-dir ./matrix \\
        --output-dir ./results_pseudotime \\
        --time-points 0,3,6 \\
        --n-pcs 15 --n-neighbors 10 \\
        --genes WSN_PB2,WSN_PB1,WSN_PA,WSN_HA,WSN_NP,WSN_NA,WSN_M,WSN_NS
"""

import argparse
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

SUPPORTED_EXT = (".xlsx", ".csv", ".tsv", ".txt")


# --------------------------------------------------------------------- input --
def parse_args():
    p = argparse.ArgumentParser(
        description="Diffusion pseudotime analysis of longitudinally sampled "
                    "single cells.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--matrix-dir", required=True,
                   help="Directory holding merged_count_<T>.* and merged_fpkm_<T>.*")
    p.add_argument("--output-dir", required=True,
                   help="Where the result tables are written")

    p.add_argument("--time-points", default="0,3,6",
                   help="Comma separated integer hours post infection")
    p.add_argument("--count-prefix", default="merged_count",
                   help="File name prefix of the count matrices")
    p.add_argument("--fpkm-prefix", default="merged_fpkm",
                   help="File name prefix of the FPKM matrices")
    p.add_argument("--sample-pattern", default=".*",
                   help="Regular expression. Only columns matching it fully are kept")

    p.add_argument("--n-pcs", type=int, default=15,
                   help="Number of principal components")
    p.add_argument("--n-neighbors", type=int, default=10,
                   help="Neighbours for the k-nearest-neighbour graph")
    p.add_argument("--diffmap-comps", type=int, default=15,
                   help="Diffusion map components")
    p.add_argument("--scale-max-value", type=float, default=10.0,
                   help="Clipping value applied when scaling")
    p.add_argument("--filter-min-count", type=int, default=1,
                   help="Drop a gene whose total count is below this value")
    p.add_argument("--filter-min-samples", type=int, default=3,
                   help="Keep a gene only when it is non-zero in at least this "
                        "many libraries")
    p.add_argument("--target-sum", type=float, default=1e6,
                   help="Library size normalisation target")
    p.add_argument("--random-state", type=int, default=0,
                   help="Seed for the stochastic steps, so that UMAP and the "
                        "neighbour graph are reproducible")

    p.add_argument("--root-cell", default="",
                   help="Name of the root library. When empty the root is "
                        "chosen automatically among the libraries sampled at "
                        "--root-time")
    p.add_argument("--root-time", type=float, default=np.nan,
                   help="Time point from which the root is taken. Defaults to "
                        "the earliest time point present")

    p.add_argument("--viral-pattern", default="^WSN",
                   help="Row names matching this pattern are viral. Viral "
                        "expression is set to zero before infection because "
                        "no viral transcript can originate from the cell at "
                        "that point. Host genes keep their true basal value")
    p.add_argument("--genes", default="WSN_PB2,WSN_PB1,WSN_PA,WSN_HA,WSN_NP,"
                                      "WSN_NA,WSN_M,WSN_NS",
                   help="Comma separated genes to correlate with pseudotime. "
                        "The default is the eight influenza virus genomic "
                        "segments")
    p.add_argument("--fpkm-offset", type=float, default=0.1,
                   help="Offset added before the log2 transform of FPKM")

    return p.parse_args()


def split_csv(value):
    return [s.strip() for s in str(value).split(",") if s.strip()]


def find_matrix(directory, prefix, time_point):
    """Locate <prefix>_<time_point>.<ext> whatever the extension is."""
    for ext in SUPPORTED_EXT:
        cand = os.path.join(directory, "%s_%s%s" % (prefix, time_point, ext))
        if os.path.isfile(cand):
            return cand
    raise FileNotFoundError(
        "No matrix found for time point %s in %s. Looked for %s_%s followed by "
        "one of %s" % (time_point, directory, prefix, time_point,
                       ", ".join(SUPPORTED_EXT)))


def read_table(path):
    """Read a gene-by-sample matrix from .xlsx, .csv or .tsv."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".xlsx":
        df = pd.read_excel(path, index_col=0)
    elif ext == ".csv":
        df = pd.read_csv(path, index_col=0)
    else:
        df = pd.read_csv(path, index_col=0, sep="\t")
    df.index = df.index.astype(str)
    df.columns = [str(c) for c in df.columns]
    return df


def parse_sample_names(samples):
    """Split <CellID>.<TimePoint> into its two parts."""
    pat = re.compile(r"^(.*)\.([0-9]+)$")
    rows = []
    bad = []
    for s in samples:
        m = pat.match(s)
        if not m:
            bad.append(s)
            continue
        rows.append({"sample": s, "cell_id": m.group(1),
                     "time": int(m.group(2))})
    if bad:
        raise ValueError(
            "These sample names do not follow the <CellID>.<TimePoint> "
            "convention: %s" % ", ".join(bad[:5]))
    return pd.DataFrame(rows).set_index("sample")


def load_data(args):
    """Merge the per-time-point matrices into count and FPKM frames."""
    time_points = [int(t) for t in split_csv(args.time_points)]
    sample_re = re.compile(args.sample_pattern)

    counts, fpkms = {}, {}
    for tp in time_points:
        c_path = find_matrix(args.matrix_dir, args.count_prefix, tp)
        f_path = find_matrix(args.matrix_dir, args.fpkm_prefix, tp)
        print("  time %dh  counts: %s" % (tp, os.path.basename(c_path)))
        print("            FPKM : %s" % os.path.basename(f_path))

        c_df = read_table(c_path)
        f_df = read_table(f_path)
        counts[tp] = c_df[[c for c in c_df.columns if sample_re.fullmatch(c)]]
        fpkms[tp] = f_df[[c for c in f_df.columns if sample_re.fullmatch(c)]]

    # Keep only genes measured at every time point, so that the longitudinal
    # comparisons are made on a common gene set.
    common = set(counts[time_points[0]].index)
    for tp in time_points[1:]:
        common &= set(counts[tp].index)
    common = sorted(common)
    print("  genes shared by all time points: %d" % len(common))

    count_df = pd.concat([counts[tp].loc[common] for tp in time_points], axis=1)
    fpkm_df = pd.concat([fpkms[tp].loc[common] for tp in time_points], axis=1)

    meta = parse_sample_names(count_df.columns.tolist())
    observed = sorted(meta["time"].unique().tolist())
    if observed != sorted(time_points):
        print("  WARNING: time points parsed from sample names (%s) differ "
              "from --time-points (%s)" % (observed, sorted(time_points)))

    return count_df, fpkm_df, meta, time_points


# ------------------------------------------------------------------ statistics --
def bh_fdr(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, dtype=float)
    m = p.size
    if m == 0:
        return np.array([])
    finite = np.isfinite(p)
    q = np.full(m, np.nan)
    if not finite.any():
        return q
    idx = np.where(finite)[0]
    sub = p[idx]
    n = sub.size
    order = np.argsort(sub)
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(1, n + 1)
    vals = sub * n / ranks
    vals_sorted = np.minimum.accumulate(vals[order][::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = vals_sorted
    q[idx] = np.clip(out, 0.0, 1.0)
    return q


# --------------------------------------------------------------------- pipeline --
def run_dpt(count_df, meta, args):
    """Normalise, embed and order the libraries."""
    import scanpy as sc
    sc.settings.verbosity = 1

    # AnnData holds cells in rows and genes in columns.
    adata = sc.AnnData(X=count_df.T.values.astype(np.float64))
    adata.obs_names = list(count_df.columns)
    adata.var_names = list(count_df.index)
    adata.obs["time"] = meta.loc[adata.obs_names, "time"].values
    adata.obs["cell_id"] = meta.loc[adata.obs_names, "cell_id"].values
    print("  input: %d libraries x %d genes" % (adata.n_obs, adata.n_vars))

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)

    sc.pp.filter_genes(adata, min_counts=args.filter_min_count)
    from scipy.sparse import issparse
    X = adata.X
    if issparse(X):
        detected = np.asarray((X > 0).sum(axis=0)).flatten()
    else:
        detected = np.asarray((X > 0).sum(axis=0)).flatten()
    adata = adata[:, detected >= args.filter_min_samples].copy()
    print("  genes after filtering: %d" % adata.n_vars)

    if adata.n_obs < 4:
        sys.exit("ERROR: at least four libraries are needed to build a "
                 "neighbour graph")

    # Highly variable gene selection is deliberately skipped. The data set is
    # small enough that every filtered gene carries information, and selecting
    # a subset would add a tunable step whose effect on the ordering is hard
    # to assess with so few cells.
    sc.pp.scale(adata, max_value=args.scale_max_value)
    n_comps = int(min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=n_comps, random_state=args.random_state)
    print("  PCA components: %d" % adata.obsm["X_pca"].shape[1])

    n_neighbors = int(min(args.n_neighbours, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors,
                    n_pcs=adata.obsm["X_pca"].shape[1],
                    random_state=args.random_state)
    print("  neighbour graph: k = %d" % n_neighbors)

    sc.tl.umap(adata, random_state=args.random_state)

    n_diff = int(min(args.diffmap_comps, adata.n_obs - 2))
    n_diff = max(n_diff, 3)
    sc.tl.diffmap(adata, n_comps=n_diff, random_state=args.random_state)
    print("  diffusion map components: %d" % (adata.obsm["X_diffmap"].shape[1] - 1))

    adata.uns["iroot"] = int(pick_root(adata, args))
    print("  root library: %s (time = %dh)"
          % (adata.obs_names[adata.uns["iroot"]],
             int(adata.obs["time"].iloc[adata.uns["iroot"]])))

    sc.tl.dpt(adata, n_branchings=0)
    print("  pseudotime range: %.4f to %.4f"
          % (adata.obs["dpt_pseudotime"].min(),
             adata.obs["dpt_pseudotime"].max()))
    return adata, n_neighbors


def pick_root(adata, args):
    """Return the index of the root library.

    The root is fixed by the recorded sampling time rather than by the
    embedding, so that pseudotime is anchored on experimental ground truth.
    Within the earliest time point the library with the smallest first
    principal component is used, which gives a deterministic choice.
    """
    if args.root_cell:
        if args.root_cell not in set(adata.obs_names):
            sys.exit("ERROR: --root-cell %s is not among the libraries"
                     % args.root_cell)
        return adata.obs_names.get_loc(args.root_cell)

    root_time = args.root_time
    if not np.isfinite(root_time):
        root_time = adata.obs["time"].min()
    mask = (adata.obs["time"] == root_time).values
    if not mask.any():
        sys.exit("ERROR: no library was sampled at time %s" % root_time)
    pc1 = adata.obsm["X_pca"][:, 0]
    candidates = np.where(mask)[0]
    return candidates[int(np.argmin(pc1[candidates]))]


# -------------------------------------------------------------------- analyses --
def viral_zeroed(values, is_viral, time_values, earliest_time):
    """Set viral expression to zero before infection.

    Any viral read in a pre-infection library comes from the inoculum rather
    than from transcription by that cell, so keeping it would place a cell on
    the trajectory according to a signal it did not produce. Host genes are
    left untouched, because their basal level is precisely the pre-infection
    state that the assay is designed to record.
    """
    out = values.copy()
    if is_viral:
        out[np.asarray(time_values) == earliest_time] = 0.0
    return out


def correlate_genes(adata, fpkm_df, args):
    """Correlate gene expression with pseudotime."""
    dpt = adata.obs["dpt_pseudotime"].values.astype(float)
    time_values = adata.obs["time"].values
    earliest = min(time_values)
    viral_re = re.compile(args.viral_pattern)
    offset = args.fpkm_offset

    genes = split_csv(args.genes)
    rows = []
    for gene in genes:
        if gene not in fpkm_df.index:
            print("  WARNING: %s is absent from the FPKM matrix, skipped" % gene)
            continue
        vals = fpkm_df.loc[gene].reindex(adata.obs_names).values.astype(float)
        is_viral = bool(viral_re.search(gene))
        vals = viral_zeroed(vals, is_viral, time_values, earliest)
        log_vals = np.log2(vals + offset)

        r, p_r = stats.pearsonr(log_vals, dpt)
        rho, p_rho = stats.spearmanr(log_vals, dpt)
        rows.append({
            "gene": gene,
            "class": "viral" if is_viral else "host",
            "pearson_r": r,
            "pearson_p": p_r,
            "spearman_rho": rho,
            "spearman_p": p_rho,
            "mean_at_earliest_time": float(np.mean(vals[time_values == earliest])),
            "max_expression": float(np.max(vals)),
        })

    if not rows:
        return pd.DataFrame(), {}
    corr = pd.DataFrame(rows)
    corr["fdr_q"] = bh_fdr(corr["pearson_p"].values)
    corr["spearman_fdr_q"] = bh_fdr(corr["spearman_p"].values)

    per_gene = {row["gene"]: row for _, row in corr.iterrows()}
    for c in ("pearson_r", "spearman_rho", "mean_at_earliest_time",
              "max_expression"):
        corr[c] = corr[c].astype(float).round(4)
    for c in ("pearson_p", "spearman_p", "fdr_q", "spearman_fdr_q"):
        corr[c] = corr[c].map(lambda v: "%.4e" % v)
    return corr, per_gene


def validate_against_time(adata):
    """Compare pseudotime with the recorded sampling time."""
    t = adata.obs["time"].values.astype(float)
    dpt = adata.obs["dpt_pseudotime"].values.astype(float)

    r, p_r = stats.pearsonr(t, dpt)
    rho, p_rho = stats.spearmanr(t, dpt)
    q = bh_fdr([p_r, p_rho])

    slope, intercept = np.polyfit(t, dpt, 1)
    val = pd.DataFrame([{
        "test": "Pearson", "statistic": r, "p_value": p_r, "fdr_q": q[0],
    }, {
        "test": "Spearman", "statistic": rho, "p_value": p_rho, "fdr_q": q[1],
    }])
    print("  pseudotime versus recorded time: r = %.4f, rho = %.4f, "
          "Spearman p = %.3e" % (r, rho, p_rho))

    # Per time point, so that the spread of the ordering can be inspected.
    by_time = (pd.DataFrame({"time": t, "dpt": dpt})
               .groupby("time")["dpt"]
               .agg(["count", "mean", "std", "min", "max"])
               .reset_index())
    by_time.columns = ["time", "n", "dpt_mean", "dpt_sd", "dpt_min", "dpt_max"]
    return val, by_time, dict(r=r, p_pearson=p_r, rho=rho, p_spearman=p_rho,
                              fdr_spearman=float(q[1]), slope=float(slope),
                              intercept=float(intercept))


# ----------------------------------------------------------------------- output --
def save_results(adata, fpkm_df, args, corr, per_gene, val, by_time, val_stats,
                 n_neighbors):
    """Write every table and the summary file."""
    out = args.output_dir
    os.makedirs(out, exist_ok=True)
    time_values = adata.obs["time"].values
    earliest = min(time_values)
    viral_re = re.compile(args.viral_pattern)
    offset = args.fpkm_offset

    # ---- per-library table, including UMAP coordinates as figure source data
    info = adata.obs.copy()
    info["UMAP1"] = adata.obsm["X_umap"][:, 0]
    info["UMAP2"] = adata.obsm["X_umap"][:, 1]
    diffmap = adata.obsm["X_diffmap"]
    for j in range(1, min(4, diffmap.shape[1])):
        info["DC%d" % j] = diffmap[:, j]
    info["PC1"] = adata.obsm["X_pca"][:, 0]

    for gene in split_csv(args.genes):
        if gene not in fpkm_df.index:
            continue
        vals = fpkm_df.loc[gene].reindex(adata.obs_names).values.astype(float)
        vals = viral_zeroed(vals, bool(viral_re.search(gene)),
                            time_values, earliest)
        info["fpkm_" + gene] = np.round(vals, 4)
        info["log2_" + gene] = np.round(np.log2(vals + offset), 4)

    info = info.reset_index().rename(columns={"index": "sample"})
    info.to_csv(os.path.join(out, "dpt_data.csv"), index=False)
    print("  written: %s" % os.path.join(out, "dpt_data.csv"))

    # ---- per-cell trajectory, already ordered by sampling time
    traj = info[["sample", "cell_id", "time", "dpt_pseudotime",
                 "UMAP1", "UMAP2"]].copy()
    traj = traj.sort_values(["cell_id", "time"])
    traj["order_within_cell"] = traj.groupby("cell_id").cumcount() + 1
    traj.to_csv(os.path.join(out, "per_cell_trajectory.csv"), index=False)
    print("  written: %s" % os.path.join(out, "per_cell_trajectory.csv"))

    # ---- validation and gene correlations
    val.to_csv(os.path.join(out, "dpt_validation.csv"), index=False)
    print("  written: %s" % os.path.join(out, "dpt_validation.csv"))
    by_time.to_csv(os.path.join(out, "pseudotime_by_timepoint.csv"), index=False)
    if len(corr):
        corr.to_csv(os.path.join(out, "gene_pseudotime_correlation.csv"),
                    index=False)
        print("  written: %s"
              % os.path.join(out, "gene_pseudotime_correlation.csv"))

    # ---- the object itself, so that no step has to be repeated for a figure
    try:
        adata.write_h5ad(os.path.join(out, "anndata.h5ad"))
        print("  written: %s" % os.path.join(out, "anndata.h5ad"))
    except Exception as exc:                      # pragma: no cover
        print("  WARNING: could not write the h5ad object: %s" % exc)

    # ---- human-readable summary
    lines = [
        "=" * 64,
        "Diffusion pseudotime summary",
        "=" * 64,
        "",
        "Reference: Haghverdi et al., Nat Methods 2016;13(10):845-848",
        "",
        "Libraries              : %d" % adata.n_obs,
        "Cells                    : %d" % adata.obs["cell_id"].nunique(),
        "Time points              : %s" % ", ".join(
            "%dh" % t for t in sorted(adata.obs["time"].unique())),
        "Genes after filtering    : %d" % adata.n_vars,
        "",
        "-" * 64,
        "Parameters",
        "-" * 64,
        "  normalisation          : counts per million, target %.3g" % args.target_sum,
        "  transform              : log1p",
        "  gene filter            : min count %d, non-zero in >= %d libraries"
        % (args.filter_min_count, args.filter_min_samples),
        "  scaling                : clipped at %g" % args.scale_max_value,
        "  PCA components         : %d" % adata.obsm["X_pca"].shape[1],
        "  neighbours             : k = %d" % n_neighbors,
        "  diffusion components   : %d" % (adata.obsm["X_diffmap"].shape[1] - 1),
        "  root library           : %s" % adata.obs_names[adata.uns["iroot"]],
        "  root selection         : %s"
        % ("given by --root-cell" if args.root_cell
           else "earliest time point, smallest PC1"),
        "  random state           : %d" % args.random_state,
        "  branching              : 0, a single linear ordering",
        "",
        "-" * 64,
        "Pseudotime versus recorded time",
        "-" * 64,
        "  Pearson r              : %.4f (p = %.3e)"
        % (val_stats["r"], val_stats["p_pearson"]),
        "  Spearman rho           : %.4f (p = %.3e, FDR = %.3e)"
        % (val_stats["rho"], val_stats["p_spearman"], val_stats["fdr_spearman"]),
        "  linear fit             : pseudotime = %.4f x time + %.4f"
        % (val_stats["slope"], val_stats["intercept"]),
    ]

    for _, row in by_time.iterrows():
        lines.append("  time %2dh  n = %2d  mean = %.4f  sd = %.4f"
                     % (row["time"], row["n"], row["dpt_mean"], row["dpt_sd"]))

    if len(corr):
        lines += [
            "",
            "-" * 64,
            "Gene expression versus pseudotime",
            "-" * 64,
            "  Viral genes are set to zero at the pre-infection time point.",
            "  Host genes keep their measured basal expression.",
            "  Values are log2(FPKM + %g). FDR is Benjamini-Hochberg across"
            % offset,
            "  the genes tested here.",
            "",
            "  %-12s %-6s %8s %8s %10s" % ("gene", "class", "r", "rho", "FDR"),
        ]
        for _, row in corr.iterrows():
            lines.append("  %-12s %-6s %8.4f %8.4f %10s"
                         % (row["gene"], row["class"], row["pearson_r"],
                            row["spearman_rho"], row["fdr_q"]))

    path = os.path.join(out, "pseudotime_summary.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  written: %s" % path)


def main():
    args = parse_args()

    if not os.path.isdir(args.matrix_dir):
        sys.exit("ERROR: --matrix-dir not found: %s" % args.matrix_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 64)
    print("Diffusion pseudotime analysis")
    print("=" * 64)

    print("\n[1/4] Loading matrices from %s" % args.matrix_dir)
    count_df, fpkm_df, meta, time_points = load_data(args)

    print("\n[2/4] Computing the pseudotime ordering")
    adata, n_neighbors = run_dpt(count_df, meta, args)

    print("\n[3/4] Validating against the recorded sampling time")
    val, by_time, val_stats = validate_against_time(adata)

    print("\n[4/4] Correlating genes with pseudotime and writing results")
    corr, per_gene = correlate_genes(adata, fpkm_df, args)
    save_results(adata, fpkm_df, args, corr, per_gene, val, by_time,
                 val_stats, n_neighbors)

    print("\nAll results are in %s" % args.output_dir)
    print("Done")


if __name__ == "__main__":
    main()
