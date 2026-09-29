#!/usr/bin/env bash
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Step 2 of 2 -- Consolidated per-library QC table
# ------------------------------------------------------------------------------
# What it does
#   Merges the first-pass statistics written by 01_align_and_quantify.sh with
#   four further per-library metrics and produces one table with one row per
#   single-cell library. This is the table behind the library-quality figures
#   of the manuscript.
#
#   Columns added on top of summary_stats.tsv
#     MappedPct    uniquely mapped reads, percent, from the STAR log
#     nCount_1M    total assigned fragments at the common down-sampled depth
#     MT_1M_Pct    percent of uniquely mapped reads on the mitochondrial
#                  contig at the common depth
#     nCount_All   total assigned fragments in the full, un-normalised library
#     nGene_All    number of genes with a non-zero count in the full library
#
# Input
#   --work-dir must be the same directory passed to 01_align_and_quantify.sh
#
# Output
#   <work-dir>/02_results/full_stats_complete.tsv
#
# Requirements
#   samtools  awk  bc
#
# Usage
#   bash 02_library_statistics.sh --work-dir /path/to/project
#
# Optional flags
#   --mt-contig NAME   mitochondrial contig name in the reference
#                      (default chrM, use MT for Ensembl-style references)
#   --mapq N           minimum MAPQ for a read to count as uniquely mapped
#                      (default 255, matching STAR unique mappers)
# ==============================================================================

set -euo pipefail

WORK_DIR=""
MT_CHR="chrM"
MAPQ=255

usage() { sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --work-dir)  WORK_DIR="$2"; shift 2 ;;
    --mt-contig) MT_CHR="$2";   shift 2 ;;
    --mapq)      MAPQ="$2";     shift 2 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "ERROR: unknown option '$1'" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$WORK_DIR" ]] || { echo "ERROR: --work-dir is required" >&2; usage; exit 1; }

STAR_DIR="$WORK_DIR/05_star"
QUANT_DIR="$WORK_DIR/06_quantify"
DS_BAM_DIR="$WORK_DIR/07_downsampled_bam"
DS_QUANT_DIR="$WORK_DIR/08_downsampled_quant"
RESULTS_DIR="$WORK_DIR/02_results"
LOGS_DIR="$WORK_DIR/03_logs"

SUMMARY_IN="$RESULTS_DIR/summary_stats.tsv"
FINAL_OUT="$RESULTS_DIR/full_stats_complete.tsv"
BAM_LIST="$QUANT_DIR/bam_list.txt"
ORIG_QUANT="$QUANT_DIR/quantification_results.txt"
LOG_FILE="$LOGS_DIR/full_stats.log"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR"

for f in "$SUMMARY_IN" "$BAM_LIST" "$ORIG_QUANT"; do
  [[ -f "$f" ]] || { echo "ERROR: expected file not found: $f" >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Locate the column of each sample inside the merged featureCounts table.
# The original version of this script derived the column from the row index,
# which silently breaks as soon as the two input files are sorted differently.
# Reading the header instead makes the join explicit and order-independent.
# ---------------------------------------------------------------------------
declare -A COL_OF
header=$(head -1 "$ORIG_QUANT")
i=0
while IFS=$'\t' read -r field; do
  i=$((i + 1))
  if [[ "$field" == */Aligned.sortedByCoord.out.bam ]]; then
    s=$(basename "$(dirname "$field")")
    COL_OF["$s"]=$i
  fi
done <<< "$header"

if [[ ${#COL_OF[@]} -eq 0 ]]; then
  echo "ERROR: no sample columns recognised in the header of $ORIG_QUANT" >&2
  echo "       Expected each column header to be a path ending in" >&2
  echo "       /Aligned.sortedByCoord.out.bam" >&2
  exit 1
fi

echo "Matched ${#COL_OF[@]} sample columns in the merged count table" >> "$LOG_FILE"

# Sample order is taken from bam_list.txt so that the output rows follow the
# alignment order rather than whatever order the first-pass table happens to use.
mapfile -t SAMPLES < <(sed 's|.*/\([^/]*\)/Aligned.sortedByCoord.out.bam|\1|' "$BAM_LIST")

{
  head -1 "$SUMMARY_IN"
} | awk -F'\t' '{print $0"\tMappedPct\tnCount_1M\tMT_1M_Pct\tnCount_All\tnGene_All"}' > "$FINAL_OUT"

# summary_stats.tsv is keyed on SampleID, so build a lookup for fast access
declare -A SUMMARY_ROW
while IFS=$'\t' read -r sid rest; do
  [[ "$sid" == "SampleID" ]] && continue
  SUMMARY_ROW["$sid"]="$rest"
done < "$SUMMARY_IN"

for s in "${SAMPLES[@]}"; do
  echo "Processing $s" >> "$LOG_FILE"

  first_pass="${SUMMARY_ROW[$s]:-}"
  if [[ -z "$first_pass" ]]; then
    echo "WARNING: $s absent from summary_stats.tsv, row filled with NA" >> "$LOG_FILE"
    first_pass=$'NA\tNA\tNA\tNA\tNA\tNA'
  fi
  down_frags=$(cut -f2 <<< "$first_pass")

  # 1. Uniquely mapped reads, percent
  star_log="$STAR_DIR/$s/Log.final.out"
  if [[ -f "$star_log" ]]; then
    mapped_pct=$(awk -F'\t' '/^Uniquely mapped reads %/ {gsub(/%/,"",$2); print $2}' "$star_log")
    if [[ -z "$mapped_pct" ]]; then
      mapped_pct=$(grep "Uniquely mapped reads %" "$star_log" | awk -F '|' '{print $2}' | tr -d ' %')
    fi
    [[ -n "$mapped_pct" ]] || mapped_pct="NA"
  else
    mapped_pct="NA"
  fi

  # 2. Assigned fragments at the common depth
  quant_1m="$DS_QUANT_DIR/${s}.quant.txt"
  if [[ -f "$quant_1m" ]]; then
    ncount_1m=$(awk 'NR>2 {sum += $7} END {printf "%.0f", sum}' "$quant_1m")
  else
    ncount_1m="NA"
  fi

  # 3. Mitochondrial fraction at the common depth
  ds_bam="$DS_BAM_DIR/${s}.downsampled.bam"
  if [[ ! -f "$ds_bam" ]]; then
    # Libraries that were already below the common depth were not re-written by
    # step 1, so the original BAM is the depth-normalised BAM for those cells.
    orig_bam="$STAR_DIR/$s/Aligned.sortedByCoord.out.bam"
    if [[ -f "$orig_bam" && "$down_frags" != "NA" && "$down_frags" -le 1000000 ]]; then
      ds_bam="$orig_bam"
      echo "Note: original BAM reused for $s (fragments <= common depth)" >> "$LOG_FILE"
    fi
  fi
  if [[ -f "$ds_bam" ]]; then
    [[ -f "${ds_bam}.bai" ]] || samtools index "$ds_bam" 2>/dev/null || true
    total_unique=$(samtools view -c -F 0x4 -q "$MAPQ" "$ds_bam" 2>/dev/null || echo 0)
    mt_unique=$(samtools view -c -F 0x4 -q "$MAPQ" "$ds_bam" "$MT_CHR" 2>/dev/null || echo 0)
    if [[ -n "$total_unique" && "$total_unique" -gt 0 ]]; then
      mt_pct=$(echo "scale=2; ${mt_unique} * 100 / ${total_unique}" | bc)
    else
      mt_pct=0
    fi
  else
    mt_pct="NA"
  fi

  # 4. Full-library totals from the merged count table
  col="${COL_OF[$s]:-}"
  if [[ -n "$col" ]]; then
    read -r ncount_all ngene_all <<< "$(awk -v c="$col" '!/^#/ && NR>1 { if ($c > 0) {sum += $c; n++} } END {printf "%d %d", sum, n}' "$ORIG_QUANT")"
  else
    ncount_all="NA"
    ngene_all="NA"
    echo "WARNING: $s not found in the merged count table header" >> "$LOG_FILE"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$s" "$first_pass" "$mapped_pct" "$ncount_1m" "$mt_pct" "$ncount_all" "$ngene_all" \
    >> "$FINAL_OUT"
done

echo "Completed. Output: $FINAL_OUT" >> "$LOG_FILE"
echo "Done. Consolidated QC table: $FINAL_OUT"
