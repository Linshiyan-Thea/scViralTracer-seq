#!/usr/bin/env bash
# ==============================================================================
# scViralTracer-seq analysis pipeline
# Step 1 of 2 -- Read QC, trimming, alignment, quantification, depth normalising
# ------------------------------------------------------------------------------
# What it does
#   1. FastQC + MultiQC on the raw reads
#   2. Adapter and quality trimming with Trim Galore
#   3. STAR alignment against a COMBINED host + viral reference
#   4. featureCounts gene-level quantification (all samples in one table)
#   5. Down-sampling of every library to a common fragment depth, followed by a
#      second featureCounts pass, so that library size does not confound the
#      per-cell comparisons reported in the manuscript
#
# Input
#   A directory of paired-end FASTQ files named
#       <SAMPLE>_1.fq.gz   <SAMPLE>_2.fq.gz
#   The suffix .fastq.gz is accepted as well. <SAMPLE> must follow the naming
#   convention described in the README, because every downstream script parses
#   the cell identifier and the sampling time point out of it.
#
# Output (all under --work-dir)
#   02_results/            MultiQC report and summary_stats.tsv
#   04_clean_data/         trimmed reads
#   05_star/<SAMPLE>/      Aligned.sortedByCoord.out.bam, Log.final.out
#   06_quantify/           quantification_results.txt, bam_list.txt
#   07_downsampled_bam/    depth-normalised BAM files
#   08_downsampled_quant/  per-sample count tables at the normalised depth
#   03_logs/               one log file per step
#
# Requirements
#   fastqc  multiqc  trim_galore  STAR  samtools  subread (featureCounts)  bc
#
# Usage
#   bash 01_align_and_quantify.sh \
#        --fastq-dir  /path/to/raw_fastq \
#        --work-dir   /path/to/project \
#        --star-index /path/to/star_index \
#        --gtf        /path/to/combined_host_virus.gtf \
#        --threads 12
#
# Optional flags
#   --tmp-dir DIR            scratch space for featureCounts (default WORK_DIR/tmp)
#   --target-fragments N     common depth for down-sampling (default 1000000)
#   --seed N                 down-sampling seed, fixed for reproducibility
#   --skip-fastqc            skip step 1
#   --skip-trim              use the raw reads directly for alignment
#
# Optional scheduler header. Uncomment and adapt if you submit to a cluster.
#   #PBS -l nodes=1:ppn=12
#   #PBS -l mem=100gb
#   #PBS -l walltime=84:00:00
# ==============================================================================

set -euo pipefail

# --------------------------------------------------------------- defaults ----
FASTQ_DIR=""
WORK_DIR=""
STAR_INDEX=""
GTF=""
TMP_DIR=""
THREADS=12
TARGET_FRAGMENTS=1000000
SEED=123
SKIP_FASTQC=0
SKIP_TRIM=0

usage() {
  sed -n '2,60p' "$0" | sed 's/^# \{0,1\}//'
}

# ------------------------------------------------------------ parse arguments -
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq-dir)        FASTQ_DIR="$2";        shift 2 ;;
    --work-dir)         WORK_DIR="$2";         shift 2 ;;
    --star-index)       STAR_INDEX="$2";       shift 2 ;;
    --gtf)              GTF="$2";              shift 2 ;;
    --tmp-dir)          TMP_DIR="$2";          shift 2 ;;
    --threads)          THREADS="$2";          shift 2 ;;
    --target-fragments) TARGET_FRAGMENTS="$2"; shift 2 ;;
    --seed)             SEED="$2";             shift 2 ;;
    --skip-fastqc)      SKIP_FASTQC=1;         shift ;;
    --skip-trim)        SKIP_TRIM=1;           shift ;;
    -h|--help)          usage; exit 0 ;;
    *) echo "ERROR: unknown option '$1'" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$FASTQ_DIR"  ]] || { echo "ERROR: --fastq-dir is required"  >&2; exit 1; }
[[ -n "$WORK_DIR"   ]] || { echo "ERROR: --work-dir is required"   >&2; exit 1; }
[[ -n "$STAR_INDEX" ]] || { echo "ERROR: --star-index is required" >&2; exit 1; }
[[ -n "$GTF"        ]] || { echo "ERROR: --gtf is required"        >&2; exit 1; }
[[ -d "$FASTQ_DIR"  ]] || { echo "ERROR: --fastq-dir not found: $FASTQ_DIR" >&2; exit 1; }
[[ -d "$STAR_INDEX" ]] || { echo "ERROR: --star-index not found: $STAR_INDEX" >&2; exit 1; }
[[ -f "$GTF"        ]] || { echo "ERROR: --gtf not found: $GTF" >&2; exit 1; }
[[ -n "$TMP_DIR"    ]] || TMP_DIR="$WORK_DIR/tmp"

# ------------------------------------------------------------------ layout ----
RESULTS_DIR="$WORK_DIR/02_results"
LOGS_DIR="$WORK_DIR/03_logs"
CLEAN_DIR="$WORK_DIR/04_clean_data"
STAR_DIR="$WORK_DIR/05_star"
QUANT_DIR="$WORK_DIR/06_quantify"
DS_BAM_DIR="$WORK_DIR/07_downsampled_bam"
DS_QUANT_DIR="$WORK_DIR/08_downsampled_quant"
SUMMARY_TSV="$RESULTS_DIR/summary_stats.tsv"

mkdir -p "$RESULTS_DIR" "$LOGS_DIR" "$CLEAN_DIR" "$STAR_DIR" "$QUANT_DIR" \
         "$DS_BAM_DIR" "$DS_QUANT_DIR" "$TMP_DIR"

# --------------------------------------------------------- discover samples ----
# Both .fq.gz and .fastq.gz are accepted. Names are sorted so that the sample
# order is identical in every table this pipeline writes.
shopt -s nullglob
SAMPLES=()
for f in "$FASTQ_DIR"/*_1.fq.gz "$FASTQ_DIR"/*_1.fastq.gz; do
  b=$(basename "$f")
  b=${b%_1.fq.gz}
  b=${b%_1.fastq.gz}
  SAMPLES+=("$b")
done
shopt -u nullglob

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no *_1.fq.gz or *_1.fastq.gz files in $FASTQ_DIR" >&2
  exit 1
fi

IFS=$'\n' SAMPLES=($(printf '%s\n' "${SAMPLES[@]}" | sort -u)); unset IFS

echo "Found ${#SAMPLES[@]} samples:"
printf '  %s\n' "${SAMPLES[@]}"

# Resolve the R1 / R2 pair of one sample, whatever the FASTQ suffix is.
read1_of() {
  local s="$1"
  for cand in "$FASTQ_DIR/${s}_1.fq.gz" "$FASTQ_DIR/${s}_1.fastq.gz"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}
read2_of() {
  local s="$1"
  for cand in "$FASTQ_DIR/${s}_2.fq.gz" "$FASTQ_DIR/${s}_2.fastq.gz"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}

# ============================================================ step 1: FastQC ==
if [[ $SKIP_FASTQC -eq 0 ]]; then
  echo "[1/5] FastQC and MultiQC"
  fastqc --threads "$THREADS" --outdir "$RESULTS_DIR" \
         $(for s in "${SAMPLES[@]}"; do read1_of "$s"; read2_of "$s"; done) \
         >> "$LOGS_DIR/qc.log" 2>&1
  multiqc "$RESULTS_DIR" -o "$RESULTS_DIR" -f >> "$LOGS_DIR/qc.log" 2>&1
else
  echo "[1/5] FastQC skipped"
fi

# ======================================================= step 2: Trim Galore ==
if [[ $SKIP_TRIM -eq 0 ]]; then
  echo "[2/5] Trim Galore (quality 25, stringency 1, length 20, paired)"
  for s in "${SAMPLES[@]}"; do
    trim_galore --quality 25 --stringency 1 --length 20 --paired \
                --cores 2 --output_dir "$CLEAN_DIR" \
                "$(read1_of "$s")" "$(read2_of "$s")" \
                >> "$LOGS_DIR/trim.log" 2>&1
  done
else
  echo "[2/5] Trimming skipped, raw reads will be aligned"
  ln -sfn "$FASTQ_DIR" "$CLEAN_DIR"
fi

# Trim Galore writes <SAMPLE>_1_val_1.fq.gz. Resolve whichever suffix appeared.
trimmed1_of() {
  local s="$1"
  for cand in "$CLEAN_DIR/${s}_1_val_1.fq.gz" "$CLEAN_DIR/${s}_1_val_1.fastq.gz"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}
trimmed2_of() {
  local s="$1"
  for cand in "$CLEAN_DIR/${s}_2_val_2.fq.gz" "$CLEAN_DIR/${s}_2_val_2.fastq.gz"; do
    [[ -f "$cand" ]] && { echo "$cand"; return 0; }
  done
  return 1
}

# ========================================================= step 3: STAR align ==
echo "[3/5] STAR alignment"
for s in "${SAMPLES[@]}"; do
  out="$STAR_DIR/$s"
  mkdir -p "$out"
  STAR --runThreadN "$THREADS" \
       --genomeDir "$STAR_INDEX" \
       --readFilesIn "$(trimmed1_of "$s")" "$(trimmed2_of "$s")" \
       --readFilesCommand zcat \
       --outFileNamePrefix "$out/" \
       --outSAMtype BAM SortedByCoordinate \
       --outSAMunmapped Within \
       --outSAMattributes Standard \
       >> "$LOGS_DIR/star.log" 2>&1
done

# ================================================ step 4: featureCounts merge ==
echo "[4/5] featureCounts on all samples"
BAM_LIST="$QUANT_DIR/bam_list.txt"
: > "$BAM_LIST"
for s in "${SAMPLES[@]}"; do
  bam="$STAR_DIR/$s/Aligned.sortedByCoord.out.bam"
  [[ -f "$bam" ]] || { echo "ERROR: missing BAM for $s" >&2; exit 1; }
  echo "$bam" >> "$BAM_LIST"
done

featureCounts -a "$GTF" -T "$THREADS" -p --countReadPairs \
              -g gene_id -t exon -M -O \
              --tmpDir "$TMP_DIR" \
              -o "$QUANT_DIR/quantification_results.txt" \
              $(cat "$BAM_LIST") >> "$LOGS_DIR/quantify.log" 2>&1

# ================================ step 5: down-sample and re-quantify per cell ==
echo "[5/5] Down-sampling to ${TARGET_FRAGMENTS} fragments and re-quantifying"
printf 'SampleID\tRawTotalFragments\tDownsampledFragments\tNumGenes\tExonicRate\tUniqueMappingRate\tMultiMappingRate\n' \
  > "$SUMMARY_TSV"

for s in "${SAMPLES[@]}"; do
  bam="$STAR_DIR/$s/Aligned.sortedByCoord.out.bam"
  star_log="$STAR_DIR/$s/Log.final.out"
  echo "  $s" >> "$LOGS_DIR/downsample.log"

  if [[ ! -f "$bam" || ! -f "$star_log" ]]; then
    echo "WARNING: BAM or STAR log missing for $s, skipped" >> "$LOGS_DIR/downsample.log"
    continue
  fi

  # Mapping rates straight from the STAR log
  raw_reads=$(grep "Number of input reads" "$star_log" | awk -F '|' '{print $2}' | tr -d ' ')
  unique_rate=$(grep "Uniquely mapped reads %" "$star_log" | awk -F '|' '{print $2}' | tr -d ' %')
  multi_rate=$(grep "% of reads mapped to multiple loci" "$star_log" | awk -F '|' '{print $2}' | tr -d ' %')

  # Proper pairs only, so that the depth unit is a fragment rather than a read
  total_frags=$(samtools view -c -f 0x2 "$bam" 2>/dev/null || echo 0)
  if [[ -z "$total_frags" || "$total_frags" -eq 0 ]]; then
    echo "WARNING: no proper pairs in $bam, skipped" >> "$LOGS_DIR/downsample.log"
    continue
  fi

  if [[ "$total_frags" -le "$TARGET_FRAGMENTS" ]]; then
    # Libraries already below the common depth are kept as they are. Forcing
    # them up is impossible and forcing them down would only throw away data.
    ds_bam="$bam"
    ds_frags="$total_frags"
    echo "    ${total_frags} fragments <= target, original BAM reused" >> "$LOGS_DIR/downsample.log"
  else
    fraction=$(awk -v t="$TARGET_FRAGMENTS" -v f="$total_frags" 'BEGIN {printf "%.6f", t/f}')
    ds_bam="$DS_BAM_DIR/${s}.downsampled.bam"
    samtools view -h -b -s "${SEED}.${fraction#0.}" -f 0x2 "$bam" > "$ds_bam" \
      2>> "$LOGS_DIR/downsample.log"
    [[ -s "$ds_bam" ]] || { echo "ERROR: empty down-sampled BAM for $s" >&2; continue; }
    samtools index "$ds_bam"
    ds_frags=$(samtools view -c -f 0x2 "$ds_bam")
  fi

  quant_file="$DS_QUANT_DIR/${s}.quant.txt"
  featureCounts -a "$GTF" -T "$THREADS" -p --countReadPairs \
                -g gene_id -t exon -M -O \
                -o "$quant_file" "$ds_bam" \
                > "$DS_QUANT_DIR/${s}.log" 2>&1
  [[ -s "$quant_file" ]] || { echo "ERROR: featureCounts failed for $s" >&2; continue; }

  # Exonic rate from the featureCounts log
  assigned=$(grep "Successfully assigned alignments" "$DS_QUANT_DIR/${s}.log" | head -1 | grep -oE '[0-9,]+' | head -1 | tr -d ',')
  total_aln=$(grep "Total alignments" "$DS_QUANT_DIR/${s}.log" | head -1 | grep -oE '[0-9,]+' | head -1 | tr -d ',')
  if [[ -n "${assigned:-}" && -n "${total_aln:-}" && "$total_aln" -gt 0 ]]; then
    exonic_rate=$(echo "scale=4; ${assigned}/${total_aln} * 100" | bc)
  else
    exonic_rate=0
  fi

  # Number of genes with a non-zero count
  num_genes=$(awk 'NR>2 && $7>0' "$quant_file" | wc -l)

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$s" "$raw_reads" "$ds_frags" "$num_genes" "$exonic_rate" "$unique_rate" "$multi_rate" \
    >> "$SUMMARY_TSV"
done

echo "Done. First-pass statistics: $SUMMARY_TSV"
echo "Next: bash 02_library_statistics.sh --work-dir $WORK_DIR"
