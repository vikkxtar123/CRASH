#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH --requeue
#SBATCH -J humann4
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH -t 120:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs metalogs

DATASET="${1:?Usage: sbatch humann4_run.sh <dataset>}"

source /path/to/config.sh
resolve_dataset "$DATASET"
load_miniforge
conda activate humann4

KNEADDATA_DIR="${DATASET_DIR}/KNEADDATA"
PROF_DIR="${DATASET_DIR}/METAPHLAN_OCT22/profiles"
HMN_DIR="${DATASET_DIR}/HUMANN4/humann_final"

THREADS=24
MAX_WORKERS=2

mkdir -p "$HMN_DIR"
export TMPDIR="$SNIC_TMP"

if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2; exit 1
fi

mapfile -t SAMPLES < <(
  find "$KNEADDATA_DIR" -maxdepth 1 -type f -name '*_paired_1.fastq.gz' -printf '%f\n' \
  | sed -E 's/_paired_1\.fastq\.gz$//' | sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
  echo "ERROR: No samples found in $KNEADDATA_DIR" >&2; exit 3
fi

echo "[$(date)] Dataset: $DATASET | ${#SAMPLES[@]} samples"
echo "[$(date)] Threads/worker: $THREADS | Workers: $MAX_WORKERS"

process_sample() {
  local SAMPLE="$1"
  local WORK="${SNIC_TMP}/hmn_${SLURM_JOB_ID}_${SAMPLE}"
  trap "rm -rf '$WORK'" RETURN

  local OUT_HMN="${HMN_DIR}/${SAMPLE}/${SAMPLE}_2_genefamilies.tsv"
  if [[ -s "$OUT_HMN" ]]; then
    echo "[$(date)] ${SAMPLE}: already complete, skipping"; return 0
  fi

  local OUT_PROFILE="${PROF_DIR}/${SAMPLE}_profile.txt"
  if [[ ! -s "$OUT_PROFILE" ]]; then
    echo "ERROR: ${SAMPLE}: profile missing   run mpa_oct22.sh first" >&2; return 1
  fi

  mkdir -p "$WORK/out"

  echo "[$(date)] ${SAMPLE}: concatenating reads to scratch"
  zcat "${KNEADDATA_DIR}/${SAMPLE}_paired_1.fastq.gz" \
       "${KNEADDATA_DIR}/${SAMPLE}_paired_2.fastq.gz" \
       > "$WORK/${SAMPLE}.fastq"

  echo "[$(date)] ${SAMPLE}: running HUMAnN4"
  humann \
    --input              "$WORK/${SAMPLE}.fastq" \
    --taxonomic-profile  "$OUT_PROFILE" \
    --output             "$WORK/out" \
    --output-basename    "$SAMPLE" \
    --threads            "$THREADS" \
    --prescreen-threshold 1.0 \
    --remove-temp-output \
    2>&1 | tee "metalogs/humann4_${SLURM_JOB_ID}_${SAMPLE}.log"

  if [[ ! -s "$WORK/out/${SAMPLE}_2_genefamilies.tsv" ]]; then
    echo "ERROR: ${SAMPLE}: HUMAnN4 outputs missing" >&2
    ls -lah "$WORK/out" >&2 || true
    return 1
  fi

  mkdir -p "${HMN_DIR}/${SAMPLE}"
  cp "$WORK/out/${SAMPLE}_"*.tsv "${HMN_DIR}/${SAMPLE}/"

  for f in "$WORK/out/${SAMPLE}_"*.tsv; do
    echo "${f#${SNIC_TMP}/}" >> "$SNIC_TMP/slurm_save_files"
  done

  echo "[$(date)] ${SAMPLE}: complete"
}

FAILED=0
PIDS=()
for SAMPLE in "${SAMPLES[@]}"; do
  process_sample "$SAMPLE" &
  PIDS+=($!)
  while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do
    sleep 10
  done
done

for PID in "${PIDS[@]}"; do
  wait "$PID" || ((FAILED++))
done

echo "[$(date)] Done. Outputs: $HMN_DIR"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples failed   check metalogs/humann4_${SLURM_JOB_ID}_*.log" >&2
  exit 1
fi
echo "[$(date)] All samples completed successfully"