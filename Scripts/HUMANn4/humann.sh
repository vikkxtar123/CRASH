#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH --requeue
#SBATCH -J humann4
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH -t 120:00:00
#SBATCH --mem=0
#SBATCH -o mpa_humann/%x_%j.out
#SBATCH -e mpa_humann/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p mpa_humann metalogs

module purge
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate humann4

RAW_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KNEADDATA"
BASE_OUT="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/HUMANN4"
PROF_DIR="$BASE_OUT/profiles"
THREADS=24        # 48 cores / 2 workers
MAX_WORKERS=2

export TMPDIR="$SNIC_TMP"

RUN_TAG="$(date +%Y%m%d)_${SLURM_JOB_ID}"
HMN_DIR="$BASE_OUT/run_${RUN_TAG}/humann"
mkdir -p "$HMN_DIR" "$PROF_DIR"

echo "[$(date)] Output: $HMN_DIR"
echo "[$(date)] Threads/worker: $THREADS  |  Workers: $MAX_WORKERS"

if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2; exit 1
fi

mapfile -t SAMPLES < <(
  find "$RAW_DIR" -maxdepth 1 -type f -name 'SRR*_paired_1.fastq.gz' -printf '%f\n' \
  | sed -E 's/_paired_1\.fastq\.gz$//' | sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
  echo "ERROR: No samples found in $RAW_DIR" >&2; exit 3
fi

echo "[$(date)] Found ${#SAMPLES[@]} samples."

process_sample() {
  local SAMPLE="$1"
  local WORK="${SNIC_TMP}/hmn_${SLURM_JOB_ID}_${SAMPLE}"
  trap "rm -rf '$WORK'" RETURN

  local OUT_PROFILE="${PROF_DIR}/${SAMPLE}_profile.txt"
  local OUT_HMN="${HMN_DIR}/${SAMPLE}/${SAMPLE}_2_genefamilies.tsv"

  if [[ -s "$OUT_HMN" ]]; then
    echo "[$(date)] ${SAMPLE}: already complete, skipping"; return 0
  fi

  if [[ ! -s "$OUT_PROFILE" ]]; then
    echo "ERROR: ${SAMPLE}: profile missing - run phase1 first" >&2; return 1
  fi

  mkdir -p "$WORK/out"

  echo "[$(date)] ${SAMPLE}: concatenating reads to scratch"
  zcat "${RAW_DIR}/${SAMPLE}_paired_1.fastq.gz" \
       "${RAW_DIR}/${SAMPLE}_paired_2.fastq.gz" \
       > "$WORK/${SAMPLE}.fastq"

  echo "[$(date)] ${SAMPLE}: running HUMAnN"
  humann \
    --input "$WORK/${SAMPLE}.fastq" \
    --taxonomic-profile "$OUT_PROFILE" \
    --output "$WORK/out" \
    --output-basename "$SAMPLE" \
    --threads "$THREADS" \
    --prescreen-threshold 1.0 \
    --remove-temp-output \
    2>&1 | tee "metalogs/humann_${SLURM_JOB_ID}_${SAMPLE}.log"

  if [[ ! -s "$WORK/out/${SAMPLE}_2_genefamilies.tsv" ]]; then
    echo "ERROR: ${SAMPLE}: HUMAnN outputs missing" >&2
    ls -lah "$WORK/out" >&2 || true
    return 1
  fi

  mkdir -p "${HMN_DIR}/${SAMPLE}"
  cp "$WORK/out/${SAMPLE}_"*.tsv "${HMN_DIR}/${SAMPLE}/"

  for f in "$WORK/out/${SAMPLE}_"*.tsv; do
    echo "${f#${SNIC_TMP}/}" >> "$SNIC_TMP/slurm_save_files"
  done

  echo "[$(date)] ${SAMPLE}: complete"
  return 0
}

# ---- Parallel main loop ---------------------------------------------
FAILED=0
PIDS=()

for SAMPLE in "${SAMPLES[@]}"; do
  process_sample "$SAMPLE" &
  PIDS+=($!)

  # Throttle to MAX_WORKERS concurrent jobs
  while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do
    sleep 10
  done
done

# Wait for all and collect exit codes
for PID in "${PIDS[@]}"; do
  wait "$PID" || ((FAILED++))
done

echo "[$(date)] Done. Outputs: $HMN_DIR"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples failed - check metalogs/humann_${SLURM_JOB_ID}_*.log" >&2
  exit 1
fi
echo "[$(date)] All samples completed successfully"