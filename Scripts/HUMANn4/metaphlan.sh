#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH --requeue
#SBATCH -J mpa_oct22
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH -t 24:00:00
#SBATCH -o mpa_oct22/%x_%j.out
#SBATCH -e mpa_oct22/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p mpa_oct22 metalogs

# ---- Environment ----------------------------------------------------
module purge
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate humann4   # contains metaphlan

# ---- Paths ----------------------------------------------------------
RAW_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KNEADDATA"
BASE_OUT="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/METAPHLAN_OCT22"
DB_DIR="/lunarc/nobackup/projects/lu2025-12-46/util/util_12-47/metaphlan_db_oct22"
INDEX="mpa_vOct22_CHOCOPhlAnSGB_202403"
DB_PKL="${DB_DIR}/${INDEX}.pkl"

THREADS="${SLURM_CPUS_PER_TASK:-12}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-4}"

# ---- New run folder -------------------------------------------------
RUN_TAG="$(date +%Y%m%d)_${SLURM_JOB_ID}"
OUT_DIR="${BASE_OUT}/run_${RUN_TAG}"
PROF_DIR="${OUT_DIR}/profiles"
mkdir -p "$OUT_DIR" "$PROF_DIR"

echo "[$(date)] Writing outputs to: $OUT_DIR"
echo "[$(date)] DB_DIR=$DB_DIR  |  INDEX=$INDEX"
echo "[$(date)] Workers: $MAX_WORKERS  |  Threads/worker: $THREADS"

# ---- Verify SNIC_TMP ------------------------------------------------
if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2
  exit 1
fi

# ---- Verify DB pkl exists -------------------------------------------
if [[ ! -s "$DB_PKL" ]]; then
  echo "ERROR: DB pkl not found: $DB_PKL" >&2
  ls -lh "$DB_DIR" | head >&2 || true
  exit 2
fi

# ---- Build sample list ----------------------------------------------
mapfile -t SAMPLES < <(
  find "$RAW_DIR" -maxdepth 1 -type f -name 'SRR*_paired_1.fastq.gz' -printf '%f\n' \
  | sed -E 's/_paired_1\.fastq\.gz$//' \
  | sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
  echo "ERROR: No samples found in $RAW_DIR (expected SRR*_paired_1.fastq.gz)" >&2
  exit 3
fi

echo "[$(date)] Found ${#SAMPLES[@]} samples."

export RAW_DIR DB_DIR INDEX THREADS PROF_DIR

# ---- Process samples in parallel ------------------------------------
for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"

  srun --exclusive -n 1 -N 1 -c "${THREADS}" bash -c '
    set -euo pipefail

    SAMPLE="'"${SAMPLE}"'"
    WORKER_ID="'"${i}"'"

    echo "[$(date)] Worker ${WORKER_ID}: Starting ${SAMPLE} on $(hostname)"

    R1="${RAW_DIR}/${SAMPLE}_paired_1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}_paired_2.fastq.gz"

    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "ERROR: Missing paired FASTQs for ${SAMPLE}" >&2
      exit 10
    fi

    OUT_PROFILE="${PROF_DIR}/${SAMPLE}_profile.txt"
    if [[ -s "$OUT_PROFILE" ]]; then
      echo "[$(date)] Worker ${WORKER_ID}: Profile exists, skipping ${SAMPLE}"
      exit 0
    fi

    WORK="${SNIC_TMP}/mpa_${SLURM_JOB_ID}_${WORKER_ID}_${SAMPLE}"
    mkdir -p "$WORK/in" "$WORK/tmp" "$WORK/out"
    export TMPDIR="$WORK/tmp"

    # Local copy to reduce Lustre I/O
    cp "$R1" "$R2" "$WORK/in/"

    echo "[$(date)] Worker ${WORKER_ID}: Running MetaPhlAn"
    metaphlan \
      "$WORK/in/${SAMPLE}_paired_1.fastq.gz,$WORK/in/${SAMPLE}_paired_2.fastq.gz" \
      --input_type fastq \
      --bowtie2db "$DB_DIR" \
      -x "$INDEX" \
      --offline \
      --no_map \
      --tmp_dir "$WORK/tmp" \
      --nproc "'"$THREADS"'" \
      -t rel_ab_w_read_stats \
      -o "$WORK/out/${SAMPLE}_profile.txt"

    if [[ ! -s "$WORK/out/${SAMPLE}_profile.txt" ]]; then
      echo "ERROR: MetaPhlAn profile missing for ${SAMPLE}" >&2
      ls -lah "$WORK/out" >&2 || true
      exit 20
    fi

    cp -f "$WORK/out/${SAMPLE}_profile.txt" "$OUT_PROFILE"

    rm -rf "$WORK"
    echo "[$(date)] Worker ${WORKER_ID}: Completed ${SAMPLE}"
  ' &> "metalogs/worker_${SLURM_JOB_ID}_${i}.log" &

  # Throttle concurrent workers
  while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do
    sleep 2
  done
done

echo "[$(date)] All workers launched, waiting..."
wait

# ---- Final check ----------------------------------------------------
FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
  if [[ ! -s "$PROF_DIR/${SAMPLE}_profile.txt" ]]; then
    echo "WARNING: ${SAMPLE} missing profile" >&2
    ((FAILED++))
  fi
done

echo "[$(date)] Completed farm. Outputs: $OUT_DIR"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples incomplete - check metalogs/worker_* logs" >&2
  exit 1
fi

echo "[$(date)] All samples completed successfully"