#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH -J mpa_s2m_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00
#SBATCH --mem=0
#SBATCH -o mpa_s2m/%x_%j.out
#SBATCH -e mpa_s2m/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p mpa_s2m

# ---- Environment ----------------------------------------------------
module purge
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate metaphlan

# ---- Paths ----------------------------------------------------------
RAW_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KNEADDATA"
BASE_OUT="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/METAPHLAN"
DB_DIR="/home/vikkxtar/lu2025-12-46/util/util_12-47/metaphlan_db_202403"
INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
DB_PKL="$DB_DIR/${INDEX}.pkl"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-6}"

# ---- New run folder -------------------------------------------------
RUN_TAG="$(date +%Y%m%d)_${SLURM_JOB_ID}"
OUT_DIR="$BASE_OUT/run_${RUN_TAG}"
PROF_DIR="$OUT_DIR/profiles"
BT2_DIR="$OUT_DIR/bowtie2"
CONS_DIR="$OUT_DIR/consensus_markers"   # sample2markers output

mkdir -p "$PROF_DIR" "$BT2_DIR" "$CONS_DIR"

echo "[$(date)] Writing outputs to: $OUT_DIR"
echo "[$(date)] DB_DIR=$DB_DIR"
echo "[$(date)] INDEX=$INDEX"
echo "[$(date)] DB_PKL=$DB_PKL"
echo "[$(date)] Workers: $MAX_WORKERS  |  Threads/worker: $THREADS"

# ---- Verify SNIC_TMP ------------------------------------------------
if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2
  exit 1
fi

# ---- Verify DB pkl exists ------------------------------------------
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

export RAW_DIR DB_DIR INDEX DB_PKL THREADS PROF_DIR BT2_DIR CONS_DIR

# ---- Process samples in parallel ------------------------------------
for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"

  srun -Q --exclusive -n 1 -N 1 -c "${THREADS}" bash -c '
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
    OUT_MAP="${BT2_DIR}/${SAMPLE}.bowtie2.bz2"
    OUT_JSON="${CONS_DIR}/${SAMPLE}.json.bz2"

    # Skip only if final outputs exist (we do NOT keep SAMs)
    if [[ -s "$OUT_PROFILE" && -s "$OUT_MAP" && -s "$OUT_JSON" ]]; then
      echo "[$(date)] Worker ${WORKER_ID}: Final outputs exist, skipping ${SAMPLE}"
      exit 0
    fi

    WORK="${SNIC_TMP}/mpa_s2m_${SLURM_JOB_ID}_${WORKER_ID}_${SAMPLE}"
    mkdir -p "$WORK/in" "$WORK/tmp" "$WORK/out"
    export TMPDIR="$WORK/tmp"

    # Local copy to reduce Lustre I/O
    cp "$R1" "$R2" "$WORK/in/"

    echo "[$(date)] Worker ${WORKER_ID}: MetaPhlAn (SAM stays local)"
    metaphlan \
      "$WORK/in/${SAMPLE}_paired_1.fastq.gz,$WORK/in/${SAMPLE}_paired_2.fastq.gz" \
      --input_type fastq \
      --db_dir "$DB_DIR" \
      --index "$INDEX" \
      --offline \
      --nproc "'"$THREADS"'" \
      --tmp_dir "$WORK/tmp" \
      --mapout "$WORK/out/${SAMPLE}.bowtie2.bz2" \
      --samout "$WORK/out/${SAMPLE}.sam.bz2" \
      -o "$WORK/out/${SAMPLE}_profile.txt"

    # Verify MetaPhlAn outputs
    if [[ ! -s "$WORK/out/${SAMPLE}_profile.txt" || ! -s "$WORK/out/${SAMPLE}.bowtie2.bz2" || ! -s "$WORK/out/${SAMPLE}.sam.bz2" ]]; then
      echo "ERROR: MetaPhlAn outputs missing for ${SAMPLE}" >&2
      ls -lah "$WORK/out" >&2 || true
      exit 20
    fi

    # Copy permanent outputs
    cp -f "$WORK/out/${SAMPLE}_profile.txt" "'"$PROF_DIR"'/" || exit 21
    cp -f "$WORK/out/${SAMPLE}.bowtie2.bz2" "'"$BT2_DIR"'/" || exit 21

    echo "[$(date)] Worker ${WORKER_ID}: sample2markers on local SAM"
    sample2markers.py \
      -d "$DB_PKL" \
      -i "$WORK/out/${SAMPLE}.sam.bz2" \
      -o "$WORK/out" \
      -f bz2 \
      -n "'"$THREADS"'" \
      --tmp "$WORK/tmp"

    # sample2markers output is .json.bz2
    if [[ ! -s "$WORK/out/${SAMPLE}.json.bz2" ]]; then
      echo "ERROR: sample2markers did not produce ${SAMPLE}.json.bz2" >&2
      ls -lah "$WORK/out" >&2 || true
      exit 30
    fi

    cp -f "$WORK/out/${SAMPLE}.json.bz2" "'"$CONS_DIR"'/" || exit 31

    # Cleanup (removes big SAM)
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

# ---- Final check -----------------------------------------------------
FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
  if [[ ! -s "$PROF_DIR/${SAMPLE}_profile.txt" || ! -s "$BT2_DIR/${SAMPLE}.bowtie2.bz2" || ! -s "$CONS_DIR/${SAMPLE}.json.bz2" ]]; then
    echo "WARNING: ${SAMPLE} missing outputs" >&2
    ((FAILED++))
  fi
done

echo "[$(date)] Completed farm. Outputs: $OUT_DIR"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples incomplete - check metalogs/worker_* logs" >&2
  exit 1
fi

echo "[$(date)] All samples completed successfully"
