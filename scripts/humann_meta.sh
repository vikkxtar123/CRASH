#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH --requeue
#SBATCH -J mpa_oct22
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=12
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs metalogs

DATASET="${1:?Usage: sbatch mpa_oct22.sh <dataset>}"

source /path/to/config.sh
resolve_dataset "$DATASET"

KNEADDATA_DIR="${DATASET_DIR}/KNEADDATA"
PROF_DIR="${DATASET_DIR}/METAPHLAN_OCT22/profiles"

MPA_DB_OCT22="/path/to/metaphlan_db_oct22"
INDEX="mpa_vOct22_CHOCOPhlAnSGB_202403"
DB_PKL="${MPA_DB_OCT22}/${INDEX}.pkl"

THREADS="${SLURM_CPUS_PER_TASK:-12}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-4}"

mkdir -p "$PROF_DIR"

if [[ ! -s "$DB_PKL" ]]; then
  echo "ERROR: DB pkl not found: $DB_PKL" >&2
  ls -lh "$MPA_DB_OCT22" | head >&2 || true
  exit 2
fi

mapfile -t SAMPLES < <(
  find "$KNEADDATA_DIR" -maxdepth 1 -type f -name '*_paired_1.fastq.gz' -printf '%f\n' \
  | sed -E 's/_paired_1\.fastq\.gz$//' \
  | sort -u
)

if (( ${#SAMPLES[@]} == 0 )); then
  echo "ERROR: No samples found in $KNEADDATA_DIR" >&2; exit 3
fi

echo "[$(date)] Dataset: $DATASET | ${#SAMPLES[@]} samples"
echo "[$(date)] DB: $INDEX | Workers: $MAX_WORKERS | Threads/worker: $THREADS"

load_miniforge
conda activate humann4

export MPA_DB_OCT22 INDEX THREADS PROF_DIR KNEADDATA_DIR SNIC_TMP SLURM_JOB_ID

for i in "${!SAMPLES[@]}"; do
  SAMPLE="${SAMPLES[$i]}"

  srun --exclusive -n 1 -N 1 -c "${THREADS}" bash -c '
    set -euo pipefail
    SAMPLE="'"${SAMPLE}"'"
    IDX="'"${i}"'"

    OUT_PROFILE="${PROF_DIR}/${SAMPLE}_profile.txt"
    if [[ -s "$OUT_PROFILE" ]]; then
      echo "[$(date)] Worker ${IDX}: profile exists, skipping ${SAMPLE}"
      exit 0
    fi

    R1="${KNEADDATA_DIR}/${SAMPLE}_paired_1.fastq.gz"
    R2="${KNEADDATA_DIR}/${SAMPLE}_paired_2.fastq.gz"
    if [[ ! -f "$R1" || ! -f "$R2" ]]; then
      echo "ERROR: Missing paired FASTQs for ${SAMPLE}" >&2; exit 10
    fi

    WORK="${SNIC_TMP}/mpa_oct22_${USER}/${SLURM_JOB_ID}/${SAMPLE}"
    mkdir -p "$WORK"/{in,tmp,out}
    export TMPDIR="$WORK/tmp"

    cp "$R1" "$R2" "$WORK/in/"

    echo "[$(date)] Worker ${IDX}: running MetaPhlAn on ${SAMPLE}"
    metaphlan \
      "$WORK/in/${SAMPLE}_paired_1.fastq.gz,$WORK/in/${SAMPLE}_paired_2.fastq.gz" \
      --input_type fastq \
      --bowtie2db   "$MPA_DB_OCT22" \
      -x            "$INDEX" \
      --offline \
      --no_map \
      --tmp_dir     "$WORK/tmp" \
      --nproc       "'"$THREADS"'" \
      -t rel_ab_w_read_stats \
      -o            "$WORK/out/${SAMPLE}_profile.txt"

    if [[ ! -s "$WORK/out/${SAMPLE}_profile.txt" ]]; then
      echo "ERROR: Profile missing for ${SAMPLE}" >&2
      ls -lah "$WORK/out" >&2 || true
      exit 20
    fi

    cp -f "$WORK/out/${SAMPLE}_profile.txt" "$OUT_PROFILE"
    rm -rf "$WORK"
    echo "[$(date)] Worker ${IDX}: completed ${SAMPLE}"
  ' &> "metalogs/mpa_oct22_${SLURM_JOB_ID}_${i}.log" &

  while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do
    sleep 2
  done
done

echo "[$(date)] All workers launched, waiting..."
wait

FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
  if [[ ! -s "${PROF_DIR}/${SAMPLE}_profile.txt" ]]; then
    echo "WARNING: ${SAMPLE} missing profile" >&2
    ((FAILED++))
  fi
done

echo "[$(date)] Done. Profiles: $PROF_DIR"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples incomplete   check metalogs/" >&2; exit 1
fi
echo "[$(date)] All ${#SAMPLES[@]} samples completed successfully"