#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH -J mpa_s2m_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

source /path/to/config.sh

# ---- Dataset ---------------------------------------------------------
# Usage: sbatch metaphlan.sh <dataset> [SRR1 SRR2 ...]
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch metaphlan.sh <dataset> [SRR1 SRR2 ...]" >&2
    exit 1
fi
resolve_dataset "$1"; shift
RAW_DIR="${DATASET_DIR}/KNEADDATA"
BASE_OUT="${DATASET_DIR}/METAPHLAN"

# ---- Environment ----------------------------------------------------
load_miniforge
conda activate metaphlan

# ---- Paths ----------------------------------------------------------
INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
DB_DIR="${MPA_DB}"
DB_PKL="${DB_DIR}/${INDEX}.pkl"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-6}"

# ---- New run folder -------------------------------------------------
RUN_TAG="$(date +%Y%m%d)_${SLURM_JOB_ID}"
OUT_DIR="${BASE_OUT}/run_${RUN_TAG}"
PROF_DIR="${OUT_DIR}/profiles"
BT2_DIR="${OUT_DIR}/bowtie2"
CONS_DIR="${OUT_DIR}/consensus_markers"
LOG_DIR="${OUT_DIR}/worker_logs"

mkdir -p "$PROF_DIR" "$BT2_DIR" "$CONS_DIR" "$LOG_DIR"

echo "[$(date)] Dataset  : ${DATASET_DIR}"
echo "[$(date)] Outputs  : ${OUT_DIR}"
echo "[$(date)] DB       : ${DB_DIR}"
echo "[$(date)] Workers  : ${MAX_WORKERS}  |  Threads/worker: ${THREADS}"

# ---- Pre-flight checks ----------------------------------------------
[[ -n "${SNIC_TMP:-}" ]] || { echo "ERROR: SNIC_TMP not set" >&2; exit 1; }
[[ -s "$DB_PKL" ]]       || { echo "ERROR: DB pkl not found: $DB_PKL" >&2; ls -lh "$DB_DIR" | head >&2; exit 2; }

# ---- Sample list -----------------------------------------------------
if [[ $# -gt 0 ]]; then
    SAMPLES=("$@")
else
    mapfile -t SAMPLES < <(
        find "$RAW_DIR" -maxdepth 1 -type f -name '*_paired_1.fastq.gz' -printf '%f\n' \
        | sed -E 's/_paired_1\.fastq\.gz$//' \
        | sort -u \
        | while read -r s; do
            if [[ -s "${PROF_DIR}/${s}_profile.txt" && \
                  -s "${BT2_DIR}/${s}.bowtie2.bz2"  && \
                  -s "${CONS_DIR}/${s}.json.bz2" ]]; then
                echo "[skip] ${s} already complete" >&2
            else
                echo "$s"
            fi
        done
    )
    echo "[$(date)] Found ${#SAMPLES[@]} samples needing processing"
fi

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: No samples found in ${RAW_DIR}" >&2; exit 3; }

export RAW_DIR DB_DIR INDEX DB_PKL THREADS PROF_DIR BT2_DIR CONS_DIR LOG_DIR

# ---- Workers --------------------------------------------------------
for i in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$i]}"

    srun -Q --exclusive -n1 -N1 -c "${THREADS}" bash -c '
        set -euo pipefail
        SAMPLE="'"${SAMPLE}"'"
        WORKER_ID="'"${i}"'"

        echo "[$(date)] Worker ${WORKER_ID}: ${SAMPLE} on $(hostname)"

        R1="${RAW_DIR}/${SAMPLE}_paired_1.fastq.gz"
        R2="${RAW_DIR}/${SAMPLE}_paired_2.fastq.gz"
        [[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing FASTQs for ${SAMPLE}" >&2; exit 10; }

        WORK="${SNIC_TMP}/mpa_s2m_${SLURM_JOB_ID}_${WORKER_ID}_${SAMPLE}"
        mkdir -p "$WORK"/{in,tmp,out}
        cd "$WORK"
        export TMPDIR="$WORK/tmp"

        cp "$R1" "$R2" "$WORK/in/"

        echo "[$(date)] Worker ${WORKER_ID}: MetaPhlAn"
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

        [[ -s "$WORK/out/${SAMPLE}_profile.txt" && \
           -s "$WORK/out/${SAMPLE}.bowtie2.bz2" && \
           -s "$WORK/out/${SAMPLE}.sam.bz2" ]] \
            || { echo "ERROR: MetaPhlAn outputs missing for ${SAMPLE}" >&2; ls -lah "$WORK/out" >&2; exit 20; }

        cp -f "$WORK/out/${SAMPLE}_profile.txt"  "'"$PROF_DIR"'/"
        cp -f "$WORK/out/${SAMPLE}.bowtie2.bz2"  "'"$BT2_DIR"'/"

        echo "[$(date)] Worker ${WORKER_ID}: sample2markers"
        sample2markers.py \
            -d "$DB_PKL" \
            -i "$WORK/out/${SAMPLE}.sam.bz2" \
            -o "$WORK/out" \
            -f bz2 \
            -n "'"$THREADS"'" \
            --tmp "$WORK/tmp"

        [[ -s "$WORK/out/${SAMPLE}.json.bz2" ]] \
            || { echo "ERROR: sample2markers missing output for ${SAMPLE}" >&2; ls -lah "$WORK/out" >&2; exit 30; }

        cp -f "$WORK/out/${SAMPLE}.json.bz2" "'"$CONS_DIR"'/"
        rm -rf "$WORK"
        echo "[$(date)] Worker ${WORKER_ID}: Completed ${SAMPLE}"
    ' &> "${LOG_DIR}/worker_${i}_${SAMPLE}.log" &

    while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do sleep 2; done
done

echo "[$(date)] All workers launched, waiting..."
wait

# ---- Final check ----------------------------------------------------
FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ ! -s "${PROF_DIR}/${SAMPLE}_profile.txt" || \
          ! -s "${BT2_DIR}/${SAMPLE}.bowtie2.bz2"  || \
          ! -s "${CONS_DIR}/${SAMPLE}.json.bz2" ]]; then
        echo "WARNING: ${SAMPLE} missing outputs   check ${LOG_DIR}" >&2
        (( FAILED++ ))
    fi
done

echo "[$(date)] Outputs: ${OUT_DIR}"
(( FAILED == 0 )) || { echo "WARNING: ${FAILED} samples incomplete" >&2; exit 1; }
echo "[$(date)] All samples completed successfully"