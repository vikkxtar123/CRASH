#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH -J deeparg_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 120:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#emailid
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

source /path/to/config.sh

# ---- Usage ----------------------------------------------------------
# sbatch deeparg_farm.sh samples_deeparg.txt
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch deeparg_farm.sh <sample_list.tsv>" >&2
    exit 1
fi
SAMPLE_LIST="$1"
[[ -f "$SAMPLE_LIST" ]] || { echo "ERROR: Sample list not found: $SAMPLE_LIST" >&2; exit 1; }

# ---- Environment ----------------------------------------------------
load_miniforge
conda activate deeparg_env

# ---- Paths ----------------------------------------------------------
DEEPARG_DB=/path/to/deeparg_db
OUTBASE=/DEEPARG_ALL
WORKER_SCRIPT="/path/to/deeparg_worker.sh"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-6}"

LOG_DIR="${OUTBASE}/worker_logs"
mkdir -p "$OUTBASE" "$LOG_DIR" logs

echo "[$(date)] Sample list : ${SAMPLE_LIST}"
echo "[$(date)] Output base : ${OUTBASE}"
echo "[$(date)] DB          : ${DEEPARG_DB}"
echo "[$(date)] Worker script: ${WORKER_SCRIPT}"
echo "[$(date)] Workers     : ${MAX_WORKERS}  |  Threads/worker: ${THREADS}"

# ---- Pre-flight checks ----------------------------------------------
[[ -n "${SNIC_TMP:-}" ]]     || { echo "ERROR: SNIC_TMP not set" >&2; exit 1; }
[[ -d "$DEEPARG_DB" ]]       || { echo "ERROR: deepARG DB not found: $DEEPARG_DB" >&2; exit 2; }
[[ -x "$WORKER_SCRIPT" ]]    || { echo "ERROR: Worker script not found/executable: $WORKER_SCRIPT" >&2; exit 3; }

# ---- Sample list ----------------------------------------------------
mapfile -t SAMPLE_LINES < "$SAMPLE_LIST"
echo "[$(date)] Total samples in list: ${#SAMPLE_LINES[@]}"

export DEEPARG_DB OUTBASE THREADS

# ---- Workers --------------------------------------------------------
WORKER_IDX=0
for LINE in "${SAMPLE_LINES[@]}"; do
    SAMPLE=$(echo "$LINE"  | cut -f1)
    KNEADDATA=$(echo "$LINE" | cut -f2)
    FINALDIR="${OUTBASE}/${SAMPLE}"

    # Skip if already complete
    if [[ -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.type"    && \
          -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.subtype" && \
          -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant" ]]; then
        echo "[skip] ${SAMPLE} already complete"
        continue
    fi

    srun -Q --exclusive -n1 -N1 -c "${THREADS}" \
        bash "${WORKER_SCRIPT}" \
            "${SAMPLE}" "${KNEADDATA}" "${WORKER_IDX}" \
            "${DEEPARG_DB}" "${OUTBASE}" \
        &> "${LOG_DIR}/worker_${WORKER_IDX}_${SAMPLE}.log" &

    (( WORKER_IDX++ )) || true
    while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do sleep 10; done
done

echo "[$(date)] All workers launched, waiting..."
wait

# ---- Final check ----------------------------------------------------
FAILED=0
for LINE in "${SAMPLE_LINES[@]}"; do
    SAMPLE=$(echo "$LINE" | cut -f1)
    FINALDIR="${OUTBASE}/${SAMPLE}"
    if [[ ! -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.type"    || \
          ! -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.subtype" || \
          ! -s "${FINALDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant" ]]; then
        echo "WARNING: ${SAMPLE} missing outputs   check ${LOG_DIR}" >&2
        (( FAILED++ )) || true
    fi
done

echo "[$(date)] Outputs: ${OUTBASE}"
(( FAILED == 0 )) || { echo "WARNING: ${FAILED} samples incomplete" >&2; exit 1; }
echo "[$(date)] All samples completed successfully"