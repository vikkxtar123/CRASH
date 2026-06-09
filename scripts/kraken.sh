#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH -J kraken2_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

# Shared paths, dataset registry, and helper functions
source /path/to/config.sh

# ---- Dataset ---------------------------------------------------------
# Usage: sbatch kraken2.sh <dataset>
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch kraken2.sh <dataset>" >&2
    exit 1
fi
resolve_dataset "$1"
RAW_DIR="${DATASET_DIR}/KNEADDATA"
OUT_DIR="${DATASET_DIR}/KRAKEN2"
mkdir -p "$OUT_DIR"

# ---- Environment ----------------------------------------------------
load_miniforge
conda activate kraken2

THREADS="${SLURM_CPUS_PER_TASK:-8}"
WORKERS="${SLURM_NTASKS_PER_NODE:-6}"
READLEN=150

echo "[$(date)] Dataset  : ${DATASET_DIR}"
echo "[$(date)] RAW_DIR  : ${RAW_DIR}"
echo "[$(date)] OUT_DIR  : ${OUT_DIR}"
echo "[$(date)] DB       : ${KRAKEN_DB}"
echo "[$(date)] Workers  : ${WORKERS} | Threads/worker: ${THREADS}"

# ---- Pre-flight checks ----------------------------------------------
[[ -n "${SNIC_TMP:-}" ]] || { echo "ERROR: SNIC_TMP not set" >&2; exit 1; }
command -v kraken2 >/dev/null || { echo "ERROR: kraken2 not found" >&2; exit 1; }
command -v bracken >/dev/null || { echo "ERROR: bracken not found" >&2; exit 1; }

# ---- Build sample list ----------------------------------------------
# Handles both SRR* (NCBI) and CRR* (CRA) naming
mapfile -t SAMPLES < <(
    find "$RAW_DIR" -maxdepth 1 -type f -name '*_paired_1.fastq.gz' -printf '%f\n' \
    | sed -E 's/_paired_1\.fastq\.gz$//' \
    | sort -u
)
echo "[$(date)] Found ${#SAMPLES[@]} samples"
(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: No samples found in ${RAW_DIR}" >&2; exit 1; }

# ---- Stage DB to local scratch once per job -------------------------
DB_LOCAL="${SNIC_TMP}/kraken_db_${SLURM_JOB_ID}"
echo "[$(date)] Staging Kraken2 DB to local disk: ${DB_LOCAL}"
mkdir -p "$DB_LOCAL"
cp -f "${KRAKEN_DB}/hash.k2d" "${KRAKEN_DB}/opts.k2d" "${KRAKEN_DB}/taxo.k2d" "$DB_LOCAL/"
cp -f "${KRAKEN_DB}"/database*mers.kmer_distrib "$DB_LOCAL/" 2>/dev/null || true
[[ -d "${KRAKEN_DB}/taxonomy" ]] && cp -a "${KRAKEN_DB}/taxonomy" "$DB_LOCAL/"
cp -f "${KRAKEN_DB}"/*.map "$DB_LOCAL/" 2>/dev/null || true

for f in hash.k2d opts.k2d taxo.k2d; do
    [[ -f "${DB_LOCAL}/${f}" ]] || { echo "ERROR: Missing ${DB_LOCAL}/${f}" >&2; exit 2; }
done
echo "[$(date)] DB staged successfully"

export RAW_DIR OUT_DIR DB_LOCAL THREADS READLEN

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

        KRAKEN_REP="${OUT_DIR}/${SAMPLE}.pluspf.report"
        BRACKEN_S="${OUT_DIR}/${SAMPLE}.bracken.S.txt"
        BRACKEN_G="${OUT_DIR}/${SAMPLE}.bracken.G.txt"

        # Skip if already complete
        if [[ -s "$KRAKEN_REP" && -s "$BRACKEN_S" && -s "$BRACKEN_G" ]]; then
            echo "[$(date)] Worker ${WORKER_ID}: outputs exist, skipping ${SAMPLE}"
            exit 0
        fi

        WORK="${SNIC_TMP}/kraken_worker_${SLURM_JOB_ID}_${WORKER_ID}"
        mkdir -p "$WORK"/{in,out}
        cp "$R1" "$R2" "$WORK/in/"

        echo "[$(date)] Worker ${WORKER_ID}: Kraken2"
        kraken2 \
            --db "$DB_LOCAL" \
            --memory-mapping \
            --threads "$THREADS" \
            --paired \
            "$WORK/in/${SAMPLE}_paired_1.fastq.gz" \
            "$WORK/in/${SAMPLE}_paired_2.fastq.gz" \
            --report "$WORK/out/${SAMPLE}.pluspf.report" \
            --output "$WORK/out/${SAMPLE}.pluspf.out"

        [[ -s "$WORK/out/${SAMPLE}.pluspf.report" ]] \
            || { echo "ERROR: Missing Kraken2 report for ${SAMPLE}" >&2; exit 11; }

        echo "[$(date)] Worker ${WORKER_ID}: Bracken"
        bracken -d "$DB_LOCAL" \
            -i "$WORK/out/${SAMPLE}.pluspf.report" \
            -o "$WORK/out/${SAMPLE}.bracken.S.txt" \
            -r "$READLEN" -l S

        bracken -d "$DB_LOCAL" \
            -i "$WORK/out/${SAMPLE}.pluspf.report" \
            -o "$WORK/out/${SAMPLE}.bracken.G.txt" \
            -r "$READLEN" -l G

        for f in \
            "$WORK/out/${SAMPLE}.pluspf.out" \
            "$WORK/out/${SAMPLE}.pluspf.report" \
            "$WORK/out/${SAMPLE}.bracken.S.txt" \
            "$WORK/out/${SAMPLE}.bracken.G.txt"; do
            [[ -s "$f" ]] || { echo "ERROR: Missing output $f" >&2; exit 12; }
        done

        cp -f "$WORK/out/${SAMPLE}.pluspf.out"    "$OUT_DIR/"
        cp -f "$WORK/out/${SAMPLE}.pluspf.report" "$OUT_DIR/"
        cp -f "$WORK/out/${SAMPLE}.bracken.S.txt" "$OUT_DIR/"
        cp -f "$WORK/out/${SAMPLE}.bracken.G.txt" "$OUT_DIR/"

        rm -rf "$WORK"
        echo "[$(date)] Worker ${WORKER_ID}: Completed ${SAMPLE}"
    ' &> "logs/kraken2_worker_${SLURM_JOB_ID}_${i}.log" &

    while (( $(jobs -p | wc -l) >= WORKERS )); do sleep 2; done
done

echo "[$(date)] All workers launched, waiting..."
wait

# ---- Final check ----------------------------------------------------
FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
    if [[ ! -s "${OUT_DIR}/${SAMPLE}.pluspf.report" || \
          ! -s "${OUT_DIR}/${SAMPLE}.bracken.S.txt" ]]; then
        echo "WARNING: ${SAMPLE} incomplete" >&2
        (( FAILED++ ))
    fi
done

echo "[$(date)] Completed ${#SAMPLES[@]} samples"
(( FAILED == 0 )) || { echo "WARNING: ${FAILED} samples failed" >&2; exit 1; }
echo "[$(date)] All samples completed successfully"