#!/usr/bin/env bash
#SBATCH -A #alloc
#SBATCH -J kneaddata
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=6
#SBATCH -t 08:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

source /path/to/config.sh

# ---- Dataset ---------------------------------------------------------
# Usage: sbatch kneaddata.sh <dataset> [SRR1 SRR2 ...]
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch kneaddata.sh <dataset> [SRR1 SRR2 ...]" >&2
    exit 1
fi
resolve_dataset "$1"; shift
RAW_DIR="${DATASET_DIR}/FASTQ"
OUT_DIR="${DATASET_DIR}/KNEADDATA"

# ---- Naming convention (CRA uses different suffixes) -----------------
if [[ "$DATASET_DIR" == *CRA* ]]; then
    R1_SUFFIX="_f1.fq.gz"
    R2_SUFFIX="_r2.fq.gz"
    SAMPLE_GLOB="CRR*_f1.fq.gz"
else
    R1_SUFFIX="_1.fastq.gz"
    R2_SUFFIX="_2.fastq.gz"
    SAMPLE_GLOB="*_1.fastq.gz"
fi
export R1_SUFFIX R2_SUFFIX

# ---- Sample list -----------------------------------------------------
if [[ $# -gt 0 ]]; then
    SAMPLES=("$@")
else
    mapfile -t SAMPLES < <(
        find "$RAW_DIR" -name "$SAMPLE_GLOB" -printf "%f\n" \
        | sed "s/${R1_SUFFIX}//" \
        | sort -u
    )
    echo "Auto-discovered ${#SAMPLES[@]} samples from ${RAW_DIR}"
fi

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: No samples found in ${RAW_DIR}" >&2; exit 1; }

# ---- Environment ----------------------------------------------------
load_miniforge
conda activate kneaddata

TRIMMOMATIC_DIR="${CONDA_PREFIX}/share/trimmomatic-0.40-0/"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

# ---- Worker ---------------------------------------------------------
run_one() {
    local SAMPLE="$1"
    local R1="${RAW_DIR}/${SAMPLE}${R1_SUFFIX}"
    local R2="${RAW_DIR}/${SAMPLE}${R2_SUFFIX}"

    [[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing FASTQs for $SAMPLE" >&2; return 2; }

    make_scratch kneaddata "$SAMPLE"
    cd "$WORK"
    export TMPDIR="$WORK/tmp"

    echo "[$(date)] $SAMPLE on $(hostname)   ${THREADS} threads"

    # Normalise names inside scratch so kneaddata call is always the same
    cp "$R1" "$WORK/in/${SAMPLE}_1.fastq.gz"
    cp "$R2" "$WORK/in/${SAMPLE}_2.fastq.gz"

    kneaddata \
        --input1 "${WORK}/in/${SAMPLE}_1.fastq.gz" \
        --input2 "${WORK}/in/${SAMPLE}_2.fastq.gz" \
        --output "${WORK}/out" \
        --scratch "${WORK}/tmp" \
        --output-prefix "$SAMPLE" \
        --reference-db "$KNEAD_DB" \
        --run-fastqc-start \
        --run-fastqc-end \
        --fastqc-options "--outdir ${WORK}/out" \
        --bowtie2-options="--very-fast" \
        --trimmomatic "$TRIMMOMATIC_DIR" \
        --max-memory 20g \
        -t "$THREADS" -p 1 \
        --remove-intermediate-output \
        --log "${WORK}/out/${SAMPLE}.log"

    # Verify before compressing
    ls -lh \
        "${WORK}/out/${SAMPLE}_paired_"*.fastq \
        "${WORK}/out/${SAMPLE}_unmatched_"*.fastq \
        "${WORK}/out/${SAMPLE}.log" >/dev/null 2>&1 \
        || { echo "ERROR: Missing expected outputs for $SAMPLE" >&2; return 3; }

    # Compress on scratch before hitting Lustre
    echo "[$(date)] Compressing outputs for $SAMPLE"
    pigz -p "$THREADS" \
        "$WORK/out/${SAMPLE}_paired_"*.fastq \
        "$WORK/out/${SAMPLE}_unmatched_"*.fastq

    mkdir -p "$OUT_DIR"
    cp "${WORK}/out/${SAMPLE}_paired_"*.fastq.gz    "$OUT_DIR/"
    cp "${WORK}/out/${SAMPLE}_unmatched_"*.fastq.gz "$OUT_DIR/"
    cp "${WORK}/out/${SAMPLE}.log"                  "$OUT_DIR/"
    [[ -d "${WORK}/tmp/fastqc" ]] && cp -r "${WORK}/tmp/fastqc" "${OUT_DIR}/fastqc_${SAMPLE}"

    register_outputs kneaddata "$SAMPLE" \
        "out/${SAMPLE}_paired_*.fastq.gz" \
        "out/${SAMPLE}_unmatched_*.fastq.gz" \
        "out/${SAMPLE}.log" \
        "tmp/fastqc/*"

    echo "[$(date)] Finished $SAMPLE"
}
export -f run_one make_scratch register_outputs
export RAW_DIR OUT_DIR KNEAD_DB TRIMMOMATIC_DIR THREADS SNIC_TMP SLURM_JOB_ID USER

# ---- Run ------------------------------------------------------------
for s in "${SAMPLES[@]}"; do
    srun --exclusive -N1 -n1 -c "$THREADS" bash -lc "run_one $s" &
done
wait
echo "[$(date)] All samples completed"