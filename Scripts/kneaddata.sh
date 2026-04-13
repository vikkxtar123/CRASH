#!/usr/bin/env bash
#SBATCH -A <your_allocation>
#SBATCH -J kneaddata
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=6
#SBATCH -t 08:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=<your_email>
#SBATCH --mail-type=END,FAIL

set -euo pipefail                  # exit on error, unset variable, or pipe failure
mkdir -p logs

# Shared paths, dataset registry, and helper functions (make_scratch, register_outputs, etc.)
# Source shared config — use the absolute path when running on the cluster
source config.sh

# ---- Dataset ---------------------------------------------------------
# Usage: sbatch kneaddata.sh <dataset> [SRR1 SRR2 ...]
# <dataset> must match a key in config.sh DATASETS (e.g. leukemia, lymphoma, cra007433)
# Optional sample IDs override auto-discovery
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch kneaddata.sh <dataset> [SRR1 SRR2 ...]" >&2
    exit 1
fi
resolve_dataset "$1"; shift        # sets $DATASET_DIR, then drops the first argument
RAW_DIR="${DATASET_DIR}/FASTQ"
OUT_DIR="${DATASET_DIR}/KNEADDATA"

# ---- Naming convention -----------------------------------------------
# CRA007433 (Chinese SRA) uses _f1/_r2 suffixes; all other datasets use _1/_2
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
# Use explicit IDs if provided, otherwise discover from FASTQ directory
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

# ---- Environment -----------------------------------------------------
load_miniforge
conda activate kneaddata

TRIMMOMATIC_DIR="${CONDA_PREFIX}/share/trimmomatic-0.40-0/"
THREADS="${SLURM_CPUS_PER_TASK:-6}"

# ---- Worker function -------------------------------------------------
# Processes one sample: copy to scratch → kneaddata → compress → copy to Lustre
run_one() {
    local SAMPLE="$1"
    local R1="${RAW_DIR}/${SAMPLE}${R1_SUFFIX}"
    local R2="${RAW_DIR}/${SAMPLE}${R2_SUFFIX}"

    [[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing FASTQs for $SAMPLE" >&2; return 2; }

    # Create isolated scratch directory for this sample ($WORK set by make_scratch)
    make_scratch kneaddata "$SAMPLE"
    cd "$WORK"
    export TMPDIR="$WORK/tmp"      # redirect tool temp files to scratch

    echo "[$(date)] $SAMPLE on $(hostname) — ${THREADS} threads"

    # Rename inputs to standard _1/_2 convention so kneaddata call is dataset-agnostic
    cp "$R1" "$WORK/in/${SAMPLE}_1.fastq.gz"
    cp "$R2" "$WORK/in/${SAMPLE}_2.fastq.gz"

    # Quality trim (Trimmomatic), host decontamination (Bowtie2 vs human DB),
    # and FastQC reports at start and end; FastQC written to scratch to avoid
    # littering the submission directory
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

    # Verify expected outputs exist before committing to Lustre
    ls -lh \
        "${WORK}/out/${SAMPLE}_paired_"*.fastq \
        "${WORK}/out/${SAMPLE}_unmatched_"*.fastq \
        "${WORK}/out/${SAMPLE}.log" >/dev/null 2>&1 \
        || { echo "ERROR: Missing expected outputs for $SAMPLE" >&2; return 3; }

    # Compress on fast scratch before writing to Lustre (avoids quota pressure)
    echo "[$(date)] Compressing outputs for $SAMPLE"
    pigz -p "$THREADS" \
        "$WORK/out/${SAMPLE}_paired_"*.fastq \
        "$WORK/out/${SAMPLE}_unmatched_"*.fastq

    # Copy final outputs to permanent Lustre location
    mkdir -p "$OUT_DIR"
    cp "${WORK}/out/${SAMPLE}_paired_"*.fastq.gz    "$OUT_DIR/"
    cp "${WORK}/out/${SAMPLE}_unmatched_"*.fastq.gz "$OUT_DIR/"
    cp "${WORK}/out/${SAMPLE}.log"                  "$OUT_DIR/"
    [[ -d "${WORK}/tmp/fastqc" ]] && cp -r "${WORK}/tmp/fastqc" "${OUT_DIR}/fastqc_${SAMPLE}"

    # Register outputs with SLURM walltime safety net so files are preserved
    # even if the job hits its time limit before clean exit
    register_outputs kneaddata "$SAMPLE" \
        "out/${SAMPLE}_paired_*.fastq.gz" \
        "out/${SAMPLE}_unmatched_*.fastq.gz" \
        "out/${SAMPLE}.log" \
        "tmp/fastqc/*"

    echo "[$(date)] Finished $SAMPLE"
}

# Export function and variables so they are visible inside srun subshells
export -f run_one make_scratch register_outputs
export RAW_DIR OUT_DIR KNEAD_DB TRIMMOMATIC_DIR THREADS SNIC_TMP SLURM_JOB_ID USER

# ---- Run -------------------------------------------------------------
# Launch one srun worker per sample; --exclusive pins each to dedicated CPUs
# All workers run concurrently up to ntasks-per-node; wait collects exit codes
for s in "${SAMPLES[@]}"; do
    srun --exclusive -N1 -n1 -c "$THREADS" bash -lc "run_one $s" &
done
wait
echo "[$(date)] All samples completed"
