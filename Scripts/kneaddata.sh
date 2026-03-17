#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH -J knead_test
#SBATCH -N 1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=12
#SBATCH -t 04:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

# ---- Environment ----------------------------------------------------
module purge
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kneaddata

# ---- Paths ----------------------------------------------------------
RAW_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/FASTQ"
OUT_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KNEADDATA"
KNEAD_DB="/home/vikkxtar/lu2025-12-46/Vik/Util/kneaddata_db"
TRIMMOMATIC_DIR="${CONDA_PREFIX}/share/trimmomatic-0.40-0/"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

SAMPLES=(SRR21311674 SRR21311675 SRR21311676)

# ---- Function -------------------------------------------------------
run_one() {
  local SAMPLE="$1"
  local R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
  local R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

  [[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing FASTQs for $SAMPLE"; return 2; }
  [[ -n "${SNIC_TMP:-}" ]] || { echo "ERROR: SNIC_TMP not set"; return 1; }

  local WORK="${SNIC_TMP}/kneaddata_${USER}/${SLURM_JOB_ID}/${SAMPLE}"
  mkdir -p "$WORK"/{in,tmp,out}

  echo "[$(date)] Processing $SAMPLE on $(hostname) with $THREADS threads"

  # Copy inputs to local disk
  cp "$R1" "$R2" "$WORK/in/"
  cd "$WORK"

  # ---- Run kneaddata ------------------------------------------------
  kneaddata \
    --input1 "$WORK/in/${SAMPLE}_1.fastq.gz" \
    --input2 "$WORK/in/${SAMPLE}_2.fastq.gz" \
    --output "$WORK/out" \
    --scratch "$WORK/tmp" \
    --output-prefix "$SAMPLE" \
    --reference-db "$KNEAD_DB" \
    --run-fastqc-start \
    --run-fastqc-end \
    --bowtie2-options="--very-fast" \
    --trimmomatic "$TRIMMOMATIC_DIR" \
    --max-memory 8g \
    -t "$THREADS" \
    -p 1 \
    --remove-intermediate-output \
    --log "$WORK/out/${SAMPLE}.log"

  # ---- Verify outputs before deleting anything ----------------------
  ls -lh \
    "$WORK/out/${SAMPLE}_paired_"*.fastq \
    "$WORK/out/${SAMPLE}_unmatched_"*.fastq \
    "$WORK/out/${SAMPLE}.log" >/dev/null 2>&1 \
    || { echo "ERROR: Missing expected outputs for $SAMPLE" >&2; return 3; }

  mkdir -p "$OUT_DIR"

  echo "[$(date)] Copying outputs for $SAMPLE"
  cp "$WORK/out/${SAMPLE}_paired_"*.fastq "$OUT_DIR/"
  cp "$WORK/out/${SAMPLE}_unmatched_"*.fastq "$OUT_DIR/"
  cp "$WORK/out/${SAMPLE}.log" "$OUT_DIR/"

  if [[ -d "$WORK/tmp/fastqc" ]]; then
    cp -r "$WORK/tmp/fastqc" "$OUT_DIR/fastqc_${SAMPLE}"
  fi

  {
    echo "kneaddata_${USER}/${SLURM_JOB_ID}/${SAMPLE}/out/${SAMPLE}_paired_*.fastq"
    echo "kneaddata_${USER}/${SLURM_JOB_ID}/${SAMPLE}/out/${SAMPLE}_unmatched_*.fastq"
    echo "kneaddata_${USER}/${SLURM_JOB_ID}/${SAMPLE}/out/${SAMPLE}.log"
    echo "kneaddata_${USER}/${SLURM_JOB_ID}/${SAMPLE}/tmp/fastqc/*"
  } >> "$SNIC_TMP/slurm_save_files"

  echo "[$(date)] Finished $SAMPLE"

}

export -f run_one
export RAW_DIR OUT_DIR KNEAD_DB TRIMMOMATIC_DIR THREADS

# ---- Parallel execution ---------------------------------------------
for s in "${SAMPLES[@]}"; do
  srun --exclusive -N1 -n1 -c "$THREADS" bash -lc "run_one $s" &
done

wait
echo "[$(date)] All samples completed"
