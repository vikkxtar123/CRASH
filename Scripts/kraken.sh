#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH -J kraken2_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -o krakenlogs/%x_%j.out
#SBATCH -e krakenlogs/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p krakenlogs

# ---- Environment ----------------------------------------------------
module purge
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken2

# ---- Paths ----------------------------------------------------------
RAW_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KNEADDATA"
OUT_DIR="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705/KRAKEN2"
DB_DIR="/home/vikkxtar/lu2025-12-46/util/util_12-47/kraken"   # <-- must contain hash.k2d opts.k2d taxo.k2d + database*mers.kmer_distrib
THREADS="${SLURM_CPUS_PER_TASK:-8}"
WORKERS=6                   
READLEN=150

mkdir -p "$OUT_DIR"

echo "[$(date)] Starting Kraken2+Bracken task farm to process KNEADDATA samples"
echo "[$(date)] Using ${WORKERS} samples at a time with ${THREADS} cores each"
echo "[$(date)] RAW_DIR=$RAW_DIR"
echo "[$(date)] OUT_DIR=$OUT_DIR"
echo "[$(date)] DB_DIR=$DB_DIR"
echo "[$(date)] READLEN=$READLEN"

# Verify required executables exist
command -v kraken2 >/dev/null
command -v bracken >/dev/null

# Verify SNIC_TMP
if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2
  exit 1
fi

# ---- Build sample list ----------------------------------------------
mapfile -t SAMPLES < <(
  find "$RAW_DIR" -maxdepth 1 -type f -name 'SRR*_paired_1.fastq.gz' -printf '%f\n' \
  | sed -E 's/_paired_1\.fastq\.gz$//' \
  | sort -u
)

echo "[$(date)] Found ${#SAMPLES[@]} samples"

# ---- Optional: stage DB once per job to local disk -------------------
# This avoids each worker hammering the shared filesystem.
DB_LOCAL="${SNIC_TMP}/kraken_db_${SLURM_JOB_ID}"
if [[ ! -d "$DB_LOCAL" ]]; then
  echo "[$(date)] Staging Kraken2 DB to local disk: $DB_LOCAL"
  mkdir -p "$DB_LOCAL"

  # Copy minimal required Kraken2 DB files
  cp -f "$DB_DIR/hash.k2d" "$DB_DIR/opts.k2d" "$DB_DIR/taxo.k2d" "$DB_LOCAL/"

  # Bracken kmer distribution files (you have multiple; copy them)
  cp -f "$DB_DIR"/database*mers.kmer_distrib "$DB_LOCAL/" 2>/dev/null || true

  # Taxonomy folder is often needed for inspect/reporting; safe to include
  if [[ -d "$DB_DIR/taxonomy" ]]; then
    cp -a "$DB_DIR/taxonomy" "$DB_LOCAL/"
  fi

  # Any mapping files if present
  cp -f "$DB_DIR"/*.map "$DB_LOCAL/" 2>/dev/null || true
fi

# Validate DB local
for f in hash.k2d opts.k2d taxo.k2d; do
  [[ -f "${DB_LOCAL}/${f}" ]] || { echo "ERROR: Missing ${DB_LOCAL}/${f}" >&2; exit 2; }
done

# Export variables for worker script
export RAW_DIR OUT_DIR DB_LOCAL THREADS READLEN WORKERS

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

    # Output naming
    KRAKEN_OUT="${OUT_DIR}/${SAMPLE}.pluspf.out"
    KRAKEN_REP="${OUT_DIR}/${SAMPLE}.pluspf.report"
    BRACKEN_S="${OUT_DIR}/${SAMPLE}.bracken.S.txt"
    BRACKEN_G="${OUT_DIR}/${SAMPLE}.bracken.G.txt"

    # Skip if already finished
    if [[ -s "$KRAKEN_REP" && -s "$BRACKEN_S" && -s "$BRACKEN_G" ]]; then
      echo "[$(date)] Worker ${WORKER_ID}: Outputs exist, skipping ${SAMPLE}"
      exit 0
    fi

    # Worker-local scratch dir
    WORK="${SNIC_TMP}/kraken_worker_${SLURM_JOB_ID}_${WORKER_ID}"
    mkdir -p "$WORK"/{in,out}

    # Copy inputs locally (reduces shared FS load)
    cp "$R1" "$R2" "$WORK/in/"

    echo "[$(date)] Worker ${WORKER_ID}: Running Kraken2 for ${SAMPLE}"

    kraken2 \
      --db "$DB_LOCAL" \
      --memory-mapping \
      --threads "$THREADS" \
      --paired "$WORK/in/${SAMPLE}_paired_1.fastq.gz" "$WORK/in/${SAMPLE}_paired_2.fastq.gz" \
      --report "$WORK/out/${SAMPLE}.pluspf.report" \
      --output "$WORK/out/${SAMPLE}.pluspf.out"

    # Basic check
    if [[ ! -s "$WORK/out/${SAMPLE}.pluspf.report" ]]; then
      echo "ERROR: Missing Kraken2 report for ${SAMPLE}" >&2
      exit 11
    fi

    echo "[$(date)] Worker ${WORKER_ID}: Running Bracken for ${SAMPLE}"

    bracken \
      -d "$DB_LOCAL" \
      -i "$WORK/out/${SAMPLE}.pluspf.report" \
      -o "$WORK/out/${SAMPLE}.bracken.S.txt" \
      -r "$READLEN" -l S

    bracken \
      -d "$DB_LOCAL" \
      -i "$WORK/out/${SAMPLE}.pluspf.report" \
      -o "$WORK/out/${SAMPLE}.bracken.G.txt" \
      -r "$READLEN" -l G

    # Verify outputs
    for f in \
      "$WORK/out/${SAMPLE}.pluspf.out" \
      "$WORK/out/${SAMPLE}.pluspf.report" \
      "$WORK/out/${SAMPLE}.bracken.S.txt" \
      "$WORK/out/${SAMPLE}.bracken.G.txt"
    do
      [[ -s "$f" ]] || { echo "ERROR: Missing expected output $f" >&2; exit 12; }
    done

    echo "[$(date)] Worker ${WORKER_ID}: Copying outputs for ${SAMPLE}"
    cp -f "$WORK/out/${SAMPLE}.pluspf.out" "$OUT_DIR/"
    cp -f "$WORK/out/${SAMPLE}.pluspf.report" "$OUT_DIR/"
    cp -f "$WORK/out/${SAMPLE}.bracken.S.txt" "$OUT_DIR/"
    cp -f "$WORK/out/${SAMPLE}.bracken.G.txt" "$OUT_DIR/"

    rm -rf "$WORK"
    echo "[$(date)] Worker ${WORKER_ID}: Completed ${SAMPLE}"
  ' &> "krakenlogs/worker_${SLURM_JOB_ID}_${i}.log" &

  # Throttle: at most WORKERS concurrent workers
  while (( $(jobs -p | wc -l) >= WORKERS )); do
    sleep 2
  done
done

echo "[$(date)] All samples submitted, waiting for completion..."
wait

# ---- Final check ----------------------------------------------------
FAILED=0
for SAMPLE in "${SAMPLES[@]}"; do
  if [[ ! -s "$OUT_DIR/${SAMPLE}.pluspf.report" || ! -s "$OUT_DIR/${SAMPLE}.bracken.S.txt" ]]; then
    echo "WARNING: ${SAMPLE} failed or incomplete" >&2
    ((FAILED++))
  fi
done

echo "[$(date)] Task farm completed: processed ${#SAMPLES[@]} samples"
if (( FAILED > 0 )); then
  echo "WARNING: ${FAILED} samples failed - check worker logs in krakenlogs/" >&2
  exit 1
fi

echo "[$(date)] All samples completed successfully"
