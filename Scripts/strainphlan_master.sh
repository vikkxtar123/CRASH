#!/usr/bin/env bash
#SBATCH -A lu2025-2-101
#SBATCH -J strainphlan_farm
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00
#SBATCH --mem=0
#SBATCH -o strain_farm/%x_%j.out
#SBATCH -e strain_farm/%x_%j.err
#SBATCH --mail-user=ja3363na-s@student.lu.se
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p strain_farm

# ---- Environment ----------------------------------------------------
module --force purge
module load SoftwareTree/Milan
module load Miniforge3/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate metaphlan

# ---- Paths ----------------------------------------------------------
export DIR2="/home/vikkxtar/lu2025-12-46/Vik/PRJNA813705"
export BASE_WORK="$DIR2/STRAINPHLAN"
export DB_DIR="/home/vikkxtar/lu2025-12-46/util/util_12-47/metaphlan_db_202403"
export INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"
export DB_PKL="$DB_DIR/${INDEX}.pkl"

export THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-6}"

# ---- New timestamped run folder -------------------------------------
export RUN_TAG="$(date +%Y%m%d)_${SLURM_JOB_ID}"
export RUN_DIR="$BASE_WORK/run_${RUN_TAG}"

export LOG_DIR="$RUN_DIR/logs"
export WORKER_LOG_DIR="$RUN_DIR/workerlogs"
export FAIL_LOG="$LOG_DIR/failures.log"

export CMARK_DIR="$BASE_WORK/consensus_markers"
export CLADEMARK_DIR="$RUN_DIR/CladeMarkers"
export OUT_DIR="$RUN_DIR/Output"

mkdir -p "$LOG_DIR" "$WORKER_LOG_DIR" "$CLADEMARK_DIR" "$OUT_DIR"
: > "$FAIL_LOG"

echo "[$(date)] RUN_DIR=$RUN_DIR"
echo "[$(date)] BASE_WORK=$BASE_WORK"
echo "[$(date)] CMARK_DIR=$CMARK_DIR"
echo "[$(date)] DB_PKL=$DB_PKL"
echo "[$(date)] Workers: $MAX_WORKERS | Threads/worker: $THREADS"

# ---- Verify SNIC_TMP ------------------------------------------------
if [[ -z "${SNIC_TMP:-}" ]]; then
  echo "ERROR: SNIC_TMP not set" >&2
  exit 1
fi
echo "[$(date)] SNIC_TMP=$SNIC_TMP"
df -h "$SNIC_TMP" || true

# ---- Verify inputs --------------------------------------------------
if [[ ! -s "$DB_PKL" ]]; then
  echo "ERROR: DB pkl not found: $DB_PKL" >&2
  exit 2
fi

if ! ls "$CMARK_DIR"/*.bz2 >/dev/null 2>&1; then
  echo "ERROR: No consensus markers found at: $CMARK_DIR/*.bz2" >&2
  exit 3
fi

cd "$RUN_DIR"

# --------------------------------------------------------------------
# 1) Extract clades present (once)
# --------------------------------------------------------------------
echo "[$(date)] Step 1: strainphlan --print_clades_only"
strainphlan \
  -s "$CMARK_DIR/"*.bz2 \
  -d "$DB_PKL" \
  --print_clades_only \
  -o "$RUN_DIR" \
  --nproc "$THREADS" \
  > "$LOG_DIR/clades.txt" 2> "$LOG_DIR/print_clades_only.err"

if [[ ! -s "$RUN_DIR/print_clades_only.tsv" ]]; then
  echo "ERROR: print_clades_only.tsv not created in $RUN_DIR" >&2
  exit 4
fi

# --------------------------------------------------------------------
# 2) Assign GTDB names to SGB clades
# --------------------------------------------------------------------
echo "[$(date)] Step 2: create output_sgb_names.tsv"

SGB2GTDB="/home/vikkxtar/.conda/envs/metaphlan/lib/python3.10/site-packages/metaphlan/utils/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB.tsv"
if [[ ! -f "$SGB2GTDB" ]]; then
  echo "ERROR: SGB2GTDB mapping file not found at: $SGB2GTDB" >&2
  exit 5
fi
cp -f "$SGB2GTDB" "$RUN_DIR/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB.tsv"

CLADES_FILE="$RUN_DIR/print_clades_only.tsv"
SGB_FILE="$RUN_DIR/mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB.tsv"
export OUTPUT_FILE="$RUN_DIR/output_sgb_names.tsv"

awk -F'\t' '
BEGIN { OFS="\t"; print "Clade","Name"; }
FNR==NR { map[$1]=$2; next; }
FNR==1  { next; }
{
  clade=$1;
  sgb=clade; sub(/^t__/,"",sgb);
  name=(sgb in map ? map[sgb] : "Not found");
  print clade, name;
}
' "$SGB_FILE" "$CLADES_FILE" > "$OUTPUT_FILE"

# Build clade list
tail -n +2 "$OUTPUT_FILE" | cut -f1 > "$RUN_DIR/clades_to_run.txt"
NCLADES=$(wc -l < "$RUN_DIR/clades_to_run.txt")
echo "[$(date)] Clades to run: $NCLADES"

# Export the worker script path so it is easy to change
WORKER_SCRIPT="/home/vikkxtar/lu2025-12-46/Vik/Scripts/strainphlan/strainphlan_worker.sh"
if [[ ! -x "$WORKER_SCRIPT" ]]; then
  echo "ERROR: worker script not found or not executable: $WORKER_SCRIPT" >&2
  echo "       Run: chmod +x strainphlan_worker.sh" >&2
  exit 6
fi

# --------------------------------------------------------------------
# 3) Farm: launch up to MAX_WORKERS workers at a time.
#    Track PIDs explicitly — wait -n does not work reliably with srun.
#    Before each launch, poll until a slot is free.
# --------------------------------------------------------------------
echo "[$(date)] Step 3: launching clade farm (max $MAX_WORKERS concurrent)"

declare -a PIDS=()

# Remove finished PIDs from the tracking array
reap_pids() {
  local alive=()
  for pid in "${PIDS[@]:-}"; do
    kill -0 "$pid" 2>/dev/null && alive+=("$pid")
  done
  PIDS=("${alive[@]:-}")
}

i=0
while read -r CLADE; do
  [[ -z "$CLADE" ]] && continue

  # Wait until a slot is free
  while true; do
    reap_pids
    (( ${#PIDS[@]} < MAX_WORKERS )) && break
    sleep 2
  done

  srun -Q --exclusive -n 1 -N 1 -c "$THREADS" \
    "$WORKER_SCRIPT" "$i" "$CLADE" \
    &> "$WORKER_LOG_DIR/worker_${SLURM_JOB_ID}_${i}.log" \
    < /dev/null &

  PIDS+=($!)
  sleep 1

  echo "[$(date)] Launched worker $i for ${CLADE} (active: ${#PIDS[@]}/$MAX_WORKERS)"
  i=$(( i + 1 )) || true
done < "$RUN_DIR/clades_to_run.txt"

echo "[$(date)] All workers launched, waiting for completion..."
wait "${PIDS[@]:-}"

# --------------------------------------------------------------------
# 4) Summary
# --------------------------------------------------------------------
FAIL_COUNT=$(wc -l < "$FAIL_LOG" | tr -d ' ')
echo "[$(date)] Failures logged: $FAIL_COUNT"
if (( FAIL_COUNT > 0 )); then
  echo "[$(date)] Failed clades:"
  cat "$FAIL_LOG"
fi

TREE_COUNT=$(find "$OUT_DIR" -type f \( -name "*.tre" -o -name "*.tree" -o -name "*.nwk" \) 2>/dev/null | wc -l | tr -d ' ')
echo "[$(date)] Trees produced: $TREE_COUNT / $NCLADES clades attempted"
echo "[$(date)] Done. Results in: $RUN_DIR"

if (( FAIL_COUNT > 0 )); then
  exit 1
fi