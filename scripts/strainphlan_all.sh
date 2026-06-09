#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH --requeue
#SBATCH -J strainphlan_all
#SBATCH -N 1
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=8
#SBATCH -t 72:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs metalogs

source /path/to/config.sh
load_miniforge
conda activate metaphlan

# ---- Consensus marker directories ------------------------------------
PRJNA_MARKERS="/PRJNA813705/STRAINPHLAN/consensus_markers"
LEUK_MARKERS="/Kulecka/Leukemia/METAPHLAN/consensus_markers"
LYMP_MARKERS="/Kulecka/Lymphoma/METAPHLAN/consensus_markers"
CRA_MARKERS="/CRA007433/METAPHLAN/consensus_markers"

# ---- Output directories ----------------------------------------------
OUT_BASE="/STRAINPHLAN_ALL"
CLADE_MARKER_DIR="${OUT_BASE}/clade_markers"
TREE_DIR="${OUT_BASE}/trees"
mkdir -p "$CLADE_MARKER_DIR" "$TREE_DIR"

# ---- DB --------------------------------------------------------------
DB_PKL="$MPA_PKL"
INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MAX_WORKERS="${SLURM_NTASKS_PER_NODE:-6}"

echo "[$(date)] Collecting all consensus markers..."

# Collect all json.bz2 files from all four datasets
mapfile -t ALL_MARKERS < <(
  find "$PRJNA_MARKERS" "$LEUK_MARKERS" "$LYMP_MARKERS" "$CRA_MARKERS" \
    -maxdepth 1 -name "*.json.bz2" | sort
)

echo "[$(date)] Total markers: ${#ALL_MARKERS[@]}"
if (( ${#ALL_MARKERS[@]} != 334 )); then
  echo "WARNING: Expected 334 markers, got ${#ALL_MARKERS[@]}" >&2
fi

[[ -n "${SNIC_TMP:-}" ]] || { echo "ERROR: SNIC_TMP not set" >&2; exit 1; }
[[ -s "$DB_PKL" ]]       || { echo "ERROR: DB pkl not found: $DB_PKL" >&2; exit 1; }

# ---- Extract clade list from merged profiles -------------------------
echo "[$(date)] Building clade list from merged profiles..."

CLADE_LIST="${OUT_BASE}/clade_list.txt"

python3 << PYEOF
import re, os

profiles = [
    "/PRJNA813705/merged_abundance.txt",
    "/Kulecka/Leukemia/leukemia_merged_profiles.tsv",
    "/Kulecka/Lymphoma/lymphoma_merged_profiles.tsv",
    "/CRA007433/cra_merged_profiles.tsv",
]

clades = set()
for pf in profiles:
    if not os.path.exists(pf):
        print(f"WARNING: profile not found: {pf}")
        continue
    with open(pf) as f:
        for line in f:
            if line.startswith("#"): continue
            clade = line.strip().split("\t")[0]
            # Only strain-level t__ entries
            if "|t__" in clade:
                clades.add(clade)

os.makedirs("${OUT_BASE}", exist_ok=True)
with open("${CLADE_LIST}", "w") as f:
    for c in sorted(clades):
        f.write(c + "\n")
print(f"Written {len(clades)} clades to ${CLADE_LIST}")
PYEOF

mapfile -t CLADES < "$CLADE_LIST"
echo "[$(date)] Clades to process: ${#CLADES[@]}"

# ---- Extract clade markers -------------------------------------------
echo "[$(date)] Extracting clade markers (this may take a while)..."

for CLADE in "${CLADES[@]}"; do
  MARKER_FILE="${CLADE_MARKER_DIR}/${CLADE}.fna"
  [[ -s "$MARKER_FILE" ]] && continue

  extract_markers.py \
    --clade      "$CLADE" \
    --output_dir "$CLADE_MARKER_DIR" \
    -d           "$DB_PKL" \
    2>/dev/null || true
done

N_MARKERS=$(find "$CLADE_MARKER_DIR" -name "*.fna" | wc -l)
echo "[$(date)] Clade marker files: $N_MARKERS"

# ---- Run StrainPhlAn per clade ---------------------------------------
echo "[$(date)] Starting StrainPhlAn farm..."

export CLADE_MARKER_DIR TREE_DIR THREADS DB_PKL SNIC_TMP SLURM_JOB_ID
ALL_MARKERS_STR="${ALL_MARKERS[*]}"
export ALL_MARKERS_STR

for i in "${!CLADES[@]}"; do
  CLADE="${CLADES[$i]}"

  srun -Q --exclusive -n 1 -N 1 -c "${THREADS}" bash -c '
    set -euo pipefail
    CLADE="'"${CLADE}"'"
    IDX="'"${i}"'"

    CLADE_OUT="${TREE_DIR}/${CLADE}"

    # Skip if tree already exists
    if ls "${CLADE_OUT}"/*.tre &>/dev/null 2>&1; then
      echo "[$(date)] Worker ${IDX}: ${CLADE} already done, skipping"
      exit 0
    fi

    MARKER="${CLADE_MARKER_DIR}/${CLADE}.fna"
    if [[ ! -s "$MARKER" ]]; then
      echo "[$(date)] Worker ${IDX}: No marker file for ${CLADE}, skipping"
      exit 0
    fi

    mkdir -p "$CLADE_OUT"

    WORK="${SNIC_TMP}/stp_${SLURM_JOB_ID}_${IDX}"
    mkdir -p "$WORK"
    export TMPDIR="$WORK"

    echo "[$(date)] Worker ${IDX}: StrainPhlAn for ${CLADE}"

    strainphlan \
      --samples        '"${ALL_MARKERS_STR}"' \
      --markers        "$MARKER" \
      --output_dir     "$CLADE_OUT" \
      --clade          "$CLADE" \
      --phylophlan_mode fast \
      --nprocs         "'"$THREADS"'" \
      -d               "'"$DB_PKL"'" \
      2>&1 | tee "metalogs/strainphlan_${SLURM_JOB_ID}_${IDX}.log"

    rm -rf "$WORK"
    echo "[$(date)] Worker ${IDX}: ${CLADE} complete"
  ' &> "metalogs/strainphlan_out_${SLURM_JOB_ID}_${i}.log" &

  while (( $(jobs -p | wc -l) >= MAX_WORKERS )); do sleep 5; done
done

echo "[$(date)] All workers launched, waiting..."
wait

# ---- Summary ---------------------------------------------------------
N_TREES=$(find "$TREE_DIR" -name "*.tre" | wc -l)
echo "[$(date)] Trees generated: $N_TREES / ${#CLADES[@]} clades attempted"
echo "[$(date)] Output: $TREE_DIR"