#!/usr/bin/env bash
# strainphlan_worker.sh
#
# Called by strainphlan_master.sh via:
#   srun -Q --exclusive -n 1 -N 1 -c $THREADS strainphlan_worker.sh <WID> <CLADE>
#
# All paths come from environment variables exported by the master:
#   LOG_DIR, FAIL_LOG, CMARK_DIR, CLADEMARK_DIR, OUT_DIR, DB_PKL, THREADS,
#   SNIC_TMP, SLURM_JOB_ID

set -euo pipefail

# ---- Arguments ------------------------------------------------------
WID="${1:?worker ID not supplied}"
CLADE="${2:?clade not supplied}"

echo "[$(date)] Worker ${WID}: starting ${CLADE} on $(hostname)"

# ---- Scratch dir on fast local disk ---------------------------------
WDIR="${SNIC_TMP}/strain_${SLURM_JOB_ID}_${WID}"
mkdir -p "$WDIR"
export TMPDIR="$WDIR"

# ---- Failure trap: log and clean up ---------------------------------
trap '{
  rc=$?
  echo "FAILED clade=${CLADE} wid=${WID} exit=${rc} time=$(date)" >> "$FAIL_LOG"
  rm -rf "$WDIR"
  exit $rc
}' ERR

# ---- Step A: extract markers for this clade -------------------------
echo "[$(date)] Worker ${WID}: extracting markers for ${CLADE}"
extract_markers.py \
  -c "$CLADE" \
  -d "$DB_PKL" \
  -o "$WDIR" \
  > "${LOG_DIR}/${CLADE}.extract.out" \
  2> "${LOG_DIR}/${CLADE}.extract.err"

mv "$WDIR/${CLADE}.fna" "$CLADEMARK_DIR/${CLADE}.fna"

# ---- Step B: run StrainPhlAn for this clade -------------------------
mkdir -p "${OUT_DIR}/${CLADE}"

echo "[$(date)] Worker ${WID}: running strainphlan for ${CLADE}"
strainphlan \
  -s "${CMARK_DIR}/"*.bz2 \
  -m "${CLADEMARK_DIR}/${CLADE}.fna" \
  -o "${OUT_DIR}/${CLADE}" \
  -c "$CLADE" \
  -d "$DB_PKL" \
  --nproc "$THREADS" \
  > "${LOG_DIR}/${CLADE}.strainphlan.out" \
  2> "${LOG_DIR}/${CLADE}.strainphlan.err"

# ---- Clean up scratch -----------------------------------------------
rm -rf "$WDIR"

echo "[$(date)] Worker ${WID}: completed ${CLADE}"