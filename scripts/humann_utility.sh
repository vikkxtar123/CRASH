#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH -J humann4_post
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 04:00:00
#SBATCH --mem=192G
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL


set -euo pipefail
mkdir -p logs

DATASET="${1:?Usage: sbatch humann4_postprocess.sh <dataset>}"

source /path/tp/config.sh
resolve_dataset "$DATASET"
load_miniforge
conda activate humann4

HMN_DIR="${DATASET_DIR}/HUMANN4/humann_final"
OUT_DIR="${DATASET_DIR}/HUMANN4/postprocess"
FLAT_DIR="${OUT_DIR}/flat"
REGROUP_DIR="${OUT_DIR}/regrouped"
mkdir -p "$OUT_DIR" "$FLAT_DIR" "$REGROUP_DIR"

# Concurrent per-sample regrouping jobs (each ~2-4 GB).
# 4 is safe at 192G; bump to 6-8 if memory headroom allows.
N_PARALLEL=4

echo "[$(date)] Dataset: $DATASET"
echo "[$(date)] Input:   $HMN_DIR"
echo "[$(date)] Output:  $OUT_DIR"

# ---- Step 0a: Flatten per-sample subfolders ---------------------------
# Structure: humann_final/SAMPLE/SAMPLE_2_genefamilies.tsv etc.
# humann_join_tables needs all files in one flat directory.
echo "[$(date)] Flattening per-sample subfolders..."

for file_type in genefamilies pathabundance pathcoverage reactions; do
  find "$HMN_DIR" -mindepth 2 -maxdepth 2 -name "*_${file_type}.tsv" \
    -exec ln -sf {} "$FLAT_DIR/" \;
done

N_GENE=$(find "$FLAT_DIR" -name "*genefamilies.tsv"  | wc -l)
N_PATH=$(find "$FLAT_DIR" -name "*pathabundance.tsv" | wc -l)
echo "[$(date)] Flat dir: ${N_GENE} genefamilies, ${N_PATH} pathabundance files"

if (( N_GENE == 0 )); then
  echo "ERROR: No genefamilies files found under $HMN_DIR" >&2
  echo "Expected: $HMN_DIR/SAMPLE/SAMPLE_2_genefamilies.tsv" >&2
  exit 1
fi

# ---- Step 0b: Per-sample regroup to KO and EggNOG ---------------------
# Doing this BEFORE join keeps peak memory low: per-sample gene-family
# tables are small, while the joined UniRef90 matrix is huge.
echo "[$(date)] Per-sample regrouping to KEGG KOs and EggNOG (parallel: $N_PARALLEL)..."

regroup_one() {
  local gf="$1"
  local base
  base=$(basename "$gf" _genefamilies.tsv)
  if [[ ! -s "${REGROUP_DIR}/${base}_ko.tsv" ]]; then
    humann_regroup_table -i "$gf" \
      -o "${REGROUP_DIR}/${base}_ko.tsv" \
      -g uniref90_ko
  fi
  if [[ ! -s "${REGROUP_DIR}/${base}_eggnog.tsv" ]]; then
    humann_regroup_table -i "$gf" \
      -o "${REGROUP_DIR}/${base}_eggnog.tsv" \
      -g uniref90_eggnog
  fi
}
export -f regroup_one
export REGROUP_DIR

for gf in "$FLAT_DIR"/*_genefamilies.tsv; do
  regroup_one "$gf" &
  while (( $(jobs -rp | wc -l) >= N_PARALLEL )); do
    wait -n
  done
done
wait

N_KO=$(find "$REGROUP_DIR" -name "*_ko.tsv"     | wc -l)
N_EG=$(find "$REGROUP_DIR" -name "*_eggnog.tsv" | wc -l)
echo "[$(date)] Regrouped: ${N_KO} KO, ${N_EG} EggNOG files"

if (( N_KO != N_GENE )) || (( N_EG != N_GENE )); then
  echo "WARNING: regroup count mismatch (expected $N_GENE)" >&2
fi

# ---- Step 1: Join per-sample tables -----------------------------------
echo "[$(date)] Joining genefamilies..."
humann_join_tables \
  --input "$FLAT_DIR" \
  --output "${OUT_DIR}/genefamilies_joined.tsv" \
  --file_name genefamilies

echo "[$(date)] Joining pathabundance..."
humann_join_tables \
  --input "$FLAT_DIR" \
  --output "${OUT_DIR}/pathabundance_joined.tsv" \
  --file_name pathabundance

echo "[$(date)] Joining pathcoverage..."
humann_join_tables \
  --input "$FLAT_DIR" \
  --output "${OUT_DIR}/pathcoverage_joined.tsv" \
  --file_name pathcoverage

echo "[$(date)] Joining KO..."
humann_join_tables \
  --input "$REGROUP_DIR" \
  --output "${OUT_DIR}/ko_joined.tsv" \
  --file_name ko

echo "[$(date)] Joining EggNOG..."
humann_join_tables \
  --input "$REGROUP_DIR" \
  --output "${OUT_DIR}/eggnog_joined.tsv" \
  --file_name eggnog

# ---- Step 2: Normalise ------------------------------------------------
echo "[$(date)] Normalising genefamilies to CPM..."
humann_renorm_table \
  --input  "${OUT_DIR}/genefamilies_joined.tsv" \
  --output "${OUT_DIR}/genefamilies_cpm.tsv" \
  --units  cpm \
  --update-snames

echo "[$(date)] Normalising pathabundance to relative abundance..."
humann_renorm_table \
  --input  "${OUT_DIR}/pathabundance_joined.tsv" \
  --output "${OUT_DIR}/pathabundance_relab.tsv" \
  --units  relab \
  --update-snames

echo "[$(date)] Normalising KO to CPM..."
humann_renorm_table \
  --input  "${OUT_DIR}/ko_joined.tsv" \
  --output "${OUT_DIR}/ko_cpm.tsv" \
  --units  cpm \
  --update-snames

echo "[$(date)] Normalising EggNOG to CPM..."
humann_renorm_table \
  --input  "${OUT_DIR}/eggnog_joined.tsv" \
  --output "${OUT_DIR}/eggnog_cpm.tsv" \
  --units  cpm \
  --update-snames

# ---- Step 3: Split stratified vs unstratified -------------------------
for f in genefamilies_cpm pathabundance_relab ko_cpm eggnog_cpm; do
  echo "[$(date)] Splitting ${f}..."
  humann_split_stratified_table \
    --input  "${OUT_DIR}/${f}.tsv" \
    --output "$OUT_DIR"
done

# ---- Step 4: Clean sample names --------------------------------------
echo "[$(date)] Cleaning sample names..."
for f in \
  "${OUT_DIR}/genefamilies_cpm_unstratified.tsv" \
  "${OUT_DIR}/pathabundance_relab_unstratified.tsv" \
  "${OUT_DIR}/ko_cpm_unstratified.tsv" \
  "${OUT_DIR}/eggnog_cpm_unstratified.tsv" \
  "${OUT_DIR}/pathcoverage_joined.tsv"; do
  if [[ -s "$f" ]]; then
    sed -i '1s/_Abundance//g; 1s/_Coverage//g; 1s/_Abundance-RPKs//g' "$f"
    echo "  Cleaned: $(basename "$f")"
  fi
done

# ---- Step 5: Cleanup intermediate dirs --------------------------------
rm -rf "$FLAT_DIR" "$REGROUP_DIR"

# ---- Summary ----------------------------------------------------------
echo "[$(date)] Post-processing complete."
echo "[$(date)] Key output files:"
for f in \
  pathabundance_relab_unstratified.tsv \
  ko_cpm_unstratified.tsv \
  eggnog_cpm_unstratified.tsv \
  genefamilies_cpm_unstratified.tsv \
  pathcoverage_joined.tsv; do
  if [[ -s "${OUT_DIR}/${f}" ]]; then
    SIZE=$(du -sh "${OUT_DIR}/${f}" | cut -f1)
    ROWS=$(wc -l < "${OUT_DIR}/${f}")
    echo "  ${f}: ${SIZE}, ${ROWS} rows"
  else
    echo "  WARNING: ${f} missing or empty"
  fi
done
  "${OUT_DIR}/pathabundance_relab_unstratified.tsv" \
  "${OUT_DIR}/ko_cpm_unstratified.tsv" \
  "${OUT_DIR}/eggnog_cpm_unstratified.tsv" \
  "${OUT_DIR}/pathcoverage_joined.tsv"; do
  if [[ -s "$f" ]]; then
    sed -i '1s/_Abundance//g; 1s/_Coverage//g; 1s/_Abundance-RPKs//g' "$f"
    echo "  Cleaned: $(basename $f)"
  fi
done

# ---- Step 6: Cleanup symlinks -----------------------------------------
rm -rf "$FLAT_DIR"

# ---- Summary ----------------------------------------------------------
echo "[$(date)] Post-processing complete."
echo "[$(date)] Key output files:"
for f in \
  pathabundance_relab_unstratified.tsv \
  ko_cpm_unstratified.tsv \
  eggnog_cpm_unstratified.tsv \
  genefamilies_cpm_unstratified.tsv \
  pathcoverage_joined.tsv; do
  if [[ -s "${OUT_DIR}/${f}" ]]; then
    SIZE=$(du -sh "${OUT_DIR}/${f}" | cut -f1)
    ROWS=$(wc -l < "${OUT_DIR}/${f}")
    echo "  ${f}: ${SIZE}, ${ROWS} rows"
  else
    echo "  WARNING: ${f} missing or empty"
  fi
done