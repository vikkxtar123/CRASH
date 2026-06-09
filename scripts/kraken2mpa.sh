#!/usr/bin/env bash
#SBATCH -A #allocation
#SBATCH -J kraken2mpa
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 1:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err
#SBATCH --mail-user=#email
#SBATCH --mail-type=END,FAIL

set -euo pipefail
mkdir -p logs

source /path/to/config.sh

# Usage: sbatch kraken2mpa.sh <dataset>
if [[ $# -lt 1 ]]; then
    echo "ERROR: Usage: sbatch kraken2mpa.sh <dataset>" >&2
    exit 1
fi
resolve_dataset "$1"
KRAKEN_DIR="${DATASET_DIR}/KRAKEN2"
OUT_DIR="${DATASET_DIR}/KRAKEN2_MPA"
mkdir -p "${OUT_DIR}/per_sample_mpa"

load_miniforge
conda activate krakentools

export KRAKEN_DIR OUT_DIR

echo "[$(date)] Dataset   : ${DATASET_DIR}"
echo "[$(date)] KRAKEN_DIR: ${KRAKEN_DIR}"
echo "[$(date)] OUT_DIR   : ${OUT_DIR}"

# ---- Approach 1: Bracken tables -------------------------------------
echo "[$(date)] Merging Bracken tables"

python3 << 'PYEOF'
import pandas as pd, os, glob

kraken_dir = os.environ["KRAKEN_DIR"]
out_dir    = os.environ["OUT_DIR"]

for level, ext in [("species", "bracken.S.txt"), ("genus", "bracken.G.txt")]:
    files = sorted(glob.glob(os.path.join(kraken_dir, f"*.{ext}")))
    print(f"Found {len(files)} Bracken {level} files")
    frames = {}
    for f in files:
        sample = os.path.basename(f).replace(f".{ext}", "")
        df = pd.read_csv(f, sep="\t").set_index("name")
        frames[sample] = df[["new_est_reads", "fraction_total_reads"]]

    counts = pd.DataFrame({s: frames[s]["new_est_reads"]          for s in frames}).fillna(0).astype(int)
    relab  = pd.DataFrame({s: frames[s]["fraction_total_reads"]   for s in frames}).fillna(0)
    counts.index.name = level
    relab.index.name  = level

    counts.to_csv(os.path.join(out_dir, f"bracken_{level}_counts.tsv"), sep="\t")
    relab.to_csv( os.path.join(out_dir, f"bracken_{level}_relab.tsv"),  sep="\t")
    print(f"  {level}: {counts.shape[0]} taxa x {counts.shape[1]} samples")

print("Bracken tables done.")
PYEOF

# ---- Approach 2: kreport2mpa ----------------------------------------
echo "[$(date)] Converting Kraken2 reports to MPA format"

CONVERTED=0
for kreport in "${KRAKEN_DIR}"/*.pluspf.report; do
    [[ -f "$kreport" ]] || continue
    sample=$(basename "$kreport" .pluspf.report)
    mpa_out="${OUT_DIR}/per_sample_mpa/${sample}.mpa.txt"

    kreport2mpa.py \
        -r "$kreport" \
        -o "$mpa_out" \
        --display-header \
        --no-intermediate-ranks && (( CONVERTED++ )) || echo "WARNING: Failed for ${sample}"
done
echo "[$(date)] Converted ${CONVERTED} files"

if (( CONVERTED > 0 )); then
    MPA_FILES=$(ls "${OUT_DIR}"/per_sample_mpa/*.mpa.txt | tr '\n' ' ')
    combine_mpa.py -i ${MPA_FILES} -o "${OUT_DIR}/merged_kraken_mpa_full.txt"

    python3 << 'PYEOF2'
import pandas as pd, os
out_dir = os.environ["OUT_DIR"]
df = pd.read_csv(os.path.join(out_dir, "merged_kraken_mpa_full.txt"), sep="\t", index_col=0)
species = df[df.index.str.contains(r"\|s__") & ~df.index.str.contains(r"\|t__")]
species.columns = [c.replace(".mpa", "") for c in species.columns]
species.to_csv(os.path.join(out_dir, "merged_kraken_mpa_species.tsv"), sep="\t")
print(f"MPA species table: {species.shape[0]} species x {species.shape[1]} samples")
PYEOF2
fi

echo "[$(date)] Done. Outputs in: ${OUT_DIR}"