# config.sh — source at the top of every script

# ---- SLURM -----------------------------------------------------------
export ACCOUNT="" #Account partiton
export MAIL_USER="" #Fill in email ID to get SLURM notifications on job completion or failure

# ---- Known datasets --------------------------------------------------
# Usage in scripts: resolve_dataset <name>  ?  sets $DATASET_DIR
declare -A DATASETS=(
    [prjna813705]="/PRJNA813705"
    [leukemia]="/KuleckaLeukemia"
    [lymphoma]="/KuleckaLymphoma"
    [cra007433]="/CRA007433"
)

resolve_dataset() {
    local key="${1,,}"   # lowercase
    DATASET_DIR="${DATASETS[$key]:-}"
    if [[ -z "$DATASET_DIR" ]]; then
        echo "ERROR: Unknown dataset '${1}'. Valid keys: ${!DATASETS[*]}" >&2
        return 1
    fi
    export DATASET_DIR
}

# ---- Shared utility paths --------------------------------------------
export UTIL="/path/to/utility/folder"
export KNEAD_DB="${UTIL}/kneaddata_db"
export MPA_DB="/path/to/metaphlan_db_202403"
export MPA_PKL="${MPA_DB}/mpa_vJun23_CHOCOPhlAnSGB_202403.pkl"
export KRAKEN_DB="/path/to/kraken"

# ---- Module helpers --------------------------------------------------
load_miniforge() {
    module purge
    module load Miniforge3/25.3.0-3
    source "$(conda info --base)/etc/profile.d/conda.sh"
}

load_parallel() {
    module purge
    module load GCCcore/12.3.0
    module load parallel/20230722
}

# ---- Scratch helper --------------------------------------------------
# Sets $WORK in the caller's scope
make_scratch() {
    local tool="$1" sample="$2"
    [[ -n "${SNIC_TMP:-}" ]] || { echo "ERROR: SNIC_TMP not set" >&2; return 1; }
    WORK="${SNIC_TMP}/${tool}_${USER}/${SLURM_JOB_ID}/${sample}"
    mkdir -p "${WORK}"/{in,tmp,out}
    export WORK
}

# ---- slurm_save_files registration -----------------------------------
register_outputs() {
    local tool="$1" sample="$2"; shift 2
    for pattern in "$@"; do
        echo "${tool}_${USER}/${SLURM_JOB_ID}/${sample}/${pattern}" \
            >> "${SNIC_TMP}/slurm_save_files"
    done
}
