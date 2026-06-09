#!/usr/bin/env bash
set -euo pipefail

SAMPLE="$1"
KNEADDATA="$2"
WORKER_ID="$3"
DEEPARG_DB="$4"
OUTBASE="$5"

source $(conda info --base)/etc/profile.d/conda.sh
conda activate deeparg_env

echo "[$(date)] Worker ${WORKER_ID}: ${SAMPLE} on $(hostname)"

R1="${KNEADDATA}/${SAMPLE}_paired_1.fastq.gz"
R2="${KNEADDATA}/${SAMPLE}_paired_2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: Missing FASTQs for ${SAMPLE}" >&2; exit 10; }

WORKDIR="${SNIC_TMP}/deeparg_${SLURM_JOB_ID}_${WORKER_ID}_${SAMPLE}"
mkdir -p "$WORKDIR"
export TMPDIR="$WORKDIR"

cp "$R1" "$R2" "$WORKDIR/"

echo "[$(date)] Worker ${WORKER_ID}: running deepARG"
deeparg short_reads_pipeline \
    --forward_pe_file "${WORKDIR}/${SAMPLE}_paired_1.fastq.gz" \
    --reverse_pe_file "${WORKDIR}/${SAMPLE}_paired_2.fastq.gz" \
    --output_file "${WORKDIR}/${SAMPLE}" \
    -d "${DEEPARG_DB}" \
    --bowtie_16s_identity 100

[[ -s "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.type"    && \
   -s "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.subtype" && \
   -s "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant" ]] \
    || { echo "ERROR: deepARG outputs missing for ${SAMPLE}" >&2; ls -lah "$WORKDIR" >&2; exit 20; }

FINALDIR="${OUTBASE}/${SAMPLE}"
mkdir -p "$FINALDIR"
cp "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant"         "$FINALDIR/"
cp "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.subtype" "$FINALDIR/"
cp "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG.merged.quant.type"    "$FINALDIR/"
cp "${WORKDIR}/${SAMPLE}.clean.deeparg.mapping.ARG"                      "$FINALDIR/"

rm -rf "$WORKDIR"
echo "[$(date)] Worker ${WORKER_ID}: Completed ${SAMPLE}"