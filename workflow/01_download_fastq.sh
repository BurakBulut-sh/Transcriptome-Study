#!/usr/bin/env bash
set -euo pipefail
source workflow/config.env

mkdir -p "${FASTQ_DIR}" "${LOG_DIR}"

RUNS="workflow/sra_run_ids.txt"
[[ -s "${RUNS}" ]] || { echo "Add run IDs to ${RUNS}"; exit 1; }

while read -r run; do
  [[ -z "${run}" || "${run}" =~ ^# ]] && continue
  echo "== ${run} =="

  "${SRA_TOOLKIT_BIN}/prefetch" "${run}"
  "${SRA_TOOLKIT_BIN}/fasterq-dump" "${run}" -O "${FASTQ_DIR}" --split-files

  # Optional cleanup of .sra if you want:
  # rm -f "${run}.sra"
done < "${RUNS}"
