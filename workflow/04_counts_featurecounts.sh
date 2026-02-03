#!/usr/bin/env bash
set -euo pipefail
source workflow/config.env

mkdir -p "${COUNTS_DIR}" "${LOG_DIR}"
[[ -n "${MOD_SUBREAD}" ]] && module load "${MOD_SUBREAD}" || true

# Produces one matrix across all BAMs (as in your notes).
out="${COUNTS_DIR}/featurecounts_matrix.tsv"

featureCounts -a "${GTF}" -g geneid -o "${out}" "${BAM_DIR}"/*.sorted.bam -T "${THREADS}" -F GTF
