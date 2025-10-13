#!/usr/bin/env bash
set -Eeuo pipefail

# ------------------------------------------------------------
# Batch Correction Runner (test harness)
# Mirrors the new CLI used by module_batch_correction_new_refactored.R
# ------------------------------------------------------------
# Default layout:
#   PROJECT_ROOT -> ".../PreprocessingPipeline"
#   DATA_DIR     -> ".../04_NORMALIZATION"
#   OUT_DIR      -> ".../05_CORRECTION"
#   METADATA     -> ".../Metadata.tsv"
#   PARAMS       -> ".../Parameters/05_CORRECTION_parameters copy.json"
#   IDS          -> ".../04_NORMALIZATION/id_mapping.tsv" (optional)
#   COUNTER      -> auto-detected (05 if data_04.tsv exists, 06 if data_05.tsv exists)
# ------------------------------------------------------------

# Helper loader path for R (portable if you export this on other machines)
export PIPELINE_HELPERS_R="${PIPELINE_HELPERS_R:-/Users/rozsurmiz/Documents/GitHub/Preprocessing/modules_R/helper_functions_R.R}"

: "${PROJECT_ROOT:="/Users/rozsurmiz/Desktop/Projects/PreprocessingPipeline/PreprocessingPipeline"}"
: "${METHOD:="limma"}"  # "limma", "combat" or "harmony"
: "${DATA_DIR:="${PROJECT_ROOT}/04_NORMALIZATION"}"
: "${OUT_DIR:="${PROJECT_ROOT}/05_CORRECTION"}"
: "${METADATA:="${PROJECT_ROOT}/Metadata.tsv"}"
: "${PARAMS:="${PROJECT_ROOT}/Parameters/05_CORRECTION_parameters.${METHOD}.json"}" #path/to/05_CORRECTION_parameters.combat.json" #/path/to/05_CORRECTION_parameters.harmony.json"
: "${IDS:="${PROJECT_ROOT}/04_NORMALIZATION/id_mapping.tsv"}"   # <-- new
: "${LOG:="${PROJECT_ROOT}/bash_log_rerun_05_correction.txt"}"
: "${MODULE:="/Users/rozsurmiz/Documents/GitHub/Preprocessing/modules_R/module_batch_correction.r"}"
: "${RSCRIPT_BIN:="Rscript"}"

# ---- Auto-detect COUNTER ----
COUNTER=""
if [[ -f "${DATA_DIR}/data_05.tsv" ]]; then
  COUNTER="06"
elif [[ -f "${DATA_DIR}/data_04.tsv" ]]; then
  COUNTER="05"
elif [[ -f "${DATA_DIR}/data_03.tsv" ]]; then
  COUNTER="04"
else
  echo "ERROR: Could not auto-detect COUNTER in ${DATA_DIR}" >&2
  exit 1
fi

# Manual override possible
if [[ -n "${COUNTER_OVERRIDE:-}" ]]; then
  COUNTER="${COUNTER_OVERRIDE}"
fi

INPUT_FILE="${DATA_DIR}/data_$(printf "%02d" $((10#${COUNTER} - 1))).tsv"
OUTPUT_FILE="${OUT_DIR}/data_${METHOD}_${COUNTER}.tsv"

# ---- Preflight checks ----
command -v "${RSCRIPT_BIN}" >/dev/null 2>&1 || { echo "ERROR: ${RSCRIPT_BIN} not found" >&2; exit 1; }
[[ -f "${MODULE}" ]]   || { echo "ERROR: MODULE not found: ${MODULE}" >&2; exit 1; }
[[ -d "${DATA_DIR}" ]] || { echo "ERROR: DATA_DIR not found: ${DATA_DIR}" >&2; exit 1; }
[[ -d "${OUT_DIR}" ]]  || { echo "ERROR: OUT_DIR not found: ${OUT_DIR}" >&2; exit 1; }
[[ -f "${METADATA}" ]] || { echo "ERROR: METADATA not found: ${METADATA}" >&2; exit 1; }
[[ -f "${PARAMS}" ]]   || { echo "ERROR: PARAMS not found: ${PARAMS}" >&2; exit 1; }
[[ -f "${INPUT_FILE}" ]] || { echo "ERROR: Expected input not found: ${INPUT_FILE}" >&2; exit 1; }

# IDS is optional: warn if set but missing, then ignore
if [[ -n "${IDS}" && ! -f "${IDS}" ]]; then
  echo "WARN: IDS mapping file not found at: ${IDS} (continuing without --ids)"
  unset IDS
fi

# ---- Safety: backup existing output if present ----
if [[ -f "${OUTPUT_FILE}" ]]; then
  ts="$(date +%Y%m%d_%H%M%S)"
  echo "NOTICE: ${OUTPUT_FILE} exists. Backing up to ${OUTPUT_FILE}.${ts}.bak"
  cp -p "${OUTPUT_FILE}" "${OUTPUT_FILE}.${ts}.bak"
fi

# ---- Show config ----
cat <<EOF
---------------------------------------------
 Batch Correction Test Runner
---------------------------------------------
PROJECT_ROOT : ${PROJECT_ROOT}
DATA_DIR     : ${DATA_DIR}
OUT_DIR      : ${OUT_DIR}
METADATA     : ${METADATA}
PARAMS       : ${PARAMS}
IDS          : ${IDS:-(none)}
COUNTER      : ${COUNTER}
LOG          : ${LOG}
MODULE       : ${MODULE}
RSCRIPT_BIN  : ${RSCRIPT_BIN}

Will read : ${INPUT_FILE}
Will write: ${OUTPUT_FILE}
---------------------------------------------
EOF

# ---- Run (pass --ids only if available) ----
RSARGS=( -d "${DATA_DIR}" -o "${OUT_DIR}" -p "${PARAMS}" -m "${METADATA}" -c "${COUNTER}" -l "${LOG}" )
if [[ -n "${IDS:-}" ]]; then
  RSARGS+=( -i "${IDS}" )
fi

"${RSCRIPT_BIN}" "${MODULE}" "${RSARGS[@]}"

echo "Done. Wrote: ${OUTPUT_FILE}"