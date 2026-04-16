#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# This is a project-specific execution example that documents step order.
# It is not intended to be a fully self-contained reproduction wrapper.

# Step 1. FASTQ -> per-lane prefix counts
STEP01_SH="${ROOT_DIR}/scripts/step01_count_prefix.sh"
FASTQ_JSON="${ROOT_DIR}/config/fastq_dict.json"
COUNT_DIR="${ROOT_DIR}/output/counts"
PREFIX_LENGTH=20
JOBS=8

bash "${STEP01_SH}" \
  --json-file "${FASTQ_JSON}" \
  --output-dir "${COUNT_DIR}" \
  --prefix-length "${PREFIX_LENGTH}" \
  --jobs "${JOBS}"

# Step 2. Per-lane counts -> merged barcode count matrix
STEP02_PY="${ROOT_DIR}/scripts/step02_merge_counts.py"
MERGED_COUNT="${ROOT_DIR}/output/count_matrix_merged.csv"
python3 "${STEP02_PY}" \
  --count-dir "${COUNT_DIR}" \
  --fastq-dict "${FASTQ_JSON}" \
  --output "${MERGED_COUNT}"

# Step 3. Barcode -> variant annotation and variant-level matrices
STEP03_PY="${ROOT_DIR}/scripts/step03_variant_mapping.py"
BARCODE_MAP_VB="${ROOT_DIR}/config/barcode_mapping_vb.csv"
BARCODE_MAP_GVVC1="${ROOT_DIR}/config/barcode_mapping_gvvc1.csv"
BARCODE_MAP_GVVC2="${ROOT_DIR}/config/barcode_mapping_gvvc2.csv"
VARIANT_ANN="${ROOT_DIR}/output/count_matrix_merged_with_variant.csv"
VARIANT_FULL="${ROOT_DIR}/output/count_matrix_variant_full.csv"
VARIANT_FILTERED="${ROOT_DIR}/output/count_matrix_variant_filtered.csv"
python3 "${STEP03_PY}" \
  --count-matrix-merged "${MERGED_COUNT}" \
  --barcode-mapping-vb "${BARCODE_MAP_VB}" \
  --barcode-mapping-gvvc1 "${BARCODE_MAP_GVVC1}" \
  --barcode-mapping-gvvc2 "${BARCODE_MAP_GVVC2}" \
  --output-with-variant "${VARIANT_ANN}" \
  --output-variant-full "${VARIANT_FULL}" \
  --output-variant-filtered "${VARIANT_FILTERED}"

# Step 4. Replicate pooling
STEP04_PY="${ROOT_DIR}/scripts/step04_pooling.py"
POOLING_SCHEME="${ROOT_DIR}/pooling_scheme.csv"
POOL_OUT="${ROOT_DIR}/output/pooling"
POOL_PREFIX="pooling_scheme_3"
python3 "${STEP04_PY}" \
  --input "${VARIANT_ANN}" \
  --output "${POOL_OUT}" \
  --prefix "${POOL_PREFIX}" \
  --pooling_scheme "${POOLING_SCHEME}"

# Step 5. MPRA statistics
STEP05_R="${ROOT_DIR}/scripts/step05_mpra_analysis.R"
CONFIG_JSON="${ROOT_DIR}/config/analysis_config.json"
MPRA_OUT="${ROOT_DIR}/output/${POOL_PREFIX}"
Rscript "${STEP05_R}" \
  --input_file "${POOL_OUT}" \
  --prefix "${POOL_PREFIX}" \
  --config_json "${CONFIG_JSON}" \
  --output_dir "${MPRA_OUT}" \
  --num_cores 8

# Step 6. Format final IGVF outputs
STEP06_PY="${ROOT_DIR}/scripts/step06_format_igvf_outputs.py"
STEP06_CONFIG="${ROOT_DIR}/config/step6_paths.template.json"
python3 "${STEP06_PY}" --config-json "${STEP06_CONFIG}"

echo "[완료] 전체 파이프라인이 끝났습니다."
