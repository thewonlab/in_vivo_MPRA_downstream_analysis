#!/bin/bash
#SBATCH -t 08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --job-name=igvf_nbc_est
#SBATCH --export=ALL

set -euo pipefail

ROOT_DIR="/proj/hyejunglab/MPRA/meta_analysis/IGVF_InvivoMPRA_analysis_step"
SCRIPT="${ROOT_DIR}/scripts/step03b_estimate_variant_nbc.py"

BARCODE_MAP_VB="/proj/hyejunglab/MPRA/RNAseq/ADDICTION/Mint/In_Vivo_MPRA_analysis/Results/bcmapping/barcode_mapping_vb_v2.csv"
BARCODE_MAP_GVVC2="/proj/hyejunglab/MPRA/RNAseq/ADDICTION/Mint/In_Vivo_MPRA_analysis/Results/bcmapping/barcode_mapping_gvvc2_v2.csv"
VARIANT_MATRIX="/proj/hyejunglab/MPRA/RNAseq/ADDICTION/Mint/In_Vivo_MPRA_analysis/Results/count_matrix/count_matrix_variant_merged_v2_cortex.csv"
OUTPUT_CSV="${ROOT_DIR}/output/variant_n_bc_estimated.csv"

mkdir -p "${ROOT_DIR}/logs"

python "${SCRIPT}" \
  --barcode-mapping-vb "${BARCODE_MAP_VB}" \
  --barcode-mapping-gvvc2 "${BARCODE_MAP_GVVC2}" \
  --variant-matrix "${VARIANT_MATRIX}" \
  --output "${OUTPUT_CSV}"
