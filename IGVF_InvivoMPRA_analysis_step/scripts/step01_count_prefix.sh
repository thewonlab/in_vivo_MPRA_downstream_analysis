#!/usr/bin/env bash
set -euo pipefail

JSON_FILE=""
OUTPUT_DIR=""
PREFIX_LENGTH=20
JOBS=""

usage() {
    cat <<'EOF'
Usage:
  bash scripts/step01_count_prefix.sh \
    --json-file <fastq_dict.json> \
    --output-dir <count_dir> \
    [--prefix-length 20] \
    [--jobs 8]

Inputs:
  fastq_dict.json format:
    {
      "SAMPLE1": ["/path/to/sample1_L001_R1.fastq.gz", ...],
      "SAMPLE2": ["/path/to/sample2_L001_R1.fastq.gz", ...]
    }

Outputs:
  <output_dir>/<sample>_<lane>_count.txt
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --json-file)
            JSON_FILE="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --prefix-length)
            PREFIX_LENGTH="$2"
            shift 2
            ;;
        --jobs)
            JOBS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

if [[ -z "${JSON_FILE}" || -z "${OUTPUT_DIR}" ]]; then
    usage >&2
    exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
    echo "Missing dependency: jq" >&2
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "Missing dependency: GNU parallel" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

if [[ -z "${JOBS}" ]]; then
    if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
        JOBS="${SLURM_CPUS_PER_TASK}"
    else
        JOBS=1
    fi
fi

TASK_LIST=$(mktemp)
trap 'rm -f "${TASK_LIST}"' EXIT

jq -r 'to_entries[] | .key as $key | .value[] | "\($key)\t\(.)"' "${JSON_FILE}" > "${TASK_LIST}"

export OUTPUT_DIR PREFIX_LENGTH

parallel -j "${JOBS}" --colsep '\t' '
    sample={1}
    path={2}
    base=$(basename "$path")
    lane=$(echo "$base" | grep -o "L00[0-9]" || true)
    if [[ -z "$lane" ]]; then
        lane="L000"
    fi
    output_file="${OUTPUT_DIR}/${sample}_${lane}_count.txt"
    echo "Processing $path -> $output_file"
    zcat "$path" | awk -v prefix_len="${PREFIX_LENGTH}" "NR % 4 == 2 {print substr(\$0, 1, prefix_len)}" | sort | uniq -c | sort -nr > "$output_file"
'
