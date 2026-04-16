# IGVF In Vivo MPRA Analysis Step

This directory documents the downstream analysis workflow used for the in vivo
MPRA IGVF submission.

It should be read primarily as a project record of how the analysis was
performed, not as a fully portable or fully reproducible software package.
Several scripts retain project-specific assumptions about sample naming,
directory structure, metadata columns, and available reference files.

## What this directory is for

- record the sequence of downstream analysis steps used in this project
- preserve cleaned copies of the main scripts used for those steps
- show how the final IGVF submission tables were derived

## What this directory is not

- a general-purpose MPRA framework
- a guaranteed end-to-end reproduction package for another environment
- a complete archive of every upstream input needed to rerun the project from scratch

## Analysis outline

The documented workflow is organized into six steps.

1. `scripts/step01_count_prefix.sh`
   - counts 5' prefix or barcode sequences from FASTQ files
   - adapted from `Count_Matrix_Kat_26_02_03/Module/A1.word_count.sh`
2. `scripts/step02_merge_counts.py`
   - merges lane-level count tables into a sample-level barcode count matrix
   - adapted from `Count_Matrix_Kat_26_02_03/Module/A2.build_count_matrix.py`
3. `scripts/step03_variant_mapping.py`
   - maps barcode rows to variant IDs and builds variant-level count matrices
   - adapted from `In_Vivo_MPRA_analysis/Module/A1_0.ct_mat.py`
   - related sidecar utility: `scripts/step03b_estimate_variant_nbc.py`
4. `scripts/step04_pooling.py`
   - pools replicate-level RNA and DNA counts according to the study pooling scheme
   - adapted from `In_Vivo_MPRA_analysis/Module/C01.pooling_replicates.py`
5. `scripts/step05_mpra_analysis.R`
   - runs activity and allelic MPRA analyses on pooled counts
   - adapted from `In_Vivo_MPRA_analysis/Module/C02.mpra_analysis.R`
6. `scripts/step06_format_igvf_outputs.py`
   - formats the project results into IGVF submission tables

## Directory contents

- `scripts/`: cleaned step-level scripts reflecting the project workflow
- `config/`: config templates used by the formatter or example wrapper
- `pooling_scheme.csv`: local copy of the pooling assignment table
- `output/`: example final IGVF output files
- `run_all_steps.sh`: project-specific example wrapper showing step order
- `input.md`: notes on upstream files used by the analysis
- `scripts/run_step03b_estimate_variant_nbc.sbatch`: separate Slurm entrypoint for optional `n_bc` estimation

## Final IGVF outputs

This workflow ultimately feeds five submission files:

1. `Reporter_Experiment.tsv.gz`
2. `Reporter_Element.tsv.gz`
3. `Reporter_Variant.tsv.gz`
4. `Reporter_Genomic_Element_Effect.bed.gz`
5. `Reporter_Genomic_Variant_Effect.bed.gz`

Example output files are included under `output/` as references.

## Script behavior and limitations

- The scripts in this directory were cleaned to make the flow easier to follow,
  but they still mirror project-specific conventions from the original analysis.
- Sample naming assumptions such as `_CORT_R` and `_CORT_D` are preserved where
  they were part of the original workflow.
- Some runtime inputs referenced by the scripts are not included in this
  directory.
- `run_all_steps.sh` is best understood as an execution example for this
  project, not as a claim that the directory is self-contained.
- `n_bc` estimation is treated as a separate sidecar calculation rather than a
  required part of the display pipeline.
- The IGVF formatter depends on finalized input paths and field mappings that
  come from this specific submission context.

## Step examples

The commands below illustrate how each step was invoked in this project style.
They are examples, not a portability guarantee.

### Step 1. Prefix counting

```bash
bash scripts/step01_count_prefix.sh \
  --json-file config/fastq_dict.json \
  --output-dir output/counts \
  --prefix-length 20 \
  --jobs 8
```

Requirements:

- `jq`
- `parallel`

### Step 2. Merge count tables

```bash
python scripts/step02_merge_counts.py \
  --count-dir output/counts \
  --fastq-dict config/fastq_dict.json \
  --output output/count_matrix_merged.csv
```

Optional flags:

- `--non-strict-sample-match`
- `--log-level DEBUG|INFO|WARNING|ERROR`

### Step 3. Variant mapping

```bash
python scripts/step03_variant_mapping.py \
  --count-matrix-merged output/count_matrix_merged.csv \
  --barcode-mapping-vb config/barcode_mapping_vb.csv \
  --barcode-mapping-gvvc1 config/barcode_mapping_gvvc1.csv \
  --barcode-mapping-gvvc2 config/barcode_mapping_gvvc2.csv \
  --output-with-variant output/count_matrix_merged_with_variant.csv \
  --output-variant-full output/count_matrix_variant_full.csv \
  --output-variant-filtered output/count_matrix_variant_filtered.csv
```

### Step 4. Replicate pooling

```bash
python scripts/step04_pooling.py \
  --input output/count_matrix_merged_with_variant.csv \
  --output output/pooling \
  --prefix pooling_scheme_3 \
  --pooling_scheme pooling_scheme.csv
```

### Step 5. MPRA analysis

```bash
Rscript scripts/step05_mpra_analysis.R \
  --input_file output/pooling \
  --prefix pooling_scheme_3 \
  --config_json config/analysis_config.json \
  --output_dir output/pooling_scheme_3 \
  --num_cores 8
```

### Step 6. IGVF formatting

```bash
python scripts/step06_format_igvf_outputs.py \
  --config-json config/step6_paths.template.json
```

## Example wrapper

```bash
bash run_all_steps.sh
```

This wrapper shows the intended order of operations for this project. It assumes
that the required project-specific config files and upstream resources already
exist.

## Python dependencies

```bash
python -m pip install -r requirements.txt
```
