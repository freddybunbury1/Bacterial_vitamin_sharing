#!/bin/bash
#SBATCH --job-name=genome_pipeline_nf
#SBATCH --account=pi-cdonnat
#SBATCH --partition=cdonnat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --output=logs/genome_pipeline_nf_%j.log
#SBATCH --error=logs/genome_pipeline_nf_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_ANNOTATION_DIR="$(dirname "$SCRIPT_DIR")"

MANIFEST="$GENOME_ANNOTATION_DIR/inputs/manifest.tsv"
PRIORITY="$GENOME_ANNOTATION_DIR/inputs/priority.tsv"
PARAMS_FILE="$SCRIPT_DIR/genome-pipeline-nf_combined_params_v2.yaml"
CONFIG="$GENOME_ANNOTATION_DIR/configs/midway3.config"
OUTDIR="${NF_OUTDIR:-$GENOME_ANNOTATION_DIR/output/nfcore}"

export NXF_WORK="${NXF_WORK:-/scratch/midway3/janast/nf_work/genome_pipeline_nf}"

mkdir -p "$OUTDIR" "$GENOME_ANNOTATION_DIR/logs" "$NXF_WORK"

module load nextflow/24.10.6

TOWER_ARGS=()
[[ -n "${TOWER_ACCESS_TOKEN:-}" ]] && TOWER_ARGS+=(-with-tower)

nextflow run Janastw/genome_pipeline_nf \
    -c "$CONFIG" \
    -params-file "$PARAMS_FILE" \
    --genome_manifest "$MANIFEST" \
    --priority_tsv "$PRIORITY" \
    --outdir "$OUTDIR" \
    "${TOWER_ARGS[@]}" \
    -with-report "$GENOME_ANNOTATION_DIR/logs/genome_pipeline_nf_report_$(date +%Y%m%d_%H%M%S).html" \
    -resume
