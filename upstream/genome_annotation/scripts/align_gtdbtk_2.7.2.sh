#!/bin/bash
#SBATCH --job-name=gtdbtk_align
#SBATCH --account=pi-cdonnat
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gtdbtk_align_%j.log
#SBATCH --error=logs/gtdbtk_align_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
CPUS="${SLURM_CPUS_PER_TASK:-64}"

INPUT_DIR="${GTDBTK_INPUT_DIR:?Set GTDBTK_INPUT_DIR to the directory containing .fna assemblies}"
OUTPUT_DIR="${GTDBTK_OUTPUT_DIR:-$REPO_ROOT/upstream/genome_annotation/output/gtdbtk}"

mkdir -p "$OUTPUT_DIR" "$REPO_ROOT/upstream/genome_annotation/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gtdbtk_v2.7.2

gtdbtk identify \
    --genome_dir "$INPUT_DIR" \
    --out_dir "$OUTPUT_DIR/identify" \
    --extension fna \
    --cpus "$CPUS"

gtdbtk align \
    --identify_dir "$OUTPUT_DIR/identify" \
    --out_dir "$OUTPUT_DIR/align" \
    --cpus "$CPUS"

gtdbtk infer \
    --msa_file "$OUTPUT_DIR/align/align/gtdbtk.bac120.user_msa.fasta.gz" \
    --out_dir "$OUTPUT_DIR/infer" \
    --cpus "$CPUS"
