#!/bin/bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH --account=pi-cdonnat
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gtdbtk_classify_%j.log
#SBATCH --error=logs/gtdbtk_classify_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT=$(git -C "$SCRIPT_DIR" rev-parse --show-toplevel)
CPUS="${SLURM_CPUS_PER_TASK:-64}"

INPUT_DIR="${GTDBTK_INPUT_DIR:?Set GTDBTK_INPUT_DIR to the directory containing .fna assemblies}"
OUTPUT_DIR="${GTDBTK_OUTPUT_DIR:-$REPO_ROOT/upstream/genome_annotation/output/gtdbtk}"

mkdir -p "$OUTPUT_DIR" "$REPO_ROOT/upstream/genome_annotation/logs"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gtdbtk_v2.7.2

gtdbtk classify_wf \
    --genome_dir "$INPUT_DIR" \
    --out_dir "$OUTPUT_DIR/classify_wf" \
    --extension fna \
    --cpus "$CPUS"
