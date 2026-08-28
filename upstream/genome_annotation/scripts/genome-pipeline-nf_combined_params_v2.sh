#!/bin/bash
set -e
source runs/token.sh
source runs/shared_paths.sh

nextflow  run /home/janast/home/projects/genome_pipeline_nf/main.nf \
  -c /home/janast/home/projects/genome_pipeline_nf/configs/midway3.config \
  -params-file runs/genome_pipeline_nf/genome-pipeline-nf_combined_params_v2.yaml \
  -with-tower \
  -with-report logs/genome-pipeline-nf_combined_params_v2_report_$(date +%Y%m%d_%H%M%S).html \
  -resume