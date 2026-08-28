#!/bin/bash

input_dir="/project/cdonnat/microtrait/data/inputs/filtered_fna"
output_dir="/project/cdonnat/microtrait/data/outputs/gtdbtk/filtered_fna_nfcore_nb3u1/2026_05_29"

gtdbtk identify --genome_dir $input_dir --out_dir $output_dir/identify --extension fna --cpus 64
gtdbtk align --identify_dir $output_dir/identify --out_dir $output_dir/align --cpus 64
gtdbtk infer --msa_file $output_dir/align/align/gtdbtk.bac120.user_msa.fasta.gz --out_dir $output_dir/infer --cpus 64