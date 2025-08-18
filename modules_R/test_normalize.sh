#!/usr/bin/env bash
set -Eeuo pipefail

# --- Variables ---
data_file="/mnt/c/MyStuff/preprocessing/Testing/Data"
param_path="/mnt/c/MyStuff/preprocessing/Testing/Data/parameters_1.json"
out_dir="/mnt/c/MyStuff/preprocessing/Testing/Data"
counter=4
log_file="/mnt/c/MyStuff/preprocessing/Testing/Data/log_file.txt"

# --- Call R script ---
"/usr/bin/Rscript" "/mnt/c/MyStuff/preprocessing/Testing/Code/normalize_R_testing_wrapper.R" \
    -d "$data_file" \
    -p "$param_path" \
    -o "$out_dir" \
    -c "$counter" \
    -l "$log_file"
