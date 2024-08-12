#!/bin/bash
# Batch process all files in a directory

# Directory containing the files to be processed
# input_dir="/home/anna_y/data/write/"
# input_dir="/home/anna_y/data/write/bmi_groups/"
# input_dir="/home/anna_y/data/write/AD_states/"
# input_dir="/home/anna_y/data/write/Class/"
# input_dir="/home/anna_y/data/write/Subclass/"
# input_dir="/home/anna_y/data/write/Subtype/"
# input_dir="/home/anna_y/data/write/msex/"
# input_dir="/home/anna_y/data/write/apoe_genotype/"

input_dir=$1
script=export_info_create_rds.sh

# Loop through each file in the directory
for file in "$input_dir"*.h5ad; do
  echo ""
  echo "Running $script on file: $file"
  # sbatch $script $file
  bash submit.sh $script $file
  # exit
done
