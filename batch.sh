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

# pattern=$1
# pattern=/home/anna_y/data/write/*/*_50k.h5ad
pattern=/home/anna_y/data/write/Subtype/*/
# script=export_info_create_rds.sh
script=nebula.R

# Loop through each file in the directory
for file in $pattern; do
  echo ""
  echo "Running $script on $file"
  # sbatch $script $file
  bash submit.sh $script $file
  # exit
done
