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

# pattern=/home/anna_y/data/write/Class/*/*.rds # for normalize_bmi
# script=normalize_bmi.R

# pattern=/home/anna_y/data/write/Class/*/ # for nebula
# script=nebula.R

pattern=/home/anna_y/results/deg_bmi_normalized/*/*/*.Clean.tsv
script=deg_results.py

# Loop through each file in the directory
for file in $pattern; do
  echo ""
  echo "Running $script on $file"

  python $script $file
  # bash submit.sh $script $file
  # exit
done
