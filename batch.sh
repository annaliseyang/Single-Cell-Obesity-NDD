#!/bin/bash
# Batch process all files

# pattern=/home/anna_y/data/write/Class/*/*.rds # for normalize_bmi
# script=normalize_bmi.R

# pattern=/home/anna_y/data/write/Class/*/ # for nebula
# script=nebula.R

# pattern=/home/anna_y/results/deg_bmi_normalized/all/*/*.Clean.tsv
# script=deg_results.py

pattern=/home/anna_y/data/write/all/*/*.h5ad
# script=deg_heatmap.py
script=deg_umap.py

# Loop through each file in the directory
for file in $pattern; do
  echo ""
  echo "Running $script on $file"

  # python $script $file
  bash submit.sh $script $file
  # exit
done
