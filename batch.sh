#!/bin/bash
# Batch process all files

# pattern=/home/anna_y/data/write/AD_states/*/*.rds
# script=normalize_bmi.R

# pattern=/home/anna_y/data/write/AD_states/*/
# script=nebula.R

# pattern=/home/anna_y/data/results/deg_bmi_normalized/AD_states/*/*.Clean.tsv
# script=deg_results.py

pattern=/home/anna_y/data/write/Subclass/*/*.h5ad
# script=deg_heatmap.py
# script=deg_umap.py
script=deg_umap_gene.py

# Loop through each file in the directory
for file in $pattern; do
  echo ""
  echo "Running $script on $file"

  # python $script $file
  bash submit.sh $script $file
  # exit
done
