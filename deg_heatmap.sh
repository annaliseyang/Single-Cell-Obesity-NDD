#!/bin/bash

group=Class
pattern=/home/anna_y/data/write/$group/*/*.PFC.h5ad # input file
# results_file=/home/anna_y/data/write/$group/*/*.PFC.h5ad
script=deg_heatmap.py

# Loop through each file in the directory
for file in $pattern; do
  echo ""
  echo "Running $script on $file"

  results_file=$file | sed 's/\.h5ad$/\.bmi_norm.msex_pmi_total_counts_age_death.Clean.tsv/'
  results_file=$(echo "$file" | sed 's/\.h5ad$/\.bmi_norm.msex_pmi_total_counts_age_death.Clean.tsv/')
  echo "Results file: $results_file"

  # python $script $file
  bash submit.sh $script $file $results_file
#   exit
done
