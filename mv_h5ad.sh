#!/bin/sh

# input_dir="/home/anna_y/data/write/Subclass/"
pattern=/home/anna_y/data/write/Subtype/*.h5ad

for file in $pattern; do
  echo ""
  echo "File: $file"
  new_dir=${file/".h5ad"/"/"}

  echo New directory: $new_dir
  mv $file $new_dir
  # exit
done
