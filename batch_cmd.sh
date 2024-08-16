# input_dir="/home/anna_y/data/write/Subclass/"
pattern=/home/anna_y/data/write/*/*_small.h5ad

for file in $pattern; do
  echo ""
  echo "File: $file"
  new_filename=${file/"_small.h5ad"/"_50k.h5ad"}
  echo New filename: $new_filename
  mv $file $new_filename
done
