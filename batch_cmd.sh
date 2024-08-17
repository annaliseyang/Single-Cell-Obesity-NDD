# input_dir="/home/anna_y/data/write/Subclass/"
# pattern=/home/anna_y/data/write/*/*_small.h5ad
pattern=/home/anna_y/results/deg_bmi_normalized/*/*/*_bmi.*

for file in $pattern; do
  echo ""
  echo "File: $file"
  new_filename=${file/"_bmi."/"."}
  echo New filename: $new_filename
  mv $file $new_filename
done
