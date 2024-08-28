# input_dir="/home/anna_y/data/write/Subclass/"
# pattern=/home/anna_y/data/results/deg_bmi_normalized/*/*/*_bmi.*
pattern=/home/anna_y/data/results/deg_bmi_normalized_v1_metascape/Class/*.zip

for file in $pattern; do
  echo ""
  echo "File: $file"
  # new_filename=${file/"_bmi."/"."}
  dir=${file/".zip"/""}
  echo Dir: $dir
  unzip $file -d $dir
  rm $file
  # echo New filename: $new_filename
  # mv $file $new_filename
done
