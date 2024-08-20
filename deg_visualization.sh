

files="/home/anna_y/data/results/deg_bmi_lv/"
for file in $files**/*.Clean.tsv; do {
    # filename=$(basename "$file")
    echo "Processing $file..."
    python deg_visualization.py $file
}
done
