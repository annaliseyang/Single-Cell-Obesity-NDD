# Define directories

var1=Subtype # Subclass, Subtype
var2=BMI_norm # bmi_lv, BMI_norm
var3=ADdiag3types # ADdiag3types, None

indir="/net/bmc-lab4/data/kellis/group/zunpeng/HumanBrainObsity/snRNA/PFC/All_$var1/"
outdir="/home/anna_y/scRNA/nebula_output_for_zunpeng/${var1}_${var2}_${var3}/"
mkdir -p $outdir
cd $indir
# Loop through .rds files and create job scripts for each cell type
for celltype in *.rds; do
    # Extract celltype name by removing specific prefixes and suffixes
    celltype_cleaned=$(echo "$celltype" | sed 's/rna4.AD427_only.//; s/.rds//')

    # Display the cleaned celltype name
    echo "$celltype_cleaned"

    # Define the output script path
    # script_path="$qs/$celltype_cleaned.sh"
    script_path="/home/anna_y/data/scripts/nebula_test/$celltype_cleaned.sh"

    # Write job script content to file
    {
        echo '#!/bin/bash'
        echo "#SBATCH -N 1"
        echo "#SBATCH -n 12"
        echo "#SBATCH -p kellis"
        echo "#SBATCH --output=slurm-%j-${celltype_cleaned}_${var1}_${var2}_${var3}.out"
        # echo "#SBATCH --array=1-335%7"
        echo "#SBATCH --job-name=bmi_$celltype_cleaned"
        # echo "source activate"
        # echo "conda activate R4.3.2"
        echo module load miniconda3/v4
        echo source activate r_env

        echo "cd $indir"
        echo "Rscript /home/anna_y/scRNA/nebula_test.R $celltype_cleaned $var2 $indir $outdir"
        echo "echo '$celltype_cleaned is Done'"
    } > "$script_path"

    # Make the script executable
    chmod +x "$script_path"
    cd /home/anna_y/scRNA/
    sbatch $script_path
done
