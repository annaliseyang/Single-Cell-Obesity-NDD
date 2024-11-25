var1=Subtype # Subclass, Subtype
var2=bmi_lv # bmi_lv, BMI_norm
var3=ADdiag3types # ADdiag3types, None
covars=("msex" "pmi" "ADdiag3types" "total_counts" "nFeaturess_RNA" "age_death")
name=${var1}_${var2}_${var3}_age_death

# indir="/net/bmc-lab4/data/kellis/group/zunpeng/HumanBrainObsity/snRNA/PFC/All_$var1/"
indir="/home/anna_y/data/write/$var1/"
outdir="/home/anna_y/scRNA/nebula_output_for_zunpeng/$name/"
mkdir -p $outdir
# cd $indir
# Loop through .rds files and create job scripts for each cell type
for input_rds in $indir*.rds; do
    # Extract celltype name by removing specific prefixes and suffixes
    # cd $celltype_dir
    # input_rds=$celltype_dir/*.rds
    echo $input_rds
    # celltype_cleaned=$(echo "$celltype" | sed 's/rna4.AD427_only.//; s/.rds//')

    # # Display the cleaned celltype name
    # echo "$celltype_cleaned"

    # Define the output script path
    # script_path="$qs/$celltype_cleaned.sh"
    script_path="$outdir/$name.sh"

    # Write job script content to file
    {
        echo '#!/bin/bash'
        echo "#SBATCH -N 1"
        echo "#SBATCH -n 12"
        echo "#SBATCH -p kellis"
        echo "#SBATCH --output=slurm-%j-$name.out"
        # echo "#SBATCH --array=1-335%7"
        echo "#SBATCH --job-name=bmi_$name"
        # echo "source activate"
        # echo "conda activate R4.3.2"
        echo module load miniconda3/v4
        echo source activate r_env

        echo "cd $indir"
        echo "Rscript /home/anna_y/scRNA/nebula_output_for_zunpeng/nebula_v1.R $input_rds $var2 $indir $outdir"
        echo "echo '$input_rds is Done'"
    } > "$script_path"

    # Make the script executable
    chmod +x "$script_path"
    cd /home/anna_y/scRNA/
    sbatch $script_path
    break
done
