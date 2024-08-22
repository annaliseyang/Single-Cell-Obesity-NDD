#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p kellis
#SBATCH --cpus-per-task=4
#SBATCH --mem=0

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=anna_y@mit.edu

#SBATCH --output=log/slurm-%j-export_info_create_rds.sh.out
#SBATCH --error=log/slurm-%j-export_info_create_rds.sh.out


script_1="export_info.py"
input_file="$1" # .h5ad filepath

script_2="create_rds.R"
input_dir=${input_file%".h5ad"} # the directory created by running export_info.py

echo "Running export_info_create_rds.sh on $input_file"

# load modules
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

# first run export_info.py
conda activate sc2024
echo Running $script_1 on $input_file
python $script_1 $input_file
conda deactivate

# move the .h5ad file to its new directory
mv $input_file $input_dir

# then run create_rds.R on the new directory
conda activate r_env
echo Running $script_2 on $input_dir
Rscript $script_2 $input_dir
conda deactivate
