#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p kellis
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G

#SBATCH --mail-type=FAIL
#SBATCH --mail-user=anna_y@mit.edu


script_1="export_info.py"
input_file="$1"

script_2="create_rds.R"
input_dir=${input_file%".h5ad"}

echo "Running export_info_create_rds.sh"

# run jobs


# module load miniconda3/v4
# source activate sc2024
# python $script_1 $input_file

# module load miniconda3/v4
# source activate r_env
# Rscript $script_2 $input_dir

# bash job.sh sc2024 $script_1 $input_file
# bash job.sh r_env $script_2 $input_dir
