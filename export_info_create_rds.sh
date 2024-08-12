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
input_dir=${input_file%".h5ad"}

echo "Running export_info_create_rds.sh on $input_file"

# run jobs

module load miniconda3/v4
source activate sc2024
echo Running $script_1 on $input_file
python $script_1 $input_file

module load miniconda3/v4
source activate r_env
echo Running $script_2 on $input_dir
Rscript $script_2 $input_dir

# bash job.sh sc2024 $script_1 $input_file
# bash job.sh r_env $script_2 $input_dir
