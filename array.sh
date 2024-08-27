#!/bin/bash
#SBATCH -p kellis
#SBATCH --array=1-472%7
#SBATCH --job-name=bmi_groups.py
#SBATCH --mem=55GB
#SBATCH --output=log/slurm-%A_%a-bmi_groups.py.out

SCRIPT=bmi_groups.py
INPUT_LIST=array_inputs.txt

# Load modules and activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate sc2024

# Get the input file based on the SLURM_ARRAY_TASK_ID
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $INPUT_LIST)

echo "Running $SCRIPT on $INPUT"

python $SCRIPT $INPUT

echo "Job Completed for $INPUT"
