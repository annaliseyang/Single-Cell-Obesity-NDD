#!/bin/bash
#SBATCH -p kellis
#SBATCH --array=1-335%7
#SBATCH --job-name=process_h5ad.py
#SBATCH --mem=0
#SBATCH --output=log/slurm-%A_%a-process_h5ad.py.out
#SBATCH --mail=FAIL

SCRIPT=process_h5ad.py
INPUT_LIST=inputs/process_h5ad.py_inputs.txt

# Load modules and activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate sc2024

# Get the input file based on the SLURM_ARRAY_TASK_ID
INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $INPUT_LIST)

echo "Running $SCRIPT on $INPUT"

python $SCRIPT $INPUT

echo "Job Completed for $INPUT"
