#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END
#SBATCH --mail-user=anna_y@mit.edu

# load environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate sc2024

script=$1
# parameters=range(0, 10)

# Loop through the list of scripts and submit each job
for $i in range(20); do
    echo "Submitting job for $script with parameter $i"
    sbatch -p kellis $script $i
done
