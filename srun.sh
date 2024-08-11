srun -p kellis --pty bash

# load environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate sc2024

python init.py
