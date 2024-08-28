#!/bin/bash

# Define paths
SCRIPT="process_h5ad.py"
INPUT_DIR="/home/anna_y/data/write/AD_states"
INPUT_LIST="inputs/${SCRIPT}_inputs.txt"
JOB_SCRIPT="array.sh"

# Find all .h5ad files and write to array_inputs.txt
find "$INPUT_DIR" -type f -name "*.h5ad" > "$INPUT_LIST"

# Calculate the number of lines (files)
NUM_LINES=$(wc -l < "$INPUT_LIST")

# Create the SLURM job script
cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#SBATCH -p kellis
#SBATCH --array=1-${NUM_LINES}%7
#SBATCH --job-name=${SCRIPT}
#SBATCH --mem=0
#SBATCH --output=log/slurm-%A_%a-${SCRIPT}.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=anna_y@mit.edu

SCRIPT=${SCRIPT}
INPUT_LIST=${INPUT_LIST}

# Load modules and activate conda environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate sc2024

# Get the input file based on the SLURM_ARRAY_TASK_ID
INPUT=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" \$INPUT_LIST)

echo "Running \$SCRIPT on \$INPUT"

python \$SCRIPT \$INPUT

echo "Job Completed for \$INPUT"
EOF

echo "SLURM job script generated: $JOB_SCRIPT"

sbatch $JOB_SCRIPT
echo "$NUM_LINES jobs queued."
