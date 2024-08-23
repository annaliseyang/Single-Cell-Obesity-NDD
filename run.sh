#!/bin/bash
# First argument is the script to be run, all other command line inputs are passed in

file=$1
shift 1
args=$@
echo Running $file...
echo Arguments: $args

# load modules
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit

if [[ "$file" == *".py" ]]; then
    # sbatch --job-name="$file" job.sh sc2024 $file $args
    conda activate sc2024
    python $file $args
elif [[ "$file" == *".R" ]]; then
    # sbatch --job-name="$file" job.sh r_env $file $args
    conda activate r_env
    Rscript $file $args
elif [[ "$file" == *".sh" ]]; then
    # sbatch -p kellis -n 4 --job-name="$file" $file $args
    bash $file $args
else
    echo "Error running $file"
fi

echo Job Completed!
