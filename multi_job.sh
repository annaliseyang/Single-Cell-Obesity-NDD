#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p kellis
#SBATCH --cpus-per-task=4
#SBATCH --mem=500G

#SBATCH --mail-type=END
#SBATCH --mail-user=anna_y@mit.edu

#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err

environment=$1
file=$2
shift 2
args=$@

# mkdir log/$SLURM_JOB_ID

echo Running: $file
echo Environment: $environment
echo Arguments: $args
echo Time of submission:
date

# load environment
module load miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate $environment

# run the file
echo Running $file...
start_time=$(date +%s)

if [[ "$environment" == "sc2024" ]]; then
    python $file $args
elif [[ "$environment" == "r_env" ]]; then
    Rscript $file $args
else
    # echo Error running $file
    echo "Environment not activated."
fi

# completion
echo Done! Time of completion:
date
end_time=$(date +%s)
time=$(( end_time - start_time ))

echo Run time: $(($time/3600)) hrs $(($time/60)) mins $(($time%60)) secs

echo "Job Completed: $SLURM_JOB_ID"
conda deactivate

# output_file="/home/anna_y/scRNA/out/slurm-$SLURM_JOB_ID.out"
# chmod 644 $output_file
# mail -s "Job Completed: $SLURM_JOB_ID" anna_y@mit.edu < $output_file
