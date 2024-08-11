#!/bin/bash
# First argument is the script to be run, all other command line inputs are passed in

file=$1
shift 1
args=$@
echo Submitting $file...
echo Arguments: $args

if [[ "$file" == *".py" ]]; then
    sbatch job.sh sc2024 $file $args
elif [[ "$file" == *".R" ]]; then
    sbatch job.sh r_env $file $args
elif [[ "$file" == *".sh" ]]; then
    sbatch -p kellis -n 4 $file $args
else
    echo "Error submitting job $file"
fi
# sbatch job.sh $env $file $args
echo $file submitted.
