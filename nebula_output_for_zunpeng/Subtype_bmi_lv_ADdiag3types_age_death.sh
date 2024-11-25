#!/bin/bash
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -p kellis
#SBATCH --output=slurm-%j-Subtype_bmi_lv_ADdiag3types_age_death.out
#SBATCH --job-name=bmi_
module load miniconda3/v4
source activate r_env
cd /home/anna_y/data/write/Subtype/
Rscript /home/anna_y/scRNA/nebula_test.R  bmi_lv /home/anna_y/data/write/Subtype/ /home/anna_y/scRNA/nebula_output_for_zunpeng/Subtype_bmi_lv_ADdiag3types_age_death/
echo ' is Done'
