#!/bin/bash
#SBATCH -J BackProj
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH -e BackProj.err
#SBATCH -o BackProj.out

echo "Starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME"
matlab -r "MyMatlabFunction($SLURM_ARRAY_TASK_ID); quit;"
#
#
#
matlab -r "master('$1', 1, $2, $3, $4, $5)"
