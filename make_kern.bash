#!/bin/bash
#SBATCH -J ConstructKernel
#SBATCH -t 2:00:00
#SBATCH -N 4
#SBATCH --mem=5G
#SBATCH -e ConstructKernel-%A-%a.err
#SBATCH -o ConstructKernel-%A-%a.out
#SBATCH --array=1,2,3
#
radius=$(echo $SLURM_ARRAY_TASK_ID | awk {'print $1*2500'})
matlab -r "addpath MATLAB-CODE/FUNCTIONS/; construct_kernel_matrix(0.01, $radius); exit"
