#!/bin/bash
#SBATCH -J BackProj
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH -e BackProj-%j.err
#SBATCH -o BackProj-%j.out

#echo "Starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME"
#matlab -r "MyMatlabFunction($SLURM_ARRAY_TASK_ID); quit;"
#
mode=1
lab_amplitude=$1
lab_wavlen=$2
lab_depth=$3
skipsta=$4
TakeDerivative=true
nderiv=1.5
ApplyGaussianSmoothingFilter=$5

label=LOOP-$lab_amplitude-$lab_wavlen-$lab_depth-$skipsta-$ApplyGaussianSmoothingFilter

matlab -r "master('$label', $mode, $lab_amplitude, $lab_wavlen, $lab_depth, $skipsta, $TakeDerivative, $nderiv, $ApplyGaussianSmoothingFilter); quit;"
