#!/bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --mem=8GB
# Default resources are 1 core with 2.8GB of memory.
# Specify a job name:
#SBATCH -J kernel_gen
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nicholas_mancinelli@brown.edu
#SBATCH -J kernel_gen
# Specify an output file
#SBATCH -o kernel_gen-%j.out
#SBATCH -e kernel_gen-%j.err
##
newdir=OUTPUT_FILES_$1
mkdir $newdir
cd $newdir
mkdir -p DATA
ln -s ../src
bash src/create_sem_synthetics.sh $2 $3 $4 $5
python src/estimate_S_arrival_time.py $5
python src/plot_seismos.py
cd ..
