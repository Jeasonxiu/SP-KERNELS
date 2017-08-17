#!/bin/bash
for LABDEP in 120000 180000; do
for LABAMP in 5000 10000 20000; do 
for LABWL in 100000 200000 400000; do 
for SKIPSTA in 1 2 4 8 16; do
for SMOOTH in true false; do
sbatch job.sbatch $LABAMP $LABWL $LABDEP $SKIPSTA $SMOOTH
done
done
done
done
done
